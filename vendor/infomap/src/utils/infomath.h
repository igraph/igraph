/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMATH_H_
#define INFOMATH_H_

#include <cmath>
#include <cstdlib>

#if defined(__GNUC__) || defined(__clang__)
#define INFOMAP_HOT __attribute__((hot))
#define INFOMAP_ALWAYS_INLINE __attribute__((always_inline))
#define INFOMAP_LIKELY(x) __builtin_expect(!!(x), 1)
#define INFOMAP_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define INFOMAP_HOT
#define INFOMAP_ALWAYS_INLINE
#define INFOMAP_LIKELY(x) (x)
#define INFOMAP_UNLIKELY(x) (x)
#endif

#if defined(INFOMAP_USE_SIMD_LOG) && defined(__ARM_NEON) && defined(__aarch64__)
#include <arm_neon.h>
#define INFOMAP_NEON_LOG 1
#endif

#if defined(INFOMAP_USE_SIMD_LOG) && defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>
#define INFOMAP_AVX2_LOG 1
#endif

namespace infomap {
namespace infomath {

  using std::log2;

  INFOMAP_HOT INFOMAP_ALWAYS_INLINE inline double plogp(double p)
  {
    return INFOMAP_LIKELY(p > 0.0) ? p * log2(p) : 0.0;
  }

#if defined(INFOMAP_AVX2_LOG)
  // Inlined SIMD log2 for four doubles via AVX2 + FMA. Same atanh-style series and
  // accuracy bound as the Neon implementation below.
  INFOMAP_HOT INFOMAP_ALWAYS_INLINE inline __m256d avx2_log2_pd(__m256d x)
  {
    __m256i bits = _mm256_castpd_si256(x);
    __m256i exp_raw = _mm256_srli_epi64(bits, 52);
    exp_raw = _mm256_and_si256(exp_raw, _mm256_set1_epi64x(0x7FFLL));
    __m256i exp_i64 = _mm256_sub_epi64(exp_raw, _mm256_set1_epi64x(1023LL));

    __m256i mant_bits = _mm256_and_si256(bits, _mm256_set1_epi64x(0x000FFFFFFFFFFFFFLL));
    mant_bits = _mm256_or_si256(mant_bits, _mm256_set1_epi64x(0x3FF0000000000000LL));
    __m256d m = _mm256_castsi256_pd(mant_bits);

    const __m256d SQRT2 = _mm256_set1_pd(1.4142135623730951);
    __m256d big = _mm256_cmp_pd(m, SQRT2, _CMP_GT_OQ);
    __m256d m_half = _mm256_mul_pd(m, _mm256_set1_pd(0.5));
    m = _mm256_blendv_pd(m, m_half, big);
    __m256i big_i = _mm256_castpd_si256(big);
    __m256i inc = _mm256_and_si256(big_i, _mm256_set1_epi64x(1LL));
    exp_i64 = _mm256_add_epi64(exp_i64, inc);

    // Exponent always fits in int32 (|e| <= 1023). Pack low 32 bits of each int64
    // lane into a 128-bit int32x4, then convert to four doubles.
    __m256i shuffled = _mm256_shuffle_epi32(exp_i64, _MM_SHUFFLE(2, 0, 2, 0));
    __m256i permuted = _mm256_permute4x64_epi64(shuffled, _MM_SHUFFLE(3, 1, 2, 0));
    __m128i exp_i32 = _mm256_castsi256_si128(permuted);
    __m256d ef = _mm256_cvtepi32_pd(exp_i32);

    __m256d one = _mm256_set1_pd(1.0);
    __m256d s = _mm256_div_pd(_mm256_sub_pd(m, one), _mm256_add_pd(m, one));
    __m256d z = _mm256_mul_pd(s, s);

    __m256d P = _mm256_set1_pd(1.0 / 19.0);
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 17.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 15.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 13.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 11.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 9.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 7.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 5.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0 / 3.0));
    P = _mm256_fmadd_pd(P, z, _mm256_set1_pd(1.0));

    __m256d log_m = _mm256_mul_pd(_mm256_mul_pd(s, P), _mm256_set1_pd(2.0));
    const __m256d INV_LN2 = _mm256_set1_pd(1.4426950408889634);
    return _mm256_fmadd_pd(log_m, INV_LN2, ef);
  }
#endif

#if defined(INFOMAP_NEON_LOG)
  // Inlined SIMD log2 for two doubles via arm64 Neon. Uses the atanh-style series
  //   log(m) = 2 * s * (1 + z/3 + z^2/5 + ... + z^9/19),   s = (m-1)/(m+1), z = s^2
  // after range reduction at sqrt(2) so |s| < 0.172. The 10-term polynomial gives
  // ~1 ulp accuracy on the reduced interval. Input must be strictly positive and
  // finite; callers are expected to mask non-positive inputs before invoking.
  INFOMAP_HOT INFOMAP_ALWAYS_INLINE inline float64x2_t neon_log2_pd(float64x2_t x)
  {
    uint64x2_t bits = vreinterpretq_u64_f64(x);
    int64x2_t exp_raw = vreinterpretq_s64_u64(vshrq_n_u64(bits, 52));
    int64x2_t exp = vsubq_s64(exp_raw, vdupq_n_s64(1023));

    uint64x2_t mant_bits = vandq_u64(bits, vdupq_n_u64(0x000FFFFFFFFFFFFFULL));
    mant_bits = vorrq_u64(mant_bits, vdupq_n_u64(0x3FF0000000000000ULL));
    float64x2_t m = vreinterpretq_f64_u64(mant_bits);

    const float64x2_t SQRT2 = vdupq_n_f64(1.4142135623730951);
    uint64x2_t big = vcgtq_f64(m, SQRT2);
    m = vbslq_f64(big, vmulq_f64(m, vdupq_n_f64(0.5)), m);
    exp = vbslq_s64(big, vaddq_s64(exp, vdupq_n_s64(1)), exp);
    float64x2_t ef = vcvtq_f64_s64(exp);

    float64x2_t one = vdupq_n_f64(1.0);
    float64x2_t s = vdivq_f64(vsubq_f64(m, one), vaddq_f64(m, one));
    float64x2_t z = vmulq_f64(s, s);

    float64x2_t P = vdupq_n_f64(1.0 / 19.0);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 17.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 15.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 13.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 11.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 9.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 7.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 5.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0 / 3.0), P, z);
    P = vfmaq_f64(vdupq_n_f64(1.0), P, z);

    float64x2_t log_m = vmulq_f64(vmulq_f64(s, P), vdupq_n_f64(2.0));
    const float64x2_t INV_LN2 = vdupq_n_f64(1.4426950408889634);
    return vfmaq_f64(ef, log_m, INV_LN2);
  }
#endif

  // Batched p * log2(p) for n independent inputs. The scalar path is bit-exact with
  // calling plogp() in a loop; SIMD paths use an inlined polynomial for log2 (~1 ulp
  // deviation from std::log2 on finite positive inputs). Inputs and outputs must
  // not alias. Non-positive inputs are masked to a zero output. Non-finite positive
  // inputs (+Inf, NaN) produce undefined results on the SIMD paths — Infomap flow
  // values are always finite, so this isn't a runtime concern, but callers passing
  // synthetic data should pre-validate.
  INFOMAP_HOT inline void plogp_batch(const double* in, double* out, int n)
  {
#if defined(INFOMAP_AVX2_LOG)
    const __m256d zero4 = _mm256_setzero_pd();
    const __m256d one4 = _mm256_set1_pd(1.0);
    int i = 0;
    for (; i + 4 <= n; i += 4) {
      __m256d x = _mm256_loadu_pd(in + i);
      __m256d pos = _mm256_cmp_pd(x, zero4, _CMP_GT_OQ);
      __m256d safe = _mm256_blendv_pd(one4, x, pos);
      __m256d lg = avx2_log2_pd(safe);
      __m256d product = _mm256_mul_pd(x, lg);
      __m256d result = _mm256_blendv_pd(zero4, product, pos);
      _mm256_storeu_pd(out + i, result);
    }
    for (; i < n; ++i) {
      out[i] = INFOMAP_LIKELY(in[i] > 0.0) ? in[i] * log2(in[i]) : 0.0;
    }
#elif defined(INFOMAP_NEON_LOG)
    const float64x2_t zero = vdupq_n_f64(0.0);
    const float64x2_t one = vdupq_n_f64(1.0);
    int i = 0;
    for (; i + 2 <= n; i += 2) {
      float64x2_t x = vld1q_f64(in + i);
      uint64x2_t pos = vcgtq_f64(x, zero);
      // For non-positive lanes, feed 1.0 to log2 (returns 0) and zero the product below.
      float64x2_t safe = vbslq_f64(pos, x, one);
      float64x2_t lg = neon_log2_pd(safe);
      float64x2_t product = vmulq_f64(x, lg);
      float64x2_t result = vbslq_f64(pos, product, zero);
      vst1q_f64(out + i, result);
    }
    for (; i < n; ++i) {
      out[i] = INFOMAP_LIKELY(in[i] > 0.0) ? in[i] * log2(in[i]) : 0.0;
    }
#else
    for (int i = 0; i < n; ++i) {
      out[i] = INFOMAP_LIKELY(in[i] > 0.0) ? in[i] * log2(in[i]) : 0.0;
    }
#endif
  }

  inline double isEqual(double a, double b, double tol = 1e-8)
  {
    return std::abs(a - b) <= tol;
  }

  /**
   * Tsallis entropy S_q of a uniform probability distribution of length n
   */
  inline double tsallisEntropyUniform(double n, double q = 1)
  {
    if (isEqual(q, 1)) {
      return std::log2(n);
    }
    return 1 / (q - 1) * (1 - pow(n, (1 - q))) / std::log(2);
  }

  /**
   * Interpolate from linear (q = 0) to log (q = 1)
   * linlog(k, 0) = k
   * linlog(k, 1) = log2(k)
   */
  inline double linlog(double k, double q = 1)
  {
    double baseCorrection = q <= 1 ? (1 - q) * std::log(2) + q : 1;
    double offsetCorrection = q <= 1 ? 1 - q : 0;
    return tsallisEntropyUniform(k, q) * baseCorrection + offsetCorrection;
  }

} // namespace infomath
} // namespace infomap

#endif // INFOMATH_H_
