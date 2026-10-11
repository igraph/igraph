/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef RANDOM_H_
#define RANDOM_H_

#include "infomath.h"

#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace infomap {
using RandGen = std::mt19937;
using StdUniformUIntDistribution = std::uniform_int_distribution<unsigned int>;

class PortableUniformUIntDistribution {
public:
  struct param_type {
    param_type(unsigned int minValue, unsigned int maxValue) : min(minValue), max(maxValue) {}

    unsigned int min;
    unsigned int max;
  };

  // Portable unbiased bounded integer mapping, independent of std::uniform_int_distribution internals.
  unsigned int operator()(RandGen& randGen, const param_type& params)
  {
    const auto min = params.min;
    const auto max = params.max;
    const std::uint64_t range = static_cast<std::uint64_t>(max) - min + 1u;
    const std::uint64_t generatorRange = static_cast<std::uint64_t>(RandGen::max()) - RandGen::min() + 1u;
    const std::uint64_t limit = generatorRange - (generatorRange % range);

    std::uint64_t value = 0;
    do {
      value = static_cast<std::uint64_t>(randGen()) - RandGen::min();
    } while (value >= limit);

    return min + static_cast<unsigned int>(value % range);
  }
};

namespace detail {
  // Detects Engine::randInt(min, max) — an engine wrapping a host RNG with its
  // own bounded-draw API (e.g. igraph_rng_get_integer). Such engines supply
  // bounded integers directly, keeping the host's platform-stable mapping.
  template <typename Engine, typename = void>
  struct HasBoundedDraw : std::false_type {};

  template <typename Engine>
  struct HasBoundedDraw<Engine, decltype(void(std::declval<Engine&>().randInt(0u, 0u)))> : std::true_type {};

  // Detects a UniformRandomBitGenerator-style callable, for a readable
  // setEngine() diagnostic when an engine matches neither shape.
  template <typename Engine, typename = void>
  struct IsUniformBitGenerator : std::false_type {};

  template <typename Engine>
  struct IsUniformBitGenerator<Engine, decltype(void(std::declval<Engine&>()()))> : std::true_type {};
} // namespace detail

template <typename UniformUIntDistribution>
class BasicRandom {
  // Default fast path: the built-in mt19937 + the templated distribution.
  RandGen m_randGen;
  UniformUIntDistribution m_uniform;

  // Optional injected engine for native embedders (e.g. igraph), issue #411.
  // Null keeps the default fast path, so default behaviour and results are
  // unchanged. An engine with its own bounded draw (randInt) is used as-is; a
  // plain bit generator is mapped with the standard distribution. Either way
  // the model is independent of the template distribution and is only
  // instantiated when setEngine() is actually used.
  struct EngineConcept {
    virtual ~EngineConcept() = default;
    virtual std::unique_ptr<EngineConcept> clone() const = 0;
    virtual void seed(unsigned int seedValue) = 0;
    virtual unsigned int randInt(unsigned int min, unsigned int max) = 0;
  };

  template <typename Engine>
  struct EngineModel : EngineConcept {
    Engine engine;
    StdUniformUIntDistribution uniform;

    explicit EngineModel(Engine engineValue) : engine(std::move(engineValue)) {}

    std::unique_ptr<EngineConcept> clone() const override
    {
      return std::unique_ptr<EngineConcept>(new EngineModel(*this));
    }

    void seed(unsigned int seedValue) override { engine.seed(seedValue); }

    unsigned int randInt(unsigned int min, unsigned int max) override
    {
      return randIntImpl(min, max, detail::HasBoundedDraw<Engine> {});
    }

    // The engine's own bounded draw wins over bit generation: a host mapping
    // (e.g. igraph_rng_get_integer) is platform-stable, where the standard
    // distribution's bit-to-integer mapping is implementation-defined.
    unsigned int randIntImpl(unsigned int min, unsigned int max, std::true_type)
    {
      return engine.randInt(min, max);
    }

    unsigned int randIntImpl(unsigned int min, unsigned int max, std::false_type)
    {
      return uniform(engine, StdUniformUIntDistribution::param_type(min, max));
    }
  };

  std::unique_ptr<EngineConcept> m_customEngine;

public:
  BasicRandom(unsigned int seed = 123) : m_randGen(seed) {}

  // Copyable so InfomapConfig (which holds a Random by value) stays copyable and
  // the engine propagates to sub-Infomaps; clone the type-erased engine.
  BasicRandom(const BasicRandom& other)
      : m_randGen(other.m_randGen),
        m_uniform(other.m_uniform),
        m_customEngine(other.m_customEngine ? other.m_customEngine->clone() : nullptr)
  {
  }

  BasicRandom& operator=(const BasicRandom& other)
  {
    if (this != &other) {
      m_randGen = other.m_randGen;
      m_uniform = other.m_uniform;
      m_customEngine = other.m_customEngine ? other.m_customEngine->clone() : nullptr;
    }
    return *this;
  }

  BasicRandom(BasicRandom&&) noexcept = default;
  BasicRandom& operator=(BasicRandom&&) noexcept = default;

  // Install a seedable RNG engine. Two shapes are accepted: an engine with its
  // own bounded draw (unsigned int randInt(min, max), used as-is so a host
  // RNG's integer mapping stays platform-stable), or a standard
  // UniformRandomBitGenerator whose bits are mapped with
  // std::uniform_int_distribution (implementation-defined sequence). randInt()
  // wins when both are present. seed(unsigned int) is required either way and
  // may be a no-op when the host owns seeding. Must be copy-constructible
  // because Random is copied to sub-Infomaps. Native-only; the language
  // bindings do not expose it (#411).
  template <typename Engine>
  void setEngine(Engine engine)
  {
    using EngineType = typename std::decay<Engine>::type;
    static_assert(std::is_copy_constructible<EngineType>::value,
                  "Random::setEngine() requires a copy-constructible engine.");
    static_assert(detail::HasBoundedDraw<EngineType>::value || detail::IsUniformBitGenerator<EngineType>::value,
                  "Random::setEngine() requires an engine with either randInt(min, max) "
                  "or a UniformRandomBitGenerator operator()().");
    m_customEngine.reset(new EngineModel<EngineType>(std::move(engine)));
  }

  void seed(unsigned int seedValue)
  {
    if (m_customEngine) {
      m_customEngine->seed(seedValue);
      return;
    }
    m_randGen.seed(seedValue);
  }

  unsigned int randInt(unsigned int min, unsigned int max)
  {
    // Hinted so the default path stays the fall-through: no engine installed
    // is the common case, and randInt() runs once per optimizer draw.
    if (INFOMAP_UNLIKELY(m_customEngine != nullptr))
      return m_customEngine->randInt(min, max);
    return m_uniform(m_randGen, typename UniformUIntDistribution::param_type(min, max));
  }

  /**
   * Get a random permutation of indices of the size of the input vector
   */
  void getRandomizedIndexVector(std::vector<unsigned int>& randomOrder)
  {
    unsigned int size = randomOrder.size();
    std::iota(randomOrder.begin(), randomOrder.end(), 0u);
    for (unsigned int i = 0; i < size; ++i)
      std::swap(randomOrder[i], randomOrder[i + randInt(0, size - i - 1)]);
  }
};

using StdRandom = BasicRandom<StdUniformUIntDistribution>;
using PortableRandom = BasicRandom<PortableUniformUIntDistribution>;
using Random = StdRandom;
} // namespace infomap

#endif // RANDOM_H_
