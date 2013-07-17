#pragma once

// Copyright (c) 2012, Sebastian Jeltsch (sjeltsch@kip.uni-heidelberg.de)
// Distributed under the terms of the GPLv2 or newer

#include <chrono>
#include <string>

#include <iostream>
#include <iomanip>
#include <map>

#define BENCH_THRESH    std::chrono::seconds(1)
#define BENCH_DURATION  std::chrono::nanoseconds
#define BENCH_CLOCK     std::chrono::steady_clock
#define BENCH_SPACING   12

#define BENCH_QUOTE(STR) #STR
#define BENCH_TO_STRING(STR) BENCH_QUOTE(STR)

namespace bench {

using namespace std::chrono;

typedef time_point<BENCH_CLOCK> time_t;

namespace detail {

inline
double reference_time(
	std::string const& name,
	double const measured_frequency)
{
	typedef std::map<std::string, double> type;
	static type _m;

	type::iterator it;
	if ((it = _m.find(name)) != _m.end())
		return it->second;

	_m[name] = measured_frequency;
	return measured_frequency;
}

inline
double duration_to_frequency(
	size_t const iterations,
	BENCH_DURATION const duration)
{
	return double(iterations) / duration_cast<nanoseconds>(duration).count() *
		nanoseconds::period::den;
}

template<typename Time>
inline
BENCH_DURATION delta_t(Time const t)
{
	return duration_cast<BENCH_DURATION>(BENCH_CLOCK::now() - t);
}

inline
void message(
	std::string const& name,
	size_t const iterations,
	BENCH_DURATION const duration)
{
	double frequency = duration_to_frequency(iterations, duration);
	double reference = reference_time(name, frequency);

	using namespace std;
	cout << setw(BENCH_SPACING) << left << name << right;
	cout.precision(2);
	cout << setw(BENCH_SPACING) << frequency   << " it/s    "
	     << setw(BENCH_SPACING) << 1/frequency << " s/it    "
	     << setw(BENCH_SPACING) << fixed << reference/frequency * 100 << "%"
		 << std::endl;
}

} // namespace detail

template <class T>
inline
void preserve(T&& val)
{
	  asm volatile("" : "+r" (val));
}

template<typename Lambda>
void run(
	std::string const& name,
	Lambda const& code,
	size_t const iterations = 1,
	size_t const offset = 0,
	time_t const time = BENCH_CLOCK::now())
{
	for (size_t ii = offset; ii < iterations ; ++ii)
		code();

	detail::message(name, iterations, detail::delta_t(time));
}

} // namespace bench

#define BENCH(NAME, ...) \
{ \
	auto __lambda = [&](){ __VA_ARGS__ }; \
	auto __t = BENCH_CLOCK::now(); \
	__lambda(); \
	auto __d = bench::detail::delta_t(__t); \
	size_t __end = BENCH_THRESH / \
		(__d.count() ? __d : std::chrono::nanoseconds(1)); \
	bench::run(BENCH_TO_STRING(NAME), __lambda, __end, 1, __t); \
}

#define BENCH_N(NAME, N, ...) \
{ \
	auto __lambda = [&](){ __VA_ARGS__ }; \
	bench::run(BENCH_TO_STRING(NAME), __lambda, N); \
}

#undef BENCH_DURATION
#undef BENCH_SPACING
