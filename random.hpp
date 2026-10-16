#ifndef RANDOM_H_
#define RANDOM_H_

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include "typedefs.hpp"

class RandomContext {
  std::mt19937_64 mEngine;

public:
  explicit RandomContext(std::uint64_t seed) : mEngine(seed) {}

  Real Uniform() {
    return std::generate_canonical<Real, std::numeric_limits<Real>::digits>(mEngine);
  }

  unsigned long Integer(unsigned long upper_exclusive) {
    if (upper_exclusive == 0) return 0;
    std::uniform_int_distribution<unsigned long> dist(0, upper_exclusive - 1);
    return dist(mEngine);
  }
};

inline RandomContext *&ActiveRandomContext() {
  static thread_local RandomContext *context = 0;
  return context;
}

inline void SetActiveRandomContext(RandomContext *context) {
  ActiveRandomContext() = context;
}

inline Real RandomUniform() {
  RandomContext *context = ActiveRandomContext();
  static thread_local std::mt19937_64 diagnostic_engine(0x7261643364ULL);
  if (!context) return std::generate_canonical<Real, std::numeric_limits<Real>::digits>(diagnostic_engine);
  return context->Uniform();
}

inline unsigned long RandomInteger(unsigned long upper_exclusive) {
  RandomContext *context = ActiveRandomContext();
  static thread_local std::mt19937_64 diagnostic_engine(0x7261643364ULL);
  if (!context) {
    if (upper_exclusive == 0) return 0;
    std::uniform_int_distribution<unsigned long> dist(0, upper_exclusive - 1);
    return dist(diagnostic_engine);
  }
  return context->Integer(upper_exclusive);
}

#endif
