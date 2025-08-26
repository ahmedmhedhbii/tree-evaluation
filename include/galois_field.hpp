/**
 * Galois Field Implementation for Tree Evaluation Algorithm
 * This file implements finite field arithmetic over GF(2^n) as required by the
 * "Tree Evaluation is in Space O(log n log log n)" paper.
 */

#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <stdexcept>

inline int leading_zeros(uint64_t x) {
#if defined(__cpp_lib_bitops) && __cpp_lib_bitops >= 201907L
  return std::countl_zero(x);
#elif defined(__GNUG__) || defined(__clang__)
  return x ? __builtin_clzll(x) : 64;
#else
  // portable fallback
  if (!x)
    return 64;
  int n = 0;
  while ((x >> (63 - n)) == 0)
    ++n;
  return n;
#endif
}

struct GF_2;
/**
 * FieldElement: Wrapper for elements in GF(2^n)
 */
struct FieldElement {
  uint64_t value;
  const GF_2 *field;

  FieldElement(uint64_t v, const GF_2 *f);
  FieldElement() : value(0), field(nullptr) {}

  // Field operations - all operations maintain the field structure
  FieldElement operator+(const FieldElement &other) const;
  FieldElement operator-(const FieldElement &other) const;
  FieldElement operator*(const FieldElement &other) const;
  FieldElement operator/(const FieldElement &other) const;
  bool operator==(const FieldElement &other) const;
};

/**
 * Irreducible Polynomials for GF(2^n) Construction are precomputed and
 * represented in octal notation.
 * source :
 * https://web.eecs.utk.edu/~jplank/plank/papers/CS-07-593/primitive-polynomial-table.txt
 */
static constexpr std::array<uint64_t, 33> IRREDUCIBLE_POLYS = {
    0,            // n=0 (unused)
    0,            // n=1 (unused)
    07,           // n=2
    013,          // n=3
    023,          // n=4
    045,          // n=5
    0103,         // n=6
    0203,         // n=7
    0435,         // n=8
    01021,        // n=9
    02011,        // n=10
    04005,        // n=11
    010123,       // n=12
    020033,       // n=13
    042103,       // n=14
    0100003,      // n=15
    0210013,      // n=16
    0400011,      // n=17
    01000201,     // n=18
    02000047,     // n=19
    04000011,     // n=20
    010000005,    // n=21
    020000003,    // n=22
    040000041,    // n=23
    0100000207,   // n=24
    0200000011,   // n=25
    0400000107,   // n=26
    01000000047,  // n=27
    02000000011,  // n=28
    04000000005,  // n=29
    010040000007, // n=30
    020000000011, // n=31
    040020000007  // n=32
};

struct GF_2 {
  uint64_t n;         // Field degree.
  uint64_t poly;      // The irreducible polynomial defining the field.
  uint64_t omega = 2; // Generator element (primitive root)
  uint64_t fieldSize; // 2^n (total number of field elements)
  uint64_t m;         // fieldSize - 1 (order of multiplicative group)

  explicit GF_2(const uint64_t n_)
      : n(n_), fieldSize(1 << n), m(fieldSize - 1), poly(IRREDUCIBLE_POLYS[n]) {
  }
  [[nodiscard]] uint64_t reduce(uint64_t x) const {
    if (x < fieldSize)
      return x;

    while (x >= fieldSize) {
      auto const shift = 64 - leading_zeros(x) - n;
      x ^= (poly << shift);
    }
    return x;
  }

  [[nodiscard]] static uint64_t add(const uint64_t x, const uint64_t y) {
    return x ^ y; // xor
  }

  [[nodiscard]] uint64_t mul(uint64_t x, uint64_t y) const { // à la russe
    if (x == 0 || y == 0)
      return 0;
    if (x == 1)
      return y;
    if (y == 1)
      return x;

    uint64_t res = 0;
    x = reduce(x);
    y = reduce(y);

    while (x != 0) {
      if (x & 1) {
        res ^= y;
      }
      x >>= 1;
      y <<= 1;
      if (y >= fieldSize) {
        y ^= poly;
      }
    }
    return reduce(res);
  }

  [[nodiscard]] uint64_t powGF(uint64_t base, uint64_t exp) const {
    if (base == 0)
      return (exp == 0) ? 1 : 0;
    if (exp == 0)
      return 1;
    if (exp == 1)
      return reduce(base);

    base = reduce(base);
    exp %= m;

    uint64_t res = 1;
    while (exp > 0) {
      if (exp & 1) {
        res = mul(res, base);
      }
      base = mul(base, base);
      exp >>= 1;
    }
    return reduce(res);
  }

  [[nodiscard]] uint64_t inv(const uint64_t x) const {
    if (x == 0)
      throw std::invalid_argument("zero is not invertible");
    return powGF(x, m - 1); // th de Fermat
  }
};
