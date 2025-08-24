#include "../include/galois_field.hpp"

inline FieldElement::FieldElement(uint64_t v, const GF_2 *f) : field(f) {
  if (!f)
    throw std::invalid_argument("Field pointer cannot be null");
  value = f->reduce(v);
}

inline FieldElement FieldElement::operator+(const FieldElement &other) const {
  if (field != other.field)
    throw std::invalid_argument("Field mismatch in addition");
  return {GF_2::add(value, other.value), field};
}

inline FieldElement FieldElement::operator-(const FieldElement &other) const {
  // In GF(2^n), subtraction is the same as addition
  return *this + other;
}

inline FieldElement FieldElement::operator*(const FieldElement &other) const {
  if (field != other.field)
    throw std::invalid_argument("Field mismatch in multiplication");
  return {field->mul(value, other.value), field};
}

inline FieldElement FieldElement::operator/(const FieldElement &other) const {
  if (field != other.field)
    throw std::invalid_argument("Field mismatch in division");

  if (value == 0 || other.value == 0) {
    return {0, field};
  }

  const uint64_t inv_val = field->inv(other.value);
  return {field->mul(value, inv_val), field};
}

inline bool FieldElement::operator==(const FieldElement &other) const {
  return field == other.field && value == other.value;
}
