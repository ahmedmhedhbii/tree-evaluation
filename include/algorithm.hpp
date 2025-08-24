#include "galois_field.hpp"

#include <vector>

struct algorithm_parameters {
  uint64_t h;
  uint64_t k;
  std::vector<std::vector<FieldElement>> f_u;
  std::vector<std::vector<std::vector<FieldElement>>> tree_functions;
  std::vector<FieldElement> tree_leaves;
  GF_2 *field;
  std::vector<std::vector<FieldElement>> registers;
  uint64_t u;
  uint64_t log_k;
};

extern algorithm_parameters algo;

GF_2 *initialize_field(const uint64_t k);
algorithm_parameters algorithm_parameters_init();
std::vector<FieldElement> get_bit_vector(uint64_t x);
FieldElement e_poly(int r, const std::vector<FieldElement> &beta);
FieldElement q_u_i(uint64_t i);
void recursive_clean();
void recursive_clean_inverse();
uint64_t clean_computation();
FieldElement naive_algorithm(uint64_t node, uint64_t level);
