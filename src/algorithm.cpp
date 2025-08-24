#include "../include/algorithm.hpp"
#include <algorithm>
#include <iostream>
#include <random>

algorithm_parameters algo;

GF_2 *initialize_field(const uint64_t k) {
  uint64_t field_degree = ceil(log2(2 * ceil(log2(k + 1)) + 2));
  GF_2 *const f = new GF_2(field_degree);
  return f;
}

algorithm_parameters algorithm_parameters_init() {
  algorithm_parameters params;
  std::cout << "Enter the height of the tree: ";
  std::cin >> params.h;
  std::cout << "Enter the size of the alphabet (k): ";
  std::cin >> params.k;

  uint64_t h = params.h;
  uint64_t k = params.k;

  params.field = initialize_field(k);

  std::vector f_u(k, std::vector<FieldElement>(k));
  std::cout << "Enter the values for the function matrix:" << std::endl;
  for (uint64_t i = 0; i < k; i++) {
    for (uint64_t j = 0; j < k; j++) {
      uint64_t val;
      std::cout << "Row " << i << ", Column " << j << ": ";
      std::cin >> val;
      f_u[i][j] = FieldElement(val, params.field);
    }
  }

  params.f_u = f_u;

  params.tree_functions.resize((1ULL << h) - 1, f_u);
  params.tree_leaves.resize(1ULL << h, FieldElement(0, params.field));

  std::cout << "Enter the value for each leaf node:" << std::endl;
  for (uint64_t i = 0; i < (1ULL << h); i++) {
    uint64_t leaf_val;
    std::cout << "Leaf node " << i << ": ";
    std::cin >> leaf_val;
    params.tree_leaves[i] = FieldElement(leaf_val, params.field);
  }

  uint64_t log_k = ceil(log2(k + 1));
  params.log_k = log_k;

  params.registers.resize(3);
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<uint64_t> distrib(0,
                                                  params.field->fieldSize - 1);
  for (int r = 0; r < 3; r++) {
    params.registers[r].resize(log_k);
    for (uint64_t i = 0; i < log_k; i++) {
      params.registers[r][i] = FieldElement(distrib(gen), params.field);
    }
  }

  params.u = 0;
  return params;
}

inline std::vector<FieldElement> get_bit_vector(uint64_t x) {
  std::vector<FieldElement> bits;
  bits.reserve(algo.log_k);
  for (uint64_t i = 0; i < algo.log_k; ++i) {
    bits.emplace_back((x >> i) & 1ULL, algo.field);
  }
  return bits;
}

FieldElement e_poly(int r, const std::vector<FieldElement> &beta) {
  FieldElement res{1, algo.field};
  FieldElement one{1, algo.field};
  for (uint64_t i = 0; i < algo.log_k; ++i) {
    res = res * (one - algo.registers[r][i] +
                 (algo.registers[r][i] + algo.registers[r][i] - one) * beta[i]);
  }
  return res;
}

FieldElement q_u_i(uint64_t i) {
  FieldElement res{0, algo.field};

  for (uint64_t alpha = 0; alpha < algo.k; ++alpha) {
    if (((alpha >> i) & 1ULL) != 1) // small optimization
      continue;
    for (uint64_t beta = 0; beta < algo.k; ++beta) {
      for (uint64_t gamma = 0; gamma < algo.k; ++gamma) {
        if (algo.tree_functions[algo.u][beta][gamma].value == alpha) {
          std::vector<FieldElement> beta_bits = get_bit_vector(beta);
          std::vector<FieldElement> gamma_bits = get_bit_vector(gamma);
          res = res + (e_poly(0, beta_bits) * e_poly(2, gamma_bits));
        }
      }
    }
  }
  return res;
}

void recursive_clean() {
  if (algo.u >= (1ULL << algo.h) - 1) {
    uint64_t leaf_index = algo.u - ((1ULL << algo.h) - 1);
    for (uint64_t i = 0; i < algo.log_k; ++i) {
      const uint64_t bit = (algo.tree_leaves[leaf_index].value >> i) & 1ULL;
      algo.registers[2][i] =
          algo.registers[2][i] + FieldElement(bit, algo.field);
    }
    return;
  }

  for (int j = 1; j <= algo.field->m; ++j) {
    auto w_j =
        FieldElement(algo.field->powGF(algo.field->omega, j), algo.field);

    for (int i = 0; i < algo.log_k; ++i) {
      algo.registers[0][i] = algo.registers[0][i] * w_j;
      algo.registers[1][i] = algo.registers[1][i] * w_j;
    }

    std::ranges::rotate(algo.registers, algo.registers.begin() + 1);

    algo.u = 2 * algo.u + 1;

    recursive_clean();

    swap(algo.registers[0], algo.registers[2]);

    algo.u += 1;

    recursive_clean();

    algo.u = (algo.u - 1) / 2;
    for (int i = 0; i < algo.log_k; ++i) {
      algo.registers[1][i] = algo.registers[1][i] - q_u_i(i);
    }

    swap(algo.registers[0], algo.registers[2]);

    algo.u = algo.u * 2 + 1;
    recursive_clean_inverse();

    swap(algo.registers[0], algo.registers[2]);

    algo.u += 1;

    recursive_clean_inverse();

    algo.u = (algo.u - 1) / 2;
    std::ranges::rotate(algo.registers, algo.registers.begin() + 2);
    swap(algo.registers[0], algo.registers[1]);

    for (int i = 0; i < algo.log_k; ++i) {
      algo.registers[0][i] = algo.registers[0][i] / w_j;
      algo.registers[1][i] = algo.registers[1][i] / w_j;
    }
  }
}

void recursive_clean_inverse() {
  if (algo.u >= (1ULL << algo.h) - 1) {
    uint64_t leaf_index = algo.u - ((1ULL << algo.h) - 1);
    for (uint64_t i = 0; i < algo.log_k; ++i) {
      const uint64_t bit = (algo.tree_leaves[leaf_index].value >> i) & 1ULL;
      algo.registers[2][i] =
          algo.registers[2][i] - FieldElement(bit, algo.field);
    }
    return;
  }

  for (int j = 1; j <= algo.field->m; ++j) {
    auto w_j =
        FieldElement(algo.field->powGF(algo.field->omega, j), algo.field);

    for (int i = 0; i < algo.log_k; ++i) {
      algo.registers[0][i] = algo.registers[0][i] * w_j;
      algo.registers[1][i] = algo.registers[1][i] * w_j;
    }

    std::ranges::rotate(algo.registers, algo.registers.begin() + 1);

    algo.u = 2 * algo.u + 1;

    recursive_clean_inverse();

    swap(algo.registers[0], algo.registers[2]);

    algo.u += 1;

    recursive_clean_inverse();

    algo.u = (algo.u - 1) / 2;
    for (int i = 0; i < algo.log_k; ++i) {
      algo.registers[1][i] = algo.registers[1][i] + q_u_i(i);
    }

    swap(algo.registers[0], algo.registers[2]);

    algo.u = algo.u * 2 + 1;
    recursive_clean();

    swap(algo.registers[0], algo.registers[2]);

    algo.u += 1;

    recursive_clean();

    algo.u = (algo.u - 1) / 2;
    std::ranges::rotate(algo.registers, algo.registers.begin() + 2);
    swap(algo.registers[0], algo.registers[1]);

    for (int i = 0; i < algo.log_k; ++i) {
      algo.registers[0][i] = algo.registers[0][i] / w_j;
      algo.registers[1][i] = algo.registers[1][i] / w_j;
    }
  }
}

uint64_t clean_computation() {
  auto registers_copy = algo.registers;

  recursive_clean();

  std::vector<FieldElement> diff;
  uint64_t res = 0;

  diff.resize(algo.log_k);

  for (uint64_t i = 0; i < algo.log_k; ++i) {
    assert(algo.registers[0][i].value == registers_copy[0][i].value);
    assert(algo.registers[1][i].value == registers_copy[1][i].value);
    res += (algo.registers[2][i] - registers_copy[2][i]).value * (1 << i);
  }

  return res;
}

FieldElement naive_algorithm(uint64_t node, uint64_t level) {
  if (level == algo.h) {
    const uint64_t leaf_index = node - ((1ULL << algo.h) - 1);
    return algo.tree_leaves[leaf_index];
  }

  FieldElement left = naive_algorithm(2 * node + 1, level + 1);
  FieldElement right = naive_algorithm(2 * node + 2, level + 1);
  return algo.tree_functions[node][left.value][right.value];
}
