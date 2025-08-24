#include "../include/algorithm.hpp"
#include <iostream>
#include <random>

void test() {
  std::random_device rd;
  std::mt19937 gen(rd());

  int total = 0, passed = 0;
  for (uint64_t h = 2; h < 6; ++h)
    for (uint64_t k = 2; k < 100; ++k) {
      for (int i = 0; i < 100; i++) {
        algorithm_parameters params;
        params.h = h;
        params.k = k;
        params.field = initialize_field(k);

        std::uniform_int_distribution<uint64_t> distrib(0, k - 1);
        params.f_u.resize(k, std::vector<FieldElement>(k)); // k * k
        for (uint64_t i = 0; i < k; ++i)
          for (uint64_t j = 0; j < k; ++j)
            params.f_u[i][j] =
                FieldElement(distrib(gen), params.field); // k * k

        params.tree_functions.resize((1ULL << h) - 1, params.f_u);

        params.tree_leaves.resize(1ULL << h);
        for (uint64_t i = 0; i < (1ULL << h); ++i)
          params.tree_leaves[i] =
              FieldElement(distrib(gen), params.field); // k * k

        params.log_k = ceil(log2(k + 1));

        params.registers.resize(3);

        for (int r = 0; r < 3; r++) {
          params.registers[r].resize(params.log_k);
          for (uint64_t i = 0; i < params.log_k; i++) {
            params.registers[r][i] = FieldElement(distrib(gen), params.field);
          }
        }

        params.u = 0;

        extern algorithm_parameters algo;
        algo = params;

        uint64_t result = clean_computation();
        FieldElement naive = naive_algorithm(0, 0);

        ++total;
        if (result == naive.value) {
          std::cout << "[Success] k=" << k << ", h=" << h
                    << " | clean=" << result << ", naive=" << naive.value
                    << std::endl;
          ++passed;
        } else {
          std::cout << "[FAIL] k=" << k << ", h=" << h << " | clean=" << result
                    << ", naive=" << naive.value << std::endl;
        }
      }
    }
  std::cout << "[SUMMARY] Passed " << passed << " / " << total << " tests."
            << std::endl;
}

int main() { test(); }
