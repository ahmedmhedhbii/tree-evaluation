#include "../include/algorithm.hpp"
#include <iostream>

extern algorithm_parameters algo;
int main() {

  algo = algorithm_parameters_init();

  const uint64_t result = clean_computation();

  const auto naive = naive_algorithm(0, 0);

  assert(result == naive.value &&
         "Result does not match the naive evaluation!!!!!");

  std::cout << "Tree evaluation result: " << result << std::endl;
  return 0;
}
