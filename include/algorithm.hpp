/**
 * Tree Evaluation Algorithm Implementation (warm-up version)
 *
 * This file implements the main algorithm from "Tree Evaluation is in Space
 * O(log n \dot log log n)".
 *
 * The algorithm evaluates a complete binary tree where:
 *  - Each internal node computes a function f_u: [k] * [k] -> [k]
 *  - Each leaf contains a value from alphabet [k] = {0, 1, ..., k-1} (and keep
 *  in mind with one node the height is 0 not 1)
 *  - The goal is to compute the value at the root using O(log n * log log n)
 *  space
 *
 * Please refer to the README for more details of the choice
 */

#include "galois_field.hpp"

#include <vector>

/**
 * Algorithm Parameters Structure
 *
 * Contains all the parameters and data structures needed for the tree
 * evaluation. This corresponds to the formal description in the paper.
 */

struct algorithm_parameters {
  uint64_t h; // tree height
  uint64_t k; // alphabet size
  std::vector<std::vector<FieldElement>>
      f_u; // table encoding of internal node functions
  std::vector<std::vector<std::vector<FieldElement>>>
      tree_functions;                    // all function tables per node
  std::vector<FieldElement> tree_leaves; // leaf values (over [k])
  GF_2 *field;                           // pointer to GF(2^n)
  std::vector<std::vector<FieldElement>> registers;
  uint64_t u;     // the u of Goldreich, current node index
  uint64_t log_k; // ⌈log_2(k+1)⌉ (the +1 is explained in the README)
};

extern algorithm_parameters algo;

/**
 * Field Initialization
 *
 * Initializes the finite field GF(2^n) where n is chosen to be large enough
 * to represent all intermediate values during computation.
 *
 * From the paper: n = ⌈log_2(2⌈(log_2(k+1)⌉ + 2) (the + 1 is
 * explained in the README) This ensures the field can represent all possible
 * polynomial coefficients.
 *
 * @param k Alphabet size
 * @return Pointer to initialized finite field
 */
GF_2 *initialize_field(const uint64_t k);

/**
 * Algorithm Parameters Initialization
 *
 * Sets up all data structures needed for the tree evaluation algorithm.
 * Prompts user for tree height, alphabet size, function definitions, and leaf
 * values.
 *
 * @return populated algorithm_parameters
 */
algorithm_parameters algorithm_parameters_init();

/**
 * Bit Vector Conversion
 *
 * Converts an integer x \in [k] to its binary representation as a vector of
 * field elements.
 *
 * The bit vector (b_0, b1, ..., b_{⌈log_2(k+1)⌉ - 1}) (for example 11 is
 * represented as (1, 1, 0, 1) and not (1, 0, 1, 1))
 *
 * @param x Integer to convert (must be < k)
 * @return Vector of field elements representing bits of x
 */
std::vector<FieldElement> get_bit_vector(uint64_t x);

/**
 * Polynomial Evaluation Function e_poly
 *
 * Computes the polynomial e (as in the paper)
 *
 * @param r Register index (0 or 2)
 * @param beta Bit vector to evaluate polynomial at f_u
 * @return Polynomial evaluation result
 */
FieldElement e_poly(int r, const std::vector<FieldElement> &beta);

/**
 * q_u_i computes the polynomial q_{u,i} as in the paper
 *
 *
 * From the paper: q_u_i represents the contribution of bit i (more compicated
 * than that) when evaluating function f_u at the current node u.
 *
 * @param i Bit position to compute
 * @return Field element representing q_{u,i}
 */
FieldElement q_u_i(uint64_t i);

/**
 * Recursive Clean Computation
 *
 * Implements the main recursive algorithm for tree evaluation (Lemma 4.5,
 * warm-up version), register programs which cleanly compute the value at the
 * root. This is the core of the space-efficient algorithm from the paper.
 *
 * Note: registers are rotated as in the Godreich paper
 */
void recursive_clean();

/**
 * Inverse Clean Computation
 *
 * -q_{u,i} is replaced by +q_{u,i}
 *
 */
void recursive_clean_inverse();

/**
 * Runs the recursive clean computation and extracts the final result
 * from the register differences.
 *
 * @return The value computed at the root of the tree
 */
uint64_t clean_computation();

/**
 * Naive Tree Evaluation (for Verification)
 *
 * Implements the straightforward recursive tree evaluation for correctness
 * checking.
 *
 * @param node Current node index
 * @param level Current level in the tree (0 = root)
 * @return Value computed at the given node
 */
FieldElement naive_algorithm(uint64_t node, uint64_t level);
