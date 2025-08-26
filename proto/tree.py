"""
Tree Evaluation Algorithm Implementation (warm-up version)

This file implements the main algorithm from "Tree Evaluation is in Space
O(log n \dot log log n)".

The algorithm evaluates a complete binary tree where:
 - Each internal node computes a function f_u: [k] * [k] -> [k]
 - Each leaf contains a value from alphabet [k] = {0, 1, ..., k-1} (and keep
 in mind with one node the height is 0 not 1)
 - The goal is to compute the value at the root using O(log n * log log n)
 space

Please refer to the README for more details of the choice
"""

import copy
import math
import random

import galois

# Global algorithm parameters (corresponds to algorithm_parameters struct in C++)
height: int = 0  # tree height (h)
k = 0  # alphabet size
tree = []  # tree structure: functions + leaf values
registers = []  # the three registers used in the algorithm
field: galois.GF = None  # pointer to GF(2^n) field
log_k = 0  # ⌈log_2(k+1)⌉ (the +1 is explained in the README)
u = 0  # the u of Goldreich, current node index
omega: galois.GF = None  # generator element (primitive root)
m = 0  # field_order - 1 (order of multiplicative group)


def reverse(lst, start, end):
    while start < end:
        lst[start], lst[end] = lst[end], lst[start]
        start += 1
        end -= 1


def rotate_left_inplace(lst, shift=1):
    n = len(lst)
    shift %= n
    if shift == 0:
        return
    reverse(lst, 0, shift - 1)
    reverse(lst, shift, n - 1)
    reverse(lst, 0, n - 1)


def e_poly(r, beta):
    """
    Polynomial Evaluation Function e_poly

    Computes the polynomial e (as in the paper)

    Args:
        r: Register index (0 or 2)
        beta: Bit vector to evaluate polynomial at f_u

    Returns:
        Polynomial evaluation result
    """
    res = field(1)
    one = field(1)
    for i in range(log_k):
        res = res * (
            one - registers[r][i] + (2 * registers[r][i] - one) * beta[i]
        )
    return res


def q_u_i(i):
    """
    q_u_i computes the polynomial q_{u,i} as in the paper

    From the paper: q_u_i represents the contribution of bit i (more complicated
    than that) when evaluating function f_u at the current node u.

    Args:
        i: Bit position to compute

    Returns:
        Field element representing q_{u,i}
    """
    global u
    res = field(0)
    for alpha in range(k):
        if ((alpha >> i) & 1) != 1:
            continue
        for beta in range(k):
            for gamma in range(k):
                beta_bits = [
                    field(int(bit)) for bit in f"{beta:0{log_k}b}"[::-1]
                ]
                gamma_bits = [
                    field(int(bit)) for bit in f"{gamma:0{log_k}b}"[::-1]
                ]
                if tree[u][beta][gamma] == alpha:
                    res = res + (e_poly(0, beta_bits) * e_poly(2, gamma_bits))
    return res


def recursive_clean():
    """
    Recursive Clean Computation

    Implements the main recursive algorithm for tree evaluation (Lemma 4.5,
    warm-up version), register programs which cleanly compute the value at the
    root. This is the core of the space-efficient algorithm from the paper.

    Note: registers are rotated as in the Godreich paper
    """
    global u
    if u >= (1 << height) - 1:
        bits = [field(int(bit)) for bit in f"{tree[u]:0{log_k}b}"[::-1]]
        for i, bit in enumerate(bits):
            registers[2][i] += bit
        return

    for j in range(1, m + 1):
        omega_j = omega**j

        for i in range(log_k):
            registers[0][i] = registers[0][i] * omega_j
            registers[1][i] = registers[1][i] * omega_j

        rotate_left_inplace(registers)

        u = u * 2 + 1

        recursive_clean()
        u = u + 1
        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean()
        u = (u - 1) // 2
        for i in range(log_k):
            registers[1][i] = registers[1][i] - q_u_i(i)

        registers[0], registers[2] = registers[2], registers[0]
        u = 2 * u + 1
        recursive_clean_inverse()

        registers[0], registers[2] = registers[2], registers[0]
        u = u + 1
        recursive_clean_inverse()

        u = (u - 1) // 2
        rotate_left_inplace(registers, 2)
        registers[0], registers[1] = registers[1], registers[0]

        for i in range(log_k):
            registers[0][i] = registers[0][i] / omega_j
            registers[1][i] = registers[1][i] / omega_j


def recursive_clean_inverse():
    """
    Inverse Clean Computation

    -q_{u,i} is replaced by +q_{u,i}
    This undoes the operations performed by recursive_clean.
    """
    global u
    if u >= (1 << height) - 1:
        bits = [field(int(bit)) for bit in f"{tree[u]:0{log_k}b}"[::-1]]
        for i, bit in enumerate(bits):
            registers[2][i] = registers[2][i] - bit
        return

    for j in range(1, m + 1):
        omega_j = omega**j

        for i in range(log_k):
            registers[0][i] = registers[0][i] * omega_j
            registers[1][i] = registers[1][i] * omega_j

        rotate_left_inplace(registers)

        u = u * 2 + 1

        recursive_clean_inverse()
        u = u + 1
        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean_inverse()
        u = (u - 1) // 2
        for i in range(log_k):
            registers[1][i] = registers[1][i] + q_u_i(i)

        registers[0], registers[2] = registers[2], registers[0]
        u = 2 * u + 1
        recursive_clean()

        registers[0], registers[2] = registers[2], registers[0]
        u = u + 1
        recursive_clean()

        u = (u - 1) // 2
        rotate_left_inplace(registers, 2)
        registers[0], registers[1] = registers[1], registers[0]

        for i in range(log_k):
            registers[0][i] = registers[0][i] / omega_j
            registers[1][i] = registers[1][i] / omega_j


def clean_computation():
    """
    Runs the recursive clean computation and extracts the final result
    from the register differences.

    Returns:
        The value computed at the root of the tree
    """
    initial_registers = copy.deepcopy(registers)
    recursive_clean()
    res = 0
    for i in range(log_k):
        assert int(registers[0][i]) == int(initial_registers[0][i])
        assert int(registers[1][i]) == int(initial_registers[1][i])
        res += int(registers[2][i] - initial_registers[2][i]) * (1 << i)
    return res, registers


def initialize_field(k):
    """
    Field Initialization

    Initializes the finite field GF(2^n) where n is chosen to be large enough
    to represent all intermediate values during computation.

    From the paper: n = ⌈log_2(2⌈(log_2(k+1)⌉ + 2) (the + 1 is
    explained in the README) This ensures the field can represent all possible
    polynomial coefficients.

    Args:
        k: Alphabet size

    Returns:
        Tuple of (field object, field_order - 1)
    """
    field_order = 1 << math.ceil(math.log2(2 * math.ceil(math.log2(k + 1)) + 2))
    field = galois.GF(field_order, repr="poly")
    return field, field_order


def initialize_tree_and_catalyst():
    """
    Algorithm Parameters Initialization

    Sets up all data structures needed for the tree evaluation algorithm.
    Prompts user for tree height, alphabet size, function definitions, and leaf
    values.

    Returns:
        Tuple of (height, k, tree, registers, field, m)
    """
    h = int(input("Enter the height of the tree: "))
    k = int(input("Enter the size of the alphabet: "))
    function = [
        [
            int(input(f"Enter the value for function at row {r}, column {c}: "))
            for c in range(k)
        ]
        for r in range(k)
    ]

    tree_functions = [function for _ in range((1 << h) - 1)]

    tree_values = [
        int(input(f"Enter the value of leaf node {i}: ")) for i in range(1 << h)
    ]
    tree = tree_functions + tree_values
    log_k = math.ceil(math.log2(k + 1))

    field, field_order = initialize_field(k)

    registers = [
        [field(random.randint(0, field.order - 1)) for _ in range(log_k)]
        for _ in range(3)
    ]  # field.Random((3, log_k))
    return h, k, tree, registers, field, field_order - 1


def naive_algorithm(node, h):
    """
    Naive Tree Evaluation (for Verification)

    Implements the straightforward recursive tree evaluation for correctness
    checking.

    Args:
        node: Current node index
        h: Current height (0 = leaf level)
    Returns:
        Value computed at the given node
    """
    if h == 0:
        return tree[node]

    left = naive_algorithm(2 * node + 1, h - 1)
    right = naive_algorithm(2 * node + 2, h - 1)

    return tree[node][left][right]


if __name__ == "__main__":
    # Initialize all global parameters
    height, k, tree, registers, field, m = initialize_tree_and_catalyst()
    log_k = math.ceil(math.log2(k + 1))
    omega = field.primitive_element
    registers_copy = copy.deepcopy(registers)

    result, final_registers = clean_computation()

    naive = naive_algorithm(0, height)

    print("*" * 80)
    assert result == naive
    print(f"Tree evaluation result: {result}")
    print(f"Initial registers:\n{registers_copy}")
    print(f"Final registers:\n{final_registers}")
