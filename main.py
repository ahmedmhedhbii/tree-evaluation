import galois
import math
import random
import copy

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

def e_poly(y_vals, beta):
    res = field(1)
    for yi, betai in zip(y_vals, beta):
        res *= (field(1) - yi + (2 * yi - field(1)) * betai)
    return res

def q_u_i(y, z, f_u, i):
    res = field(0)
    for alpha in range(k):
        for beta in range(k):
            for gamma in range(k):
                alpha_bits = [field(int(bit)) for bit in f'{alpha:0{log_k}b}'[::-1]]
                beta_bits  = [field(int(bit)) for bit in f'{beta:0{log_k}b}'[::-1]]
                gamma_bits = [field(int(bit)) for bit in f'{gamma:0{log_k}b}'[::-1]]
                if alpha_bits[i] == field(1) and f_u[beta][gamma] == alpha:
                    res += e_poly(y, beta_bits) * e_poly(z, gamma_bits)
    return res

def recursive_clean(node, h):
    if h == 0:
        bits = [field(int(bit)) for bit in f'{tree[node]:0{log_k}b}'[::-1]]
        for i, bit in enumerate(bits):
            registers[2][i] += bit #
        return

    left_child = 2 * node + 1
    right_child = 2 * node + 2

    for j in range(1, m + 1):
        omega_j = omega ** j
        for i in range(log_k):
            registers[0][i] *= omega_j
            registers[1][i] *= omega_j

        rotate_left_inplace(registers)
        recursive_clean(left_child, h - 1)

        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean(right_child, h - 1)
        for i in range(log_k):
            registers[1][i] -= q_u_i(registers[0], registers[2], tree[node], i)

        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean_inverse(left_child, h - 1)

        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean_inverse(right_child, h - 1)

        rotate_left_inplace(registers, 2)
        registers[0], registers[1] = registers[1], registers[0]

        for i in range(log_k):
            registers[0][i] /= omega_j
            registers[1][i] /= omega_j

def recursive_clean_inverse(node, h):
    if h == 0:
        bits = [field(int(bit)) for bit in f'{tree[node]:0{log_k}b}'[::-1]]
        for i, bit in enumerate(bits):
            registers[2][i] -= bit
        return

    left_child = 2 * node + 1
    right_child = 2 * node + 2

    for j in range(1, m + 1):
        omega_j = omega ** j
        for i in range(log_k):
            registers[0][i] *= omega_j
            registers[1][i] *= omega_j

        rotate_left_inplace(registers)
        recursive_clean_inverse(left_child, h - 1)
        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean_inverse(right_child, h - 1)
        for i in range(log_k):
            registers[1][i] += q_u_i(registers[0], registers[2], tree[node], i)
        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean(left_child, h - 1)
        registers[0], registers[2] = registers[2], registers[0]
        recursive_clean(right_child, h - 1)
        rotate_left_inplace(registers, 2)
        registers[0], registers[1] = registers[1], registers[0]
        for i in range(log_k):
            registers[0][i] /= omega_j
            registers[1][i] /= omega_j


def clean_computation(h):
    initial_registers = copy.deepcopy(registers)
    recursive_clean(0, h)
    res = 0
    for i in range(log_k):
        res += int(registers[2][i] - initial_registers[2][i]) * (2 ** i)
    return res, registers


def initialize_field(k):
    field_order = 1 << math.ceil(math.log2(2 * math.ceil(math.log2(k)) + 2))
    field = galois.GF(field_order)
    return field, field_order


def initialize_tree_and_catalyst():
    h = int(input("Enter the height of the tree: "))
    k = int(input("Enter the size of the alphabet: "))
    num_nodes = (1 << (h + 1)) - 1
    function = [[int(input(f"Enter the value for function at row {r}, column {c}: ")) for c in range(k)]
                for r in range(k)]
    tree_functions = [function for _ in range((1 << h) - 1)]

    tree_values = [int(input(f"Enter the value of leaf node {i}: ")) for i in range(1 << h)]
    tree = tree_functions + tree_values

    log_k = math.ceil(math.log2(k + 1))

    field, field_order = initialize_field(k)
    registers = [[field(random.randint(0, field.order - 1)) for _ in range(log_k)] for _ in
                 range(num_nodes)]

    return h, k, tree, registers, field, field_order - 1


def the_not_so_cool_algorithm(node, level):
    if level == height:
        return tree[node]

    left = the_not_so_cool_algorithm(2 * node + 1, level + 1)
    right = the_not_so_cool_algorithm(2 * node + 2, level + 1)

    return tree[node][left][right]

if __name__ == "__main__":
    height, k, tree, registers, field, m = initialize_tree_and_catalyst()
    log_k = math.ceil(math.log2(k + 1))
    omega = field.primitive_element
    registers_copy = copy.deepcopy(registers)
    result, final_registers = clean_computation(height)

    print(result)
    naive = the_not_so_cool_algorithm(0, 0)
    print(naive)

    print("*" * 80)
    assert result == naive
    print(f"Tree evaluation result: {result}")
    print(f"Initial registers:\n{registers_copy}")
    print(f"Final registers:\n{final_registers}")


