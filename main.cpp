#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <ranges>
#include <random>

using namespace std;


struct GF_2;


struct FieldElement {
    uint64_t value;
    const GF_2* field;

    FieldElement(uint64_t v, const GF_2* f);
    FieldElement() : value(0), field(nullptr) {}

    FieldElement operator+(const FieldElement& other) const;
    FieldElement operator-(const FieldElement& other) const;
    FieldElement operator*(const FieldElement& other) const;
    FieldElement operator/(const FieldElement& other) const;
    bool operator==(const FieldElement& other) const;
};


struct GF_2 {
    uint64_t n;                 // Field degree.
    uint64_t poly;              // The irreducible polynomial for GF(2^m).
    uint64_t omega;             // The generator element.
    uint64_t fieldSize;         // 2^n
    uint64_t m;                 // fieldSize - 1

    [[nodiscard]] explicit GF_2(const uint64_t n) : n(n){
        fieldSize = 1ULL << n;
        m = fieldSize - 1;

        omega = 2; // if the irreducible polynomial is primitive
        switch (n) { // trinomial polynomials (3 elements) to accelerate the code not all of them though
            case 2: poly = 07; break;
            case 3: poly = 013; break;
            case 4: poly = 023; break;
            case 5: poly = 045; break;
            case 6: poly = 0103; break;
            case 7: poly = 0203; break;
            case 8: poly = 0435; break;
            case 9: poly = 01021; break;
            case 10: poly = 02011; break;
            case 11: poly = 04005; break;
            case 12: poly = 010123; break;
            case 13: poly = 020033; break;
            case 14: poly = 042103; break;
            case 15: poly = 0100003; break;
            case 16: poly = 0210013; break;
            case 17: poly = 0400011; break;
            case 18: poly = 01000201; break;
            case 19: poly = 02000047; break;
            case 20: poly = 04000011; break;
            case 21: poly = 010000005; break;
            case 22: poly = 020000003; break;
            case 23: poly = 040000041; break;
            case 24: poly = 0100000207; break;
            case 25: poly = 0200000011; break;
            case 26: poly = 0400000107; break;
            case 27: poly = 01000000047; break;
            case 28: poly = 02000000011; break;
            case 29: poly = 04000000005; break;
            case 30: poly = 010040000007; break;
            case 31: poly = 020000000011; break;
            case 32: poly = 040020000007; break;
            default: throw invalid_argument("Unsupported for now");
        }
    }

    [[nodiscard]] uint64_t reduce(uint64_t x) const {
        if (x == 0) return 0;
        while(x >= fieldSize) {
            // __builtin_clzll returns the number of leading zeros in a 64-bit integer.
            auto const shift = 64 - __builtin_clzll(x) - n;
            x ^= (poly << shift);
        }
        return x;
    }

    [[nodiscard]] static uint64_t add(const uint64_t x, const uint64_t y) {
        return x ^ y; // xor
    }

    [[nodiscard]] uint64_t mul(uint64_t x, uint64_t y) const { // à la russe
        uint64_t z = 0;
        while (x > 0) {
            if (x & 1) {
                z ^= y;
            }
            x >>= 1;
            y <<= 1;
            if (y & fieldSize) {
                y ^= poly;
            }
        }
        return z;
    }

    [[nodiscard]] uint64_t powGF(uint64_t x, uint64_t exp) const {
        uint64_t result = 1;
        exp %= m;

        while (exp > 0) {
            if (exp & 1)
                result = mul(result, x);
            x = mul(x, x);
            exp >>= 1;
        }
        return result;
    }

     [[nodiscard]] uint64_t inv(const uint64_t x) const {
        if (x == 0)
            throw invalid_argument("zero is not invertible");
        return powGF(x, fieldSize - 2); // th de Fermat
    }
};


FieldElement::FieldElement(uint64_t v, const GF_2* f) : field(f) {
    if (!f) throw invalid_argument("Field pointer cannot be null");
    if (v >= f->fieldSize)
        throw invalid_argument("Element exceeds field size");
    value = v;
}


FieldElement FieldElement::operator+(const FieldElement& other) const {
    if (field != other.field)
        throw invalid_argument("Field mismatch in addition");
    return { GF_2::add(value, other.value), field };
}


FieldElement FieldElement::operator-(const FieldElement& other) const {
    return *this + other;
}


FieldElement FieldElement::operator*(const FieldElement& other) const {
    if (field != other.field)
        throw invalid_argument("Field mismatch in multiplication");
    return {FieldElement(field->mul(value, other.value), field)};
}

FieldElement FieldElement::operator/(const FieldElement& other) const {
    if (field != other.field)
        throw invalid_argument("Field mismatch in division");
    const uint64_t inv_val = field->inv(other.value);
    return {FieldElement(field->mul(value, inv_val), field)};
}

bool FieldElement::operator==(const FieldElement& other) const {
    return field == other.field && value == other.value;
}


struct algorithm_parameters {
    uint64_t h;
    uint64_t k;
    vector<vector<FieldElement>> f_u;
    vector<vector<vector<FieldElement>>> tree_functions;
    vector<FieldElement> tree_leaves;
    GF_2* field;
    vector<vector<FieldElement>> registers;
    uint64_t u;
    uint64_t log_k;
};


GF_2* initialize_field(const uint64_t k) {
    uint64_t field_degree = ceil(log2(2 * ceil(log2(k)) + 2));
    auto const f = new GF_2(field_degree);
    return f;
}


algorithm_parameters algorithm_parameters_init() {
    algorithm_parameters params;
    cout << "Enter the height of the tree: ";
    cin >> params.h;
    cout << "Enter the size of the alphabet (k): ";
    cin >> params.k;

    uint64_t h = params.h;
    uint64_t k = params.k;

    params.field = initialize_field(k);
    vector f_u(k, vector<FieldElement>(k));
    cout << "Enter the values for the function matrix:" << endl;
    for (uint64_t i = 0; i < k; i++) {
        for (uint64_t j = 0; j < k; j++) {
            uint64_t val;
            cout << "Row " << i << ", Column " << j << ": ";
            cin >> val;
            f_u[i][j] = FieldElement(val, params.field);
        }
    }

    params.f_u = f_u;

    params.tree_functions.resize((1ULL << h) - 1, f_u);
    params.tree_leaves.resize(1ULL << h, FieldElement(0, params.field));
    cout << "Enter the value for each leaf node:" << endl;
    for (uint64_t i = 0; i < (1ULL << h); i++) {
        uint64_t leaf_val;
        cout << "Leaf node " << i << ": ";
        cin >> leaf_val;
        params.tree_leaves[i] = FieldElement(leaf_val, params.field);
    }

    uint64_t log_k = ceil(log2(k + 1));
    params.log_k = log_k;

    params.registers.resize(3);
    for (int r = 0; r < 3; r++) {
        params.registers[r].resize(log_k);
        for (uint64_t i = 0; i < log_k; i++) {
            params.registers[r][i] = FieldElement(rand() % params.field->fieldSize, params.field);
        }
    }

    params.u = 0;
    return params;
}


algorithm_parameters algo = algorithm_parameters_init();

vector<FieldElement> get_bit_vector(uint64_t x) {
    vector<FieldElement> bits;
    bits.reserve(algo.log_k);
    for (uint64_t i = 0; i < algo.log_k; ++i) {
        bits.emplace_back((x >> i) & 1ULL, algo.field);
    }
    return bits;
}

FieldElement e_poly(int r, const vector<FieldElement>& beta) {
    FieldElement res(1, algo.field);
    FieldElement one(1, algo.field);
    for (uint64_t i = 0; i < algo.log_k; ++i) {
        res = res * (one - algo.registers[r][i] + (algo.registers[r][i] + algo.registers[r][i] - one) * beta[i]);
    }
    return res;
}


FieldElement q_u_i(const uint64_t node, const int i) {
    FieldElement res(0, algo.field);
    for (uint64_t alpha = 0; alpha < algo.k; ++alpha) {
        if (((alpha >> i) & 1ULL) != 1)
            continue;
        for (uint64_t beta = 0; beta < algo.k; ++beta) {
            for (uint64_t gamma = 0; gamma < algo.k; ++gamma) {
                if (algo.tree_functions[node][beta][gamma].value == alpha) {
                    vector<FieldElement> beta_bits = get_bit_vector(beta);
                    vector<FieldElement> gamma_bits = get_bit_vector(gamma);
                    res = res + (e_poly(0, beta_bits) * e_poly(2, gamma_bits));
                }
            }
        }
    }
    return res;
}


void recursive_clean_inverse(uint64_t h);


void recursive_clean(uint64_t h) {
    if (h == 0) {
        uint64_t leaf_index = algo.u - ((1ULL << algo.h) - 1);
        for (uint64_t i = 0; i < algo.log_k; ++i) {
            const uint64_t bit = (algo.tree_leaves[leaf_index].value >> i) & 1ULL;
            algo.registers[2][i] = algo.registers[2][i] + FieldElement(bit, algo.field);
        }
        return;
    }

    for (int j = 1; j <= algo.field->m; ++j) {
        auto w_j = FieldElement(algo.field->powGF(algo.field->omega, j), algo.field);

        for (int i = 0; i < algo.log_k; ++i) {
            algo.registers[0][i] = algo.registers[0][i] * w_j;
            algo.registers[1][i] = algo.registers[1][i] * w_j;
        }

        ranges::rotate(algo.registers, algo.registers.begin() + 1);

        algo.u = 2 * algo.u + 1;

        recursive_clean(h - 1);

        swap(algo.registers[0], algo.registers[2]);

        algo.u += 1;

        recursive_clean(h - 1);

        for (int i = 0; i < algo.log_k; ++i) {
            algo.registers[1][i] = algo.registers[1][i] - q_u_i((algo.u - 1) / 2, i);
        }

        swap(algo.registers[0], algo.registers[2]);

        algo.u -= 1;
        recursive_clean_inverse(h - 1);

        swap(algo.registers[0], algo.registers[2]);

        algo.u += 1;

        recursive_clean_inverse(h - 1);

        algo.u = (algo.u - 1) / 2;
        ranges::rotate(algo.registers, algo.registers.begin() + 2);
        swap(algo.registers[0], algo.registers[1]);

        for (int i = 0; i < algo.log_k; ++i) {
            algo.registers[0][i] = algo.registers[0][i] / w_j;
            algo.registers[1][i] = algo.registers[1][i] / w_j;
        }
    }
}

void recursive_clean_inverse(uint64_t h) {
    if (h == 0) {
        uint64_t leaf_index = algo.u - ((1ULL << algo.h) - 1);
        for (uint64_t i = 0; i < algo.log_k; ++i) {
            const uint64_t bit = (algo.tree_leaves[leaf_index].value >> i) & 1ULL;
            algo.registers[2][i] = algo.registers[2][i] - FieldElement(bit, algo.field);
        }
        return;
    }

    for (int j = 1; j <= algo.field->m; ++j) {
        auto w_j = FieldElement(algo.field->powGF(algo.field->omega, j), algo.field);

        for (int i = 0; i < algo.log_k; ++i) {
            algo.registers[0][i] = algo.registers[0][i] * w_j;
            algo.registers[1][i] = algo.registers[1][i] * w_j;
        }

        ranges::rotate(algo.registers, algo.registers.begin() + 1);

        algo.u = 2 * algo.u + 1;

        recursive_clean_inverse(h - 1);

        swap(algo.registers[0], algo.registers[2]);

        algo.u += 1;

        recursive_clean_inverse(h - 1);

        for (int i = 0; i < algo.log_k; ++i) {
            algo.registers[1][i] = algo.registers[1][i] + q_u_i((algo.u - 1) / 2, i);
        }

        swap(algo.registers[0], algo.registers[2]);

        algo.u -= 1;
        recursive_clean(h - 1);

        swap(algo.registers[0], algo.registers[2]);

        algo.u += 1;

        recursive_clean(h - 1);

        algo.u = (algo.u - 1) / 2;
        ranges::rotate(algo.registers, algo.registers.begin() + 2);
        swap(algo.registers[0], algo.registers[1]);

        for (int i = 0; i < algo.log_k; ++i) {
            algo.registers[0][i] = algo.registers[0][i] / w_j;
            algo.registers[1][i] = algo.registers[1][i] / w_j;
        }
    }
}

uint64_t clean_computation(uint64_t h) {
    auto registers_copy = algo.registers;

    recursive_clean(h);

    vector<FieldElement> diff;
    uint64_t res = 0;

    diff.resize(algo.log_k);

    for (uint64_t i = 0; i < algo.log_k; i++) {
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


int main() {
    const uint64_t h = algo.h;
    const uint64_t result = clean_computation(h);

    const auto naive = naive_algorithm(0, 0);

    assert(result == naive.value && "Result does not match the naive evaluation!!!!!");
    cout << "Tree evaluation result: " << result << endl;
    return 0;
}
