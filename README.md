# Tree Evaluation in Space $\mathcal{O}(\log n \cdot \log \ \log n)$ 

This repository contains implementations of the tree evaluation algorithm from the paper **"Tree Evaluation is in Space $\mathcal{O}(\log n \cdot \log \ \log n)$"** (warm-up version) by James Cook & Ian Mertz. 

## Overview 
The algorithm evaluates a complete binary tree where: 

- Each internal node computes a function $f_u: [k] × [k] \rightarrow [k]$. 
- Each leaf contains a value from alphabet $[k] = \{0, 1, \cdots, k-1 \}$. 
- The goal is to compute the value at the root using $\mathcal{O}(\log n \cdot \log \ \log n)$ space. 
- Tree height is measured such that a single node has height 0. 

**Note**: 
- the variable log_k is $\lceil \log (k+1) \rceil$ and not $\lceil \log k \rceil$.
- **Registers**: Rotated as in Goldreich paper, three-register system (x,y,z).
- **Bit vectors**: Little-endian format (11 $\rightarrow$ (1,1,0,1)).

## Quick Start 

### Python

```bash
pip install galois==0.4.3
python proto/tree.py
```

### C++ (C++20 or newer)
#### macOS / Linux:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/tree
```


### Windows (macOS / Linux too):
```bash
g++ -std=c++20 -O3 -Iinclude \
    src/main.cpp src/galois_field.cpp src/algorithm.cpp \ # no need for tests/test.cpp
    -o tree.exe
```
## Example of usage:
```txt
Enter the height of the tree: 2
Enter the size of the alphabet (k): 3
Enter the values for the function matrix:
Row 0, Column 0: 1
Row 0, Column 1: 2
Row 0, Column 2: 0
Row 1, Column 0: 0
Row 1, Column 1: 1
Row 1, Column 2: 2
Row 2, Column 0: 2
Row 2, Column 1: 0
Row 2, Column 2: 1
Enter the value for each leaf node:
Leaf node 0: 1
Leaf node 1: 0
Leaf node 2: 2
Leaf node 3: 1
Tree evaluation result: 1
```
Tree representation:
```txt
                                    +--------------+
                                    |              |
                     +--------------+      f_u     +------------------+
                     |              |              |                  |
                     |              +--------------+                  |
                     |                                                |
              +------+-----+                                  +-------+-------+
              |      f_u   |                                  |      f_u      |
      +-------+            +---------+                +-------+               +----------+
      |       +------------+         |                |       +---------------+          |
+-----+-----+                 +------+------+  +------+------+                    +------+------+
|     1     |                 |      0      |  |      2      |                    |      1      |
+-----------+                 +-------------+  +-------------+                    +-------------+

```
