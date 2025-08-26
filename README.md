# tree-evaluation

for the python version:
pip install requirements.txt

python proto/tree.py

on macos, linux:
cmake -DCMAKE_BUILD_TYPE=Release -B build/ && cmake --build build/


# Use Ninja (if installed)
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -B ./build && cmake --build ./build

# Or use MinGW Makefiles
cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -B ./build && cmake --build ./build

# Or use MSYS Makefiles (in MSYS2 environment)
cmake -G "MSYS Makefiles" -DCMAKE_BUILD_TYPE=Release -B ./build && cmake --build ./build
