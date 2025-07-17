# Directed Graph Edge Collapser

This contains 2 header files:
* `DIR_EDGE_COLLAPSER.h` - Main version of the algorithm
* `DIR_EDGE_COLLAPSER_MEM.h` - Version 2 which uses bitset, need to specify the maximum number of vertices in the line `constexpr size_t MAX_VERTICES = 10000;`. Highly memory efficient, mostly used for high density graphs. Takes more time than v1.

The main files:
* `test_with_mem.cpp` - This runs the algo once, displaying the memory usage statistics.
* `test_run_multiple.cpp` - This runs the algo till no more edges can be collapsed.

Compile the files using:
```
g++ -O3 -march=native -std=c++17 -Wall -o DirEdgeCollapser <file>.cpp
```

Run the compiled file as:
```
./DirEdgeCollapser <input>.flag <output>.flag
```
