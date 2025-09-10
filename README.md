# Directed Graph Edge Collapser

This contains 2 header files:
* `DIR_EDGE_COLLAPSER.h` - Main version of the algorithm
* `DIR_EDGE_COLLAPSER_MEM.h` - Version 2 which uses bitset, need to specify the maximum number of vertices in the line `constexpr size_t MAX_VERTICES = 700;`. Highly memory efficient, mostly used for high density graphs. Takes more time than v1.

The main files:
* `diredgecollapser.cpp` - This runs the algo once, displaying the memory usage statistics.
* `diredgecollapser_mult.cpp` - This runs the algo till no more edges can be collapsed.
* `diredgecollapser_mem_mult.cpp` - Multiple iterations with memory optimization.

Compile the files using:
```
g++ -O3 -march=native -std=c++17 -Wall -o DirEdgeCollapser <file>.cpp
```

Run the compiled file as:
```
./DirEdgeCollapser <input>.flag <output>.flag
```
