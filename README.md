Compute the substring complexity and related measures
===============
Author: Nicola Prezza (nicola.prezza@gmail.com)

Based on the algorithm described in the paper:

Tomasz Kociumaka, Gonzalo Navarro, Nicola Prezza, Toward a Definitive Compressibility Measure for Repetitive Sequences. IEEE Trans. Inf. Theory 69(4): 2074-2092 (2023)

### Brief description

This software computes the measure $\delta = max_k (d_k/k)$, where $d_k$ is the number of distinct factors of length $k$ of the input string. The tool outputs also related statistics, such as $argmax_k (d_k/k)$, $d_k/k$, and $d_k$ for small values of $k$.

### Download

To clone the repository, run:

> git clone https://github.com/nicolaprezza/substring-complexity

### Compile

You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

We use cmake to generate the Makefile. Create a build folder in the main folder:

> mkdir build

run cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  delta file

This command will output the statistics on the input string contained in the input file.