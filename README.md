Compute/approximate the substring complexity and related measures
===============
Author: Nicola Prezza (nicola.prezza@gmail.com), Davide Cenzato

### Brief description

This software computes/approximates the measure $\delta = max_k (d_k/k)$, where $d_k$ is the number of distinct factors of length $k$ of the input string. 

The tool `delta` computes the exact measure and outputs also related statistics, such as $argmax_k (d_k/k)$, $d_k/k$, and $d_k$ for small values of $k$. Based on the algorithm described in:

*Tomasz Kociumaka, Gonzalo Navarro, Nicola Prezza, Toward a Definitive Compressibility Measure for Repetitive Sequences. IEEE Trans. Inf. Theory 69(4): 2074-2092 (2023)*

The tool `delta-stream` computes an approximation of $\delta$ using sublinear working space. Paper: work in progress.
 
### Download

To clone the repository, run:

> git clone --recursive https://github.com/regindex/substring-complexity

### Compile

You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

We use cmake to generate the Makefile. Create a build folder in the main folder:

> mkdir build

run cmake:

> cd build; cmake ..

and compile:

> make

### Run

To compute the exact measures using linear working space, run

>  delta file

To compute an approximation of delta on a stream using sublinear working space and store the sketch, run

>  delta-stream -s -o "output_path" < file

or

> some_command_generating_output | delta-stream -s -o "output_path"

To compute an approximation of delta given a stored sketch, run

>  delta-stream -d "sketch_path"

To merge two sketches and store the result, run

>  delta-stream -m "sketch_path_1" "sketch_path_2" -o "merged_sketch_path"

To compute the approximation of the delta NCD of two sketches, run

>  delta-stream -c "sketch_path_1" "sketch_path_2"
