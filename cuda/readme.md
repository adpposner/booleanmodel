# Companion simulation model for solution verification

The current directory contains a GPU-based (specifically CUDA) simulation to test some of the results from the solution model. It has been tested to work with NVIDIA's CUDA Toolkit v10, built using arch=sm_61 and tested on both and NVIDIA GTX 1060 and NVIDIA GTX 1080ti. 
The CUDA toolkit is available [here](https://developer.nvidia.com/cuda-downloads); one must obviously have a compatible graphics card as well (configured in `main.cu` as device 0 - change to whatever device is appropriate).

Building this project should not take much more than a simple `make` command, given that all of its dependencies (e.g. CURAND) are part of the CUDA toolkit.

**Note** - make sure to check and set the macros in `util.h` appropriately, specifically `NVPO` and `NUM_NETS_PER_PARAMSET`. The former specifies the amount of testing of each network for a given parameter set, while the latter corresponds to the *number* of generated/tested networks. 
Higer values will tend to produce higher precision results at the expense of compute time.

## Usage
The system is meant to be used to reproduce results from the direct solution, and is configured accordingly. 
The program is invoked using:
> `./main path_to_input output_prefix num_paramsets`
Here, `path_to_input` should be a valid path to an *output* file from the direct solution - in the parent directory, `sample_output.txt` can be used as a valid input. 
`output_prefix` will designate the path/prefix of the output files, of which there will be two. 
The first, `[output-prefix].params` should just contain a list of the parameter sets and values of nMicro, analogous to the first 11 columns of the input file.
The second, `[output-prefix].matrices` contains flattened average transition matrices. Each row corresponds to a full matrix.
The last value `num_paramsets` will tell the program how many lines of the input file to simulate (or until EOF)

This format is clearly different from the original direct solution file format. Invoking many of the linear algebra subroutines and such are inconvenient in the C programming language, to say the least. 
Because this model is used more for verification than for obtaining results, the linear algebra is relegated to a simple python script, which is well-suited to handle the relatively low throughput of the simulation-based model.

One can simply enter (from the working directory):
> `./simresults.py output-prefix`
using the value of `output-prefix` from before. 
This will produce a file with the extension `[output-prefix].processed.txt`, which is directly comparable (i.e. same columns) to the original output data tested.

### Caveats
GPUs are great, but getting the same kind of resolution in transition matrices corresponding to the direct solution is just not possible in a reasonable amount of time. While this program can certainly be optimized plenty, the smallest nonzero entry in any transition matrix obtained from the simulation is 1/(NVPO*NUM_NETS_PER_PARAMSET), obviously not the most fine-grained set of numbers. 
Furthermore, though I haven't confirmed it, the commodity GPUs only really work with single-precision, and it is very possible that the pseudorandomn number generators from the CURAND library may not produce sufficiently uniform floats to really get at "the answer." 
Nonetheless, it still does back up the solution better than expected. Any suggestions are welcome, as I am not a CUDA guru of any sort.

### Questions
Some amount of computer programming knowledge (as well as appropriate hardware) is required in order to compose and test the simulation model, although it should be easier than the direct solution if one has the CUDA toolkit. 
If you hit any roadblocks, contact me and I will help to the best of my ability. Thanks!
