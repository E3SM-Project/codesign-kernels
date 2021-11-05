The code is intended to be used to experiment with nested loop structures to investigate
performance on CPUs and GPUs.

To build a version, you simply type
make target 
where target is a particular supported set of compiler options. For 
example, on Summit we support the targets nvhpc-cpu and nvhpc-acc 
which are the Nvidia fortran compilers and the cpu refers to a cpu-only
build while the acc version supports OpenACC for gpu offloading.
Similarly, there are targets for Cray wrappers, xlf and OpenMP Offloading.
You can look within the Makefile to see other supported targets and
add additional targets for particular machines.

The resulting executable is named "nested".

When running, the executable reads various model size options from an input
namelist file called nested.nml. It includes the model sizes nCells, nEdges,
nVertLevels. For a typical MPAS mesh, nVertLevels ranges from 60-100 with 100
being more typical. Similarly on a node, there are typically a few thousand
cells and the number of Edges is ~3-4x the number of cells. The nAdv parameter
refers to the max number of cells that contribute to an edge flux for
high-order advection and is typically 10. Finally, the nIters parameter
allows the nested loop kernels to be run nIters times to give a larger
run-time for better timing results.

This kernel code is actually coded as a single-task MPI code currently,
both because many HPC systems are expecting MPI jobs and because the
plan was to add some MPI parallelism later to explore running with
various combinations of CPU and GPU. The executable must be launched
with the usual parallel launch supported on the machine (eg mpirun, srun,
jsrun).  

If you prefer to run without MPI, you can simply comment out the
MPI_Init and MPI_Finalize calls and adjust the make target accordingly
to build an MPI-free serial version.

