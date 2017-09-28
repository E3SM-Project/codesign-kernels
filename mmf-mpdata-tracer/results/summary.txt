## Results

Platform: Titan
Date: 09/28/2017.

These are times obtained for serial runs on Titan compute node using various compilers. Only, pgiacc version used the GPU.
Run using: aprun ./advect

advect.gnu.4.9.3
 CPU Timing:    1.7239000000000001E-002
advect.gnu.7.1.0
 CPU Timing:    2.2919999999999999E-002
advect.intel.17.0.0.098
 CPU Timing:    5.903000000000000E-003   <- Intel seems better for serial CPU
advect.pgi.17.7
 CPU Timing:    6.9068999999999997E-003
advect.pgiacc.17.7
 CPU Timing:    7.0648000000000004E-003
 OpenACC-1 Timing:    4.8589999999999999E-004
 OpenACC-2 Timing:    3.0899999999999998E-004
 OpenACC-1 Timing:    3.9070000000000001E-004
 OpenACC-2 Timing:    2.9399999999999999E-004
 
 A detailed nvprof report is present in the directory.
