% NL-Bayes image denoising.

# ABOUT

* Author    : Marc Lebrun <marc.lebrun.ik@gmail.com>
* Copyright : (C) 2013 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides an implementation of the NL-Bayes image denoising.

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS. 

- Compilation. 
Automated compilation requires the make program.

- Library. 
This code requires the libpng library.

- Image format. 
Only the PNG format is supported. 
 
-------------------------------------------------------------------------
Usage:
1. Download the code package and extract it. Go to that directory. 

2. Compile the source code (on Unix/Linux/Mac OS). 
There are two ways to compile the code. 
(1) RECOMMENDED, with Open Multi-Processing multithread parallelization 
(http://openmp.org/). Roughly speaking, it accelerates the program using the 
multiple processors in the computer. Run
make OMP=1

OR
(2) If the compiler does not support OpenMp, run
make

3. Run NL_Bayes image denoising.
./NL_Bayes
You can decide to use the "Homogeneous Criteria" in the first, second or both steps, 
by using UseArea1 = true and/or UseArea2 = true.
You can also wanted to see the bias (result of the algorithm applied
to the original image), in this case use compute_bias = true.
 
Example, run
./NL_Bayes input.png 10 1 ImNoisy.png ImDenoised.png ImBasic.png ImDiff.png ImBias.png ImBiasBasic ImDiffBias.png 1 0 1

4. Results are available in the file "measures.txt".

# ABOUT THIS FILE

Copyright 2013 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
