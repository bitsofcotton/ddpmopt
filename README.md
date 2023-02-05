# ddpmopt
Convert meaning-less image into context meaning-full flavour.
This isn't useful to get full context meaning-full image by once, but we can apply this for dft-ed curvatures bitmaps and so on. So in general use, we should implement some another layers based on this.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Tips on malloc options
Some of the implementation needs to run them with specifying malloc options.
(cf. &gt;&gt;&gt; on OpenBSD)
This is because we need huge number of allocations/frees to run.

Also, we need to do ulimit or edit /etc/login.conf for large malloc use cases required by larger than medium sized input.

# Tips on shared memory
If we run them with openmp, we need large shared memory size.
They usually configurable by sysctl on unix-like systems.

# Tips on recommended input condition
We recommend to use large pgm input with rough calculation or small ppm input with accurate calculation. Otherwise, we get meaning-less result.

# Usage:
    ./ddpmopt 0  <in0.ppm> ... > cache.txt
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <out0.ppm> ... < cache.txt
    ./ddpmopt -0 <out0.ppm> ... < cache.txt
    ./predg <in0.ppm> ...
    ./qredg <inout0.ppm> ...

