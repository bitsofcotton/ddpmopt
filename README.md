# ddpmopt
Convert meaning-less image into context meaning-full flavour.
This isn't useful to get full context meaning-full image by once, but we can apply this for dft-ed curvatures bitmaps and so on. So in general use, we should implement some another layers based on this.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Tips on calculation order
ddpmoptp and ddpmoptq needs O((input mem size) * plen^1.5 * log^2((input mem size) * plen)).
We can 1line when we met ddpmoptp.
Either of them needs denlarge in gokicheck after calculation.

# Tips on malloc options
Some of the implementation needs to run them with specifying malloc options.
(cf. &gt;&gt;&gt; on OpenBSD)
This is because we need huge number of allocations/frees to run.

# Usage:
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <out0.ppm> ... < cache.txt
    ./ddpmoptp <in0.ppm> ...
    ./ddpmoptq <inout0.ppm> ...
    ./1line 0 <inout0.ppm> ...
    ./1line -<number of colums> <inout0.ppm> ...

