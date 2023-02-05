# ddpmopt
Convert meaning-less image into context meaning-full flavour.
This isn't useful to get full context meaning-full image by once, but we can apply this for dft-ed curvatures bitmaps and so on. So in general use, we should implement some another layers based on this.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Tips on calculation order
ddpmoptp and ddpmoptq needs O((input mem size) * plen^1.5 * log^2((input mem size) * plen)).

# Tips on malloc options
Some of the implementation needs to run them with specifying malloc options.
(cf. &gt;&gt;&gt; on OpenBSD)
This is because we need huge number of allocations/frees to run.

Also, we need to do ulimit or edit /etc/login.conf for large malloc use cases required by larger than medium sized input.

# Tips on shared memory
If we run them with openmp, we need large shared memory size.
They usually configurable by sysctl on unix-like systems.

# Tips on accuracy
From some numerical test, 32 bit is a slight rough to ddpmopt when input number is small. 64 bit is better but slightly not enough.
From ddpmopt meaning, ddpmoptp should use condorcet jury condition with ddpmoptpm.
The modern middle range laptop computer should solve 64x64 multiple color image with 64bit condorcet jury ddpmoptpm in several minutes with OpenMP, but they are rough enough because PRNG with optimization is not better combination to solve some images, so we should shrink output by half or some.

# No use method
We work with small enough image better, but if we gather them and revert to large image, it's worse than original image prediction.
Either, if we work with imagemagick's resize, despeckle, resize chain, they make not better than original material by recognized by human (this is not tested with quantity).

# Tips on recommended input condition
We recommend to use large pgm input with rough calculation or small ppm input with accurate calculation. Otherwise, we get meaning-less result.

# Usage:
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <out0.ppm> ... < cache.txt
    ./ddpmoptp <in0.ppm> ...
    ./ddpmoptq <inout0.ppm> ...
    ./predg <in0.ppm> ...

