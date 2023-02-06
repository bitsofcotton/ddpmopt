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

Also, if we run these programs with openmp, we need large shared memory size.
They are usually configurable by sysctl on unix-like systems.

# Tips on small images
If we treat small images instead of the originals, they seems to be better by looking, but if we treat original images as some slided parts of them and invert originals, they doesn't enhance at all in these programs.

# Tips on predg/qredg
We take class P prediction on each pixel of makeProgramInvariant.
This is because we make the hypothesis there exists hidden invariant on whole of the set of images in complexity-wise meaning.

# Usage:
    ./ddpmopt 0  <in0.ppm> ... > cache.txt
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <out0.ppm> ... < cache.txt
    ./predg <in0.ppm> ...
    ./qredg <inout0.ppm> ...

