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
We take class Prange prediction on each pixel of makeProgramInvariant.
This is because we make the hypothesis there exists hidden invariant on whole of the set of images for each image in complexity-wise meaning.

# Tips on context structure we don't suppose
We don't suppose the structure on the inputs that effects the whole program differential equation coefficient or its distribution who describes the structure.
So we only make the hypothesis above, we should shrink the images after \[pq\]redg or ddpmopt. This is because only the each pixel context on whole invariant is determined.

# Tips on goki check cc denlarge penlarge chain
Some chain causes non expected result they seems to be not broken in meaning.
So we should use per one task per one method.
Also, we should use huge number of the input datas because of non context-full but calculation over-{} return result gained otherwise.
However, we recommend you and us to use denlarge sharpen chain on goki_check_cc after predg/qredg.

# Usage:
    ./ddpmopt 0  <in0.ppm> ... > cache.txt
    ./ddpmopt 00 <in0.ppm> ... > cache.txt
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <inout0.ppm> ... < cache.txt
    ./ddpmopt -0 <inout0.ppm> ... < cache.txt
    ./predg <in0.ppm> ...
    ./qredg <inout0.ppm> ...

# Real close
2023/03/01
2023/03/09 bug fix after close #1.
2023/03/13 bug fix after close #2.
2023/03/13 integrate all files into lieonn.hh after close #3.
2023/03/18 merge latest p0, after close #4.
2023/03/20 qred makes only row-direction prediction.

