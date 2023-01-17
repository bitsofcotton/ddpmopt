# ddpmopt
Convert meaning-less image into context meaning-full flavour.
This isn't useful to get full context meaning-full image by once, but we can apply this for dft-ed curvatures bitmaps and so on. So in general use, we should implement some another layers based on this.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Usage:
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <out0.ppm> ... < cache.txt
    # Following needs O(mem size * plen^1.5 * log^2(mem size * plen))
    # elementary and pred function.
    # However, the modern CPU cache implementation isn't effective to
    # this source, so the calculation clock is near 
    # baseband mem clock the cache mishit random access / const.
    # Also, we need gokicheck denlarge command after them.
    ./ddpmoptp <in0.ppm> ...
    ./ddpmoptq <inout0.ppm> ...

