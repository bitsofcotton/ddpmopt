# ddpmopt
Convert meaning-less image into context meaning-full flavour.
This isn't useful to get full context meaning-full image by once, but we can apply this for dft-ed curvatures bitmaps and so on. So in general use, we should implement some another layers based on this.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Usage:
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <out0.ppm> ... < cache.txt

