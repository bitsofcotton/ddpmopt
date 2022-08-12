# ddpmopt
Convert meaning-less image into context meaning-full flavour.
This isn't useful to get full context meaning-full image by once, but we can apply this for dft-ed curvatures bitmaps and so on. So in general use, we should implement some another layers based on this.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Usage:
    ./tools <step-size> <block-size> <recursion-size> <meaning-less-overwrited.ppm> <in0.ppm> ... > cache.txt
    ./tools -<step-size> <block-size> <recursion-size> <meaning-less-overwrited.ppm> < cache.txt
    # block-size    < 0 for no add ranndom.
    # recursio-size < 0 for zero clear input.
