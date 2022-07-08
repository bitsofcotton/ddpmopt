# ddpmopt
Convert meaning-less image into context meaning-full.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Usage:
    ./tools <step-size> <block-size> <recursion-size> <meaning-less-overwrited.ppm> <in0.ppm> ... > cache.txt
    ./tools -<step-size> <block-size> <recursion-size> <meaning-less-overwrited.ppm> < cache.txt
    # block-size    < 0 for no add ranndom.
    # recursio-size < 0 for zero clear input.
