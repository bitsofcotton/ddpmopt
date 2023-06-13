# ddpmopt
Convert meaning-less image into context meaning-full flavour.
This isn't useful to get full context meaning-full image by once, but we can apply this for dft-ed curvatures bitmaps and so on. So in general use, we should implement some another layers based on this.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them.

# Tips on malloc options
Some of the implementation needs to run them with specifying malloc options.
(cf. jjj&gt;&gt;&gt; on OpenBSD)
This is because we need huge number of allocations/frees to run.

Also, we need to do ulimit or edit /etc/login.conf for large malloc use cases required by larger than medium sized input.

Also, if we run these programs with openmp, we need large shared memory size.
They are usually configurable by sysctl on unix-like systems.
However, openmp costs large enough beat with single core even in multicore because of their (and our) memory architecture.

# Denoise
We recommend you and us to use goki denlarge+? 1 ..., this is because to denoise analytic noises we cannot reduce by calculate them once.
Then, we recoomend one of the two methods after denlarge+?: (i) mogrify -despeckle -despeckle -equalize ... or (ii) use some of the enlarge method alike stable diffusion.

# The ddpmopt binary glitch
We need rewrite to beat with ddpmopt binary as up to 6 or 11 elements they represents per each stage caused by s2m2s condition.
We need 3x3 tile or so to represent them.

# Usage:
    ./ddpmopt 0  <in0.ppm> ... > cache.txt
    ./ddpmopt 00 <in0.ppm> ... > cache.txt
    ./ddpmopt +  <in0.ppm> ... > cache.txt
    ./ddpmopt ++ <in0.ppm> ... > cache.txt
    ./ddpmopt -  <inout0.ppm> ... < cache.txt
    ./ddpmopt -0 <inout0.ppm> ... < cache.txt
    ./topt <dim>? < in.txt > cache.txt
    cat cache.txt - | ./topt -<num>
    ./predg <in0.ppm> ...
    ./qredg <inout0.ppm> ...

# Real close
2023/03/01
2023/03/09 bug fix after close #1.
2023/03/13 bug fix after close #2.
2023/03/13 integrate all files into lieonn.hh after close #3.
2023/03/18 merge latest p0, after close #4.
2023/03/20 qred makes only row-direction prediction, after close #5.
2023/03/24 ddpmopt.cc +0 00 -0 command fix, code clean, after close #6.
2023/03/29 merge p0 pnext, after close #7.
2023/03/31 persistent prediction to get average 1/2 * 2/3 pixels.
2023/04/02 merge p2 result.
2023/04/03 better simple predv.
2023/04/04 update readme.
2023/04/05 fix makeProgramInvariant scale accuracy stability.
2023/04/19 add topt.cc.
2023/04/21 shape up around makeProgramInvariant/revertProgramInvariant, algorithm changed.
2023/04/23 qredg.cc prediction/original norm fix.
2023/05/18 predv function change to better analytical, update README.md.
2023/06/08 predv normalization fix.
2023/06/11 topt.cc output normalization fix. update readme.
2023/06/13 update readme.

