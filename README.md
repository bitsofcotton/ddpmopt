# ddpmopt
Apply some of the filter to input stream.
We can use this for bitsofcotton/i2g filtered images.

This should behaves deterministic ones, however, some of the condition causes non deterministic flavours even ever in -D_FLOAT_BITS_=... condition.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them but different flavoured one, we only focus to apply each pixel context to color image into monochrome one, which have the structure completely depends on filters' multiple meaning or complexity.

# Tips on malloc options
We need to do ulimit or edit /etc/login.conf for large malloc use cases required by larger than medium sized input.

Using this with mimalloc or so can increase memory usage with multi thread on some systems.

# Tips around c++ compilers
Some of the lieonn.hh operator \>\> class doesn't works as expected, might be compilers' bug.

# Tips around input sizes on this
We get function entropy depend result with [pq]redg, this is because we treat each pixel as equally treated ones.
So we should adjust in/output by them. eg. goki_check_cc:cleans\[lc\]\? before and after doing them.
Around better size of this, please refer bitsofcotton/p8 for extends.

# When we have broken result on some of the [pq]redg output
The P01 predictor makes the hypothesis the structure is continuous enough and timing aligned input.
If in the best case, it's around 2/3, if in the worst case, it's around 1/3, so sometimes negated image also works well. However, in the worst on prediction case, it's 1/2, we cannot predict this with such of the timing value.
So in the case some of the gulf appears to the next step or in the case some timing attack counter measures appears on input stream, some step previous/next seems to break.
This is specification of this implementation, so is intended to be so.

We can specify skip number to [pq]redg.

# Usage:
    ./predg(32|64)?(mp)? <in0.ppm> ...
    ./qredg(32|64)?(mp)? <in0.ppm> ...
    ./ddpmopt(32|64)?(mp)? + <in0in.ppm> <in0out.pgm> ... > cache.txt
    ./ddpmopt(32|64)?(mp)? - <in0.ppm> ... < cache.txt

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
2023/06/13 update readme. ddpmopt.cc large change.
2023/06/14 ddpmopt.cc fix now works, update readme.
2023/06/15 update readme, some fixes around predg, qredg, ddpmopt.
2023/06/24 fix qredg.cc autoLevel in lieonn.
2023/06/27 ddpmopt.cc sz 3 to 2 change, update readme.
2023/07/07 [pq]redg predict along with middle image/line.
2023/07/08 invariant causes +1. qredg.cc fix crash.
2023/07/12 update readme.
2023/08/02 update topt, lieonn as crush to the last.
2023/08/03 topt, ddpmopt fix on apply. update readme.
2023/08/04 update readme. fix predv result norm.
2023/08/26 ddpmopt - option calculation change, update readme. update [pq]redg for recursive ones to beat with geometric average limit.
2023/09/04 auto configure predg/qredg param, conservative.
2023/09/11 auto configure predg/qredg size, aggressive ga param.
2023/09/12 fix last broken predv func and qredg.cc.
2023/09/13 change ddpmopt retry and retry on geometric average.
2023/09/19 autogamma after doing predg. update readme.
2023/09/19 predvResizeMat resize size fix to most reasonable one.
2023/09/22 ddpmopt change not to use crush but with linearInvariant.
2023/09/24 fix last up, don't know why they worked well without crash on last debug.
2023/09/25 change output size strategy, not using resize, preferring complement to predict.
2023/10/03 update readme.
2023/10/05 update readme, should close except for some of the ddpmopt for pairs of the images.
2023/10/18 update readme.
2023/10/19 mipmap impl, update readme.
2023/10/20 revert, mipmap doesn't work well.
2023/10/22 instead of complement trace of input, we should do shrink after output is done. This is from some of the numerical tests, so whole image orthogonality isn't matter even when input data is small enough, instead of them, we try to vanish some prediction errors with geometric mean.
2023/10/23 ddpmopt strategy algorithm in/output large change.
2023/10/26 update readme.
2023/10/27 update readme. close.
2023/10/30 copy structure reliably with randtools meaning.
2023/11/13 update readme.
2023/11/19 revert each line complement condition, force to use complement.
2023/11/20 only complement with 2, they smoothes output enough with our measure (&gt; 2/3).
2023/11/21 update readme.
2023/11/23 fix known tips, there was at most double error.
2023/11/25 implement and revert the test to complement before to revertProgramInvariant, they doesn't improve well differed to shrinking.
2023/12/06 fix ddpmopt.cc as makeProgramInvariant, revertProgramInvariant to better compatible with [0,1].
2023/12/07 fix ddpmopt.cc makeProgramInvariant double apply serious bug in crush. Also, update readme as compatible with in/output. realclose.
2023/12/15 use the tactics not to apply twice make/revertProgramInvariant on prediction, the invariant is already taken, however, this can causes P^t 1 == 0 condition on linearInvariant, we don't fix them.
2024/03/19 half p8 compatible change on prediction P1I to P01 with predg... qredg... binary.
2024/04/02 fix P01, vanish predvc.
2024/04/04 only use large accuracy on calculating pnextcache, but this is broken with cache naming.
2024/04/07 rgb2xyz, xyz2rgb on first/last of predg/qredg if color is 3ch.
2024/04/09 fix maeProgramInvariant x_k==0 to wrap into x_k:=1, refix pnext with lower digits.
2024/04/10 we count the function entropy enough beat with in/output on [pq]redg.
2024/04/12 update readme.
2024/04/14 take a median after predv before revertProgramInvariant.
2024/04/15 update readme, it is the specification of this some output to be broken.

