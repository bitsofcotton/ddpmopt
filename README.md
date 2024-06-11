# ddpmopt
Apply some of the filter to input stream.
We can use this for bitsofcotton/i2g filtered images.

This should behaves deterministic ones, however, these are very sensitive on accuracy, so only lowest bit has a difference condition, we get whole image different condition.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them but different flavoured one, we only focus to apply each pixel context to color image into monochrome one, which have the structure completely depends on filters' multiple meaning or complexity.

# Tips on malloc options
We need to do ulimit or edit /etc/login.conf for large malloc use cases required by larger than medium sized input.

Using this with mimalloc or so can increase memory usage with multi thread on some systems.

# Tips around c++ compilers
Some of the lieonn.hh operator \>\> class doesn't work as expected, might be compilers' bug.

# Tips around input sizes on this
We get function entropy depend result with [pq]redg, this is because we treat each pixel as equally treated ones without any of the internal states.
So we should adjust in/output by them. eg. goki_check_cc:clean\[lcC\] before and after doing them.
Around better size of this, please refer bitsofcotton/p8 for extends.

However, cleanl input and despeckle or negate output also could work.

# P01 hypothesis
The P01 predictor makes the hypothesis the structure is continuous enough.

In the case there's brand new observation on the states on each pixel/image context should have the next image condition, we fail with this predictor.
This is specification of this implementation, so is intended to be so.

N.B.
However, the purpose what is the brand new observation on the stream is obscure, so it's some of the quantity counting on continuous structure.
So to predict, separating meaning into many much of the dimension causes {0,...,0} cannot grow some another values, so it's not a better tactics for us.
So blending some of the meaning into some dictionary based or statistics based representation dimensions and predict with them, then, re-separate them is better for us. However, the blended image nor separated meaning is the matter if such of occasions work with to be correctly working.
These predg/qredg has a context of each pixel/each image with function entropy depend size image, so to improve this, we need some of the context or some of the dictionaries or some of the pattern to complement with.
The randtools compression/decompression depends on factor matrix sparsity, however, they depends on observation/heat could be after of the calculation, so we should have some of the dictionary based functions, however, they're already done by deep learning recent days.

So tips on predicting graphics set, if there's some of the continuity we can recognize, we get better result for us however, if there's no continuity we can recognize, we will get some continuous result on the place we cannot recognize, so to recognize the one is harder and harder than former ones.

So around them, we can re-emphasis original meanings they satisfies and enlarge them by ongoing deep learning methods from some of the groups on inet s.t. GPTn-based graphics translators. So we recommend to use them s.t. img2img on them after predicting graphics set by ours.

However, our predictor is being jammed condition now, even the base system/compiler is correct case, there's some of the saturation on function pattern max number concern switch by intent can occur.
This is because we test P01 predictor with bitsofcotton/p1:p1.cc with -28, 28 argv with bitsofcotton/p2:cr.py with s for each by raw .cc file and predict some of the PRNG on this machine, even the first argument often goes well if another predictors are filtering condition, output can be controlled as to be return to origin point condition, so p1 28 \| p1 -28 \| cr.py s pseudo condition after applying the raw .cc filter can work well on my machie before this post.
So after posting this, the condition can be broken from our experiences.

So if this is the jamming from clang++ nor base system nor function saturation to P01 predictor, the result can be broken any which in the meaning on the prediction. We observe this especially some graphics set prediction with only the first predicted graphics going to be most of the null condition graphics. We didn't discriminate whether this is from the graphics theirselves condition with our predictor nor, it's from the intension from the jammer intent selves.

# Usage:
    ./predg(32|64)?(mp)? <in0.ppm> ...
    ./qredg(32|64)?(mp)? <in0out.ppm> ...
    ./ddpmopt(32|64)?(mp)? + <in0in.ppm> <in0out.pgm> ... > cache.txt
    ./ddpmopt(32|64)?(mp)? - <in0.ppm> ... < cache.txt
    ./tcont [xyit] <in0.ppm> ...
    cp `./tcont i <in0.ppm> ... | sort | head -n ... | tr '\n' ' '` outdir

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
2024/04/17 update readme.
2024/04/18 add tcont.cc, real close with omake.
2024/04/18 won't update without lieonn.hh change.
2024/04/29 update readme.
2024/05/05 p01 class fix step-step to 1-step correction.
2024/05/07 correct output depth limit to have a meaning. update readme.
2024/06/01 add predg0... qredg0... compile option, these using latest results on \{p0,p1,p2\}.
2024/06/02 revert. it's nonsense.
2024/06/05 fix p01 crash in rare cases.
2024/06/07 fix number of predictions to reasonable one. add another implementation on python predg.cc, only QR decomposition is differ but this has a better results?? update readme.
2024/06/09 factorize into each bit and predict with them. leave with this but this have color intensity == {0,1} confusion bug.
2024/06/09 fix last bug. average step skips. add readme.
2024/06/11 update readme.
2024/06/12 update readme.

