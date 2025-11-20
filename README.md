# ddpmopt (closed)
Apply some of the filter or predict input stream.
We can use this for bitsofcotton/i2g filtered images.

# Context
There exists Denoising Diffusion Probabilistic Models (DDPM; Ho et al. 2020). So this is another try on them but different flavoured one, we only focus to apply each pixel context to color image into monochrome one, which have the structure completely depends on filters' multiple meaning or complexity.

# Tips on malloc options
We need to do ulimit or edit /etc/login.conf for large malloc use cases required by larger than medium sized input.

Using this with mimalloc or so can increase memory usage with multi thread on some systems.

# Tips on predictors
Implanted comments into lieonn.hh .

# Usage:
    # copy color structure
    ./ddpmoptp?(mp)? + <in0out.pgm> <in0in.ppm> ... > cache.txt
    # apply color structure
    ./ddpmoptp?(mp)? - <in0.ppm> ... < cache.txt
    # predict following image
    ./ddpmoptp?(mp)? p <in0.ppm> ...
    # reverse whole pixel context (each bit input)
    ./ddpmoptp?(mp)? w <in0.ppm> <in0.ppm-4.ppm> ... <addition-4.ppm>
    # predict down scanlines
    ./ddpmoptp?(mp)? q <in0out.ppm> ...
    # show continuity
    ./ddpmoptp?(mp)? [xyit] <in0.ppm> ...
    # some of the volume curvature like transform
    ./ddpmoptp?(mp)? c <in0.ppm> ...
    # test input series of graphics predictable or not into test.ppm
    ./ddpmoptp?(mp)? T <in0.ppm> ...

# Leave
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
2024/06/13 update readme.
2024/06/14 update readme.
2024/06/15 conclude 2024/06/11-2024/06/14 conditions readme.
2024/06/16 revert P210 to original, then, P01, P0 pred temporarily, update readme.
2024/06/17 merge p2 logic with p10 class.
2024/06/18 code cleaning, update readme.
2024/06/18 speed remedy.
2024/06/19 add restriction for getting average on PprogressionOnce010n.
2024/06/20 our machine is infected, take a most logically valid in our program predictions.
2024/06/21 fix fatal error on PprogressionOnce::next, they doesn't use predictors.
2024/06/21 revert and brush up, add fiocursed.cc series, brush readme.
2024/06/22 update readme, updates around pprogression causes only 1 step after/before prediction.
2024/06/22 p01 fatal fix. make/revert program invariant change friendly to predictions.
2024/06/23 large change around class instance initializer, also have progression short range fix.
2024/06/23 fatal fix around last update.
2024/06/24 fatal fix around last update, rotate predMat sloppy case.
2024/06/26 some of the assertion fix, update readme.
2024/06/26 update readme, fix around Ppersistent buf.full condition.
2024/06/27 fix predv last norm condition calculations. predMat bugs might be fixed.
2024/06/29 update readme and comments.
2024/06/30 re-insert periods with better stable method. update readme.
2024/07/06 Ppersistent now use maximum length for predictions. Also readme update.
2024/07/07 code cleaning. merge Pprogression improve but no affects.
2024/07/08 internal state range strategy change, use all of the input to reduce. update readme.
2024/07/09 revert bitwise prediction causes whole image invariant works same as theoretical ones, however, each pixel context isn't enough on prediction but is enough on whole image condition information amount as better weighted. update readme.
2024/07/10 revert [pq]redg.cc as no each bit condition, instead of this, use goki_check_cc:test.py:bit command.
2024/07/20 update readme, might our system is infected.
2024/08/18 update -\[369\] predictors for recursive but equivalent.
2024/09/03 update \[pq\]redg...p.. for auto tuned recursive but for tiny images.
2024/09/04 update last up with proper recursive value.
2024/09/05 omit error output in zeroFix.
2024/09/06 update and fix readme.
2024/09/09 merge p1/pp3.cc result, change only output forward pred ones.
2024/09/10 merge p1/pp3.cc result, re-re-leave.
2024/09/10 fix pnoise meaning. update readme, re-re-re-leave.
2024/09/12 update readme. leave.
2024/09/22 append dft hack, add readme, releave.
2024/09/23 fix _PREDV_==3 _PREDV_DFT_ case unit value. Changed preddg... to predfg... to make some readability. Fix around comments and readme.md.
2024/09/24 brush up, eliminate exhaust of the resource to get tiny output improve in finite and up to aleph_0 condition.
2024/09/25 elim dead code, update readme. leave.
2024/09/26 improve heap resource efficiency.
2024/09/27 refactoring predv, predvp0. the target stream we predict is now concrete, so PP0 type changed.
2024/09/28 add step after predictions.
2024/09/29 update readme. leave the repository really.
2024/10/05 add predga.cc and compile option.
2024/10/26 fix predvall meaning after p2/README.md:Tips on reseed.
2024/10/31 add predgs.cc.
2024/11/01 fix last predgs.cc index map, update readme, accuracy is not enough for SVD.
2024/11/02 update readme. elim predgs.cc, predga.cc.
2024/11/03 update readme. rerere-leave here. predg.cc a cmd list up fix.
2024/11/12 delete tips on reseeding, reseeding is not so harder. replaced flip, flop template function in lieonn suitable with gcc however 128bit long double isn't compile on our main pc environment.
2024/11/14 integreate all commands on this repository into ddpmopt.cc but the binary is very fat.
2024/11/17 add w command.
2024/11/19 improved lieonn.hh:taylor command speed and accuracy this causes q command better works. update readme. something error occured first upload of this change on github.com. this change leads us to pnext r variable doubles.
2024/11/20 update readme for recent knowns.
2024/11/20 update readme.
2024/11/30 add c command. update readme.
2024/12/02 taylor improvement, taylor function reclose with this.
2024/12/03 w command fix also readme.md update.
2024/12/04 fix w command output, backward had a glitch, so eliminated. update readme.
2024/12/05 backport p1 | p0 results, brush up code, replace [0a] command to p command, update readme.
2024/12/07 c command fix.
2024/12/09 changed to output only a single prediction. w command crash fix, memory efficiency improve. q command crash fix.
2024/12/11 fix readme w command usage.
2024/12/13 leave here, might return here.
2024/12/26 use montecarlo method instead of doing each step average on prediction.
2024/12/28 fix 'q', 'w' commands with last change.
2024/12/29 update 'w' commands suitable with predv1 method impementation is predv4.
2024/12/30 update readme.
2025/01/06 eliminate condorcet's jury method, they've no effects.
2025/01/06 update readme for compatible with latest goki_check_cc.
2025/01/27 improve pred... memory usage without predv4.
2025/02/01 fix readme memory usage notation.
2025/02/05 predv function to get better prediction - real value distribution by PRNG tests.
2025/02/15 add PP0 as PSVD ... as a dead code, they doesn't improve output enough on our machines with small number of inputs.
2025/02/16 fix predv4 alignments affects all of outputs w command.
2025/02/18 revert using predvp0 to using predv, they might come from infection.
2025/02/20 move include comments into lieonn.hh . update reamde.md fix meaning on predictions we will re freeze with this.
2025/02/22 not optimal but better looking q command output size with specifying step to predictor.
2025/02/23 add readme.md notes.
2025/03/01 add readme.md note around DFT.
2025/03/03 add T command for test. revert subtraction to multiplication and sgn method to have gokibin bit preprocessed inputs.
2025/03/04 apply T command tests into original p, w, q command. either revert to original p, w, q commands with renewing T command test.
2025/03/05 our invariant condition is being attacked, we use 2 of dimension output but in fact we need at least 4 dimension output for all.
2025/03/06 yellow output is lead by small input number, also some readme fix we often don't need entropy feeding control.
2025/03/09 brush up lieonn.hh phase periodical jamming matters. we only make hypothesis PRNG we use isn't match the predictor/original stream phase period they have.
2025/03/11 add and fix readme. close.
2025/03/12 brush up readme, freeze.
2025/03/13 add PQ command, update readme around 4 of candidate results.
2025/03/22 close with this Readme.md.
2025/04/01 add readme.md because we're in infected condition, also close because of the condition.
2025/04/17 auto tune dimension in F_2 case other than 4 dimension to target.
2025/04/18 qQ command strategy change.
2025/04/19 rebrush up lieonn.hh easy to read whole, fix ind2vd lt, gt mis exchange.
2025/05/16 backport p2 result causes single output per each.
2025/05/20 we eliminate PQ command delta in/output because of backporting p0/p0p.cc result.
2025/05/21 T command extend, update readme, reclose.
2025/05/21 slim down, w command out change.
2025/05/23 retarget cultivated input stream but will close soon.
2025/05/25 code cleaning, select FEED_MUCH as a single implementation default.
2025/06/08 rework into possible thin layer but enough layers from p2 result.
2025/06/10 persistent uint32_t use option with _PERSISTENT_ compile option but slightly use int32_t for index op and void* for pointer.
2025/06/11 compat compile option to gcc4.2.1.
2025/06/12 compat compile option with one variant of gcc2.95.3.
2025/06/17 fix deep template function reverse computation. update readme. close.
2025/06/19 fix deep reverse computation with logical one, we don't trust numerical test on this machine. close.
2025/06/20 merge latest p2 result includes operator >> on simple float accuracy fix around PERSISTENT option.
2025/06/21 add README.md target result section.
2025/06/22 change T command output to better reasonable one. output 2 of the image because we target binary valued result.
2025/06/22 merge concept ok. leave.
2025/06/23 code clean, flush. update readme, comment, a little speed remedy.
2025/06/25 readme.md move into lieonn.hh comment implanted. also implement some stopping layeres.
2025/06/26 add T- command, update lieonn.hh comment.
2025/06/28 refactor/fix around lieonn. re-compat with gcc2953.
2025/06/29-30 refactor/investigate/fix around lieonn. add pAbsentMajority. T,p,q command chg.
2025/07/01 various bug fixes, some speed remedy, p012next strategy change.
2025/07/02 add upper counter measure for jammers however they're extremely heavy.
2025/07/02-03 refresh vs. jammer conditions, we add new pFeedLebesgue function they caused us better structure. debug ok.
2025/07/04 speed remedy, debug ok, comment diet.
2025/07/06 brush up lieonn. bug fixes. revertByProgramInvariant important fix.
2025/07/10 add pPersistentQ to shift gulfs. code cleaning.
2025/07/12-13 persistent debug, normalize pCbrtMarkov input, lieonn.hh refactoring.
2025/07/14-16 some addition to measureament condition, debug and slim up comments in lieonn.
2025/07/17-19 blending PRNG, param meaning change, T command shrink down.
2025/07/20 brush up, debug, also comment on lieonn.hh. update readme.
2025/07/24 stacking bricks causes this result, parameter auto configure, brush up.
2025/07/25 measureable condition ok. exclude _P_JAM_ things to compile option.
2025/07/26-27 levi stream condition for abstract on lieonn.hh. we should have embryonic condition close later on here.
2025/07/28 update readme, pPolish now don't eliminates unstable region.
2025/08/01 add upper layers works might be well.
2025/08/02-03 sectional improvement.
2025/08/04-07 merge latest p2 result.
2025/08/11 merge latest p2 result.
2025/08/12-15 merge latest p2 result.
2025/08/16 merge latest p2 result.
2025/08/16 merge latest p2 result.
2025/08/17-23 merge latest p2 result either our computer is infected (but this is high probability).
2025/08/25 merge latest p2 result as combine candidates.
2025/08/30 merge latest p2 result as separating again but append measureable cond.
2025/09/01 offset bug fix on pAppendMeasure.
2025/09/05 fix to compile with latest g++, pAppendMeasure fix (should be done).
2025/09/16 w cmd order change.
2025/09/25 fix p01next, pRS bugs, don't know why they worked well.
2025/10/03 should be compiler's bug. -&gt; something index matter, ok.
2025/10/04 merge latest p2 result causes ok result.
2025/11/03 check point.
2025/11/12 fix last up by p2 logic.
2025/11/12 better result reducing harmful layers.
2025/11/15 better thin calculation layers. code clean.
2025/11/17 merge latest p2 result.
2025/11/21 merge latest p2 result.

