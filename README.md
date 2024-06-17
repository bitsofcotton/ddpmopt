# bitsofcotton/p0
Generic predictor for Riemann measurable input streams and their variants (on aleph_1 cf. C).

There's a plenty of the room to make this into n-variable predictor, but this repository won't implement such of them.

# Contexts
There exists discrete fourier transform on given (same interval) series (This is well described on everywhere.).
And, if we make DFT and IDFT on them, there exists differential on them in DFT meaning.
(In continuous meaning, this is described a little on the books I refered,
 but, in discrete meaning, I can't find preceding results but might be exists.)

# How to use:
    P0maxRank<double> p(status);
    // pnext<T>(range, step) for riemann measurable prediction vector.
    ...
      xnext = p.next(x);

# How to use (commandline):
    ./p0 <status>? <progression>? < stream.txt
    # 0 > status for blurred next step prediction, needs 2 <= \|status\|.
    # 0 < status for P0maxRank plain prediction.
    # 0 > progression for next 1 step prediction with some of the deltas.
    # 0 < progression for next |progression| step prediction.

# Check status length is valid for accuracy or not
    ./p0r <len> <step> <r>
    # last line shows accuracy is valid or not. if valid, it's near \[1, 0, ...\].

# Proof on Riemann measurable condition.
If original function is in C1, there exists F(z,&theta;) := f(z+\bar{z})+i\*f(z-\bar{z})\*tan(&theta;) in C1 on z in C.
So each &theta;, exists F: holomorphic function at some axis that real axis is same as f.
So &theta; = &pi;/4 \pm &epsilon;, if we have a area, with gutzmer's inequation and edit integrate path,
f have laurent series and as the upper bound of coefficients,
we can cut them as error in numerical calculation on some cut off.
If f has a laurent series around {x| |x-a|<1} and no essential singular point on there,
taylor series description exists in weak differential meaning because 1+x+x^2+...
converges and x^k -&gt; 0 and k -&gt; \infty on there.  
There's possibility to have area with lagrange multipliers with limit operations, but if then, it has to use the condition of f
other than the condition below C0. And if then, tan(&theta;) should not have pole or infinite value.
exp(Sum log(z-|z|cis(&pi;/4+t_k)))\*(Sum((f(z+\bar{z})+i\*f(z-\bar{z})\*tan(&pi;/4+t\_k))/(z-|z|cis(&pi;/4+t\_k)) =: f(z+\bar{z})\*g(z)+f(z-\bar{z})\*h(z) = f(z+\bar{z})\*(g\_real(z+\bar{z})+i\*g\_imag(z-\bar{z}))+f(z-\bar{z})\*(h\_real(z+\bar{z})+h\_imag(z+\bar{z})) == (f(z+\bar{z})\*g\_real(z+\bar{z}) - f(z-\bar{z})\*h\_imag) + imaginary part, so replace on (1+i)t, F(x):=(f(x)\*g\_real(x)-f(x)\*h\_imag(x))==f(x)\*(g\_real(x)-h\_imag(x)), F is holomorphic on some small range (f in C0) and Re F = f(x)\*some G(x), G(x) not depends f, so f(x) / G(x) =: f2(x), same method, f(x) will be holomorphic on some small range.

And, we can weaken this condition with cauchy's integrate theorem on ja.wikipedia.org link (doi:10.1090/S0002-9947-1900-1500519-7), C1 condition to C0 condition (uses epsilon delta on multiple of different differential value with limit needs smaller than continuous condition).  

N.B. if there's no C0 condition, with the condition x_next:=integrate^x_now(f(x) - some &alpha;), the prediction is valid because each of them are structure of the sum between first point and each of them. So if we can define integrate(original f) (and if it's continuous enough), the prediction is valid if original stream probability can be riemann integrable (and it's continuous enough).

N.B. DFT differential itself is ideal for trigometric function multiply sums. So applying tilt itself returns curved result, but this is reasonable one on the range with IDFT * DFT condition.

N.B. p0 uses northPole class to hack essential point whether or not if near the defined region.

# Tips on diff/integrate.
This p0 uses weak differential, so integrate(diff) != id because of complex part is omitted (but Re diff(diff) is no matter if we calculate in R or C).

# Tips on recursive on diff/integrate.
Frequency space prediction is equivalent to differential/integratial prediction in this. But it is equivalent to plain prediction in this.

# Tips on prediction vector recursive.
P0 calculates left differential by periodical one. To vanish right hand differential, we take sum of each range on next step vector. Whether this works well or not depends on left hand side distribution. But with some of differential sequence meaning, if we have some differential sequence average as 0 value, we can use them bared prediction.

# Tips start with complex function.
northPole class is formal laurent series essential point hack, if not both side is essential, multiplication inverse series works well, otherwise, the essential point hack enforces them into north pole near the defined region.
We need to northPole twice because only once sometimes doesn't converges bothside coefficients to 0 in general. But with the range they converges (non infinite values on any of samplepoints we have on the range), we only need to northPole once.
Because of this, we can handle any function which we can riemann integrable. Since we take probability, so if lebesgue integrate on original function is in C0, we can handle any function.

# Tips on some higher series.
There's a plenty of a space to extend this algorithm with &uparrow;, &downarrow; definition (they concludes exp(exp(...)) in C, inverse function is defined in randtools meaning inverse limit, but it's also in H if second operand is in C\\R.), they causes 3-variable function causes all of (x, f(x), status) handled on the series. With them, we can handle n-variable on reasonable function (because (n-1)-status block causes only the (x, f(x), status) triplet with extreme accuracy.).
But with this 1-variable predictor : (x, f(x)) any correspondence, if the parameter is large enough to input stream (as all of the status is on the taylor size), it's equivalent to the one.

-&gt; Knuth tower symbol seems to make some symmetric result on some of imaginary axises. This causes only the symmetric series on some of taylor series on H (with original series as j axis plotted). So this causes as the almost same prediction result via quarternion discrete fourier transform on syntax meaning, and they causes only the frequency-phase phase shifted results on the transformation. So we can define the differential on them on 2 of axises, but they seems to non information-rich one. And on information-rich meaning, we can apply this p0 to frequency series on complex space, they shall be equivalent to the one.

# Tips on o-minimal.
With referring wikipedia article: https://ja.wikipedia.org/wiki/%E5%AE%9F%E9%96%89%E4%BD%93 (2022/03/19), every weakly o-minimal structure is real closed field causes dim A = 1, 2, 4, 8 (https://ja.wikipedia.org/wiki/%E5%A4%9A%E5%85%83%E4%BD%93 (2022/03/19)), so if we can define {f(x) \| x in R} as semi-ordered one, octanion is enough dimension to calculate. If f has a zero divisor other than zero, the prediction fails in this.

# Tips on sedenion.
We calculate next step in sedenion, so with this, we can treat 2 dimension of semi-order as pseudo (x, status)-pair on input.

In advance and plenty of the space, there's a little possibility with this condition as to be fast and small space with operators on sedenion.

# Tips on appropriate handle higher dimension.
The real part and imginary part simple DFT is equivalent to plain prediction in P0. So any of summation based and probability based, sectional measurement based also in them. So to avoid this, P0DFT applies the lower prediction on abs and arg, this is independent from finite linear combination. -&gt; this is reverted because of p2/p.cc use.

# Tips on circular structure.
We take (summation ratio) - 1 multiple times. This causes handle the period as period on calculation because i-axis plotted sum f'/f goes near to log(f) and once goes log(f) + i &pi;/2, with applying some ratio twice or more goes arg(z) depend causes the period but not aligned ones on -1 origin point.

In this meaning, random permutation of sin(x) series is a difficult one, but with some small range permutation series applied by this method causes ok. Otherwise, the range we specify is small enough, so this fails.

# Tips on lebesgue measurable pseudo condition.
We gain sectional measurement expected prediction value by plain predictor, but if so, we should average some of input and skip them because of middle part of the frequency to be interpreted uncontinuous condition on original function lebesgue integrate.

# Tips on lebesgue measurable.
There exists non lebesgue measuarable set if axiom of choice with the set of greater or equal to aleph_1. This can be described in computer with covering a piece of paper completely with writing, even in algorithm meaning. But if so, the algorithm will be non-deterministic nor only the small finite series we can get. Otherwise, I have no idea the functions described in upto aleph_0 with constitutive method which is on non lebesgue measuarable set. So it would be greately appreciated if you know such set especially described in projection of infinite group combinations.

# Tips on one function recursion.
The tips below some lines are depends on one function recursion data series, but they're described with appending status block each points and only right multiply counter with x\+ := x \*\_(some group) status \*\_(some another group) counter and some f, ... f(x) ... f(x\+) ... data series.
To invert this method, we should take invariant on f with some appliable dimension.
Some of the preceding result on some fields shows non special differential equation with proper boundary condition into one dimension projection can be rewrited into some stable approximation block invariant with small error but distant points accumulate large error. This can be described with taylor matrix recursion causes unstable each elements, but with long enough before range with small next step condition, this is reduced. And any of the x^k is described with some range vector in this. (so arctan twice before to apply???)
So original function is just continuous or one-function recursion, they have a invariant structure on some short ranges.

# Tips on invariant size.
A longer invariants includes shorter invariant in trivial.
The invariant dimension size we need on original global function depends on the categorization of the structure of them.

# Tips on short invariant.
But with randtools meaning, we always have invariants with some of categorized ones in also small dimension. This causes linear summation combination on prediction of each invariants. Therefore, we sometimes have better results even in smaller invariant dimensions. But if the combination coefficient slowly moves in general, so no stable result be gained in smaller dimensions in general.
A recursive apply on this causes apply them into large number of clustered invariant categorization. So if there's a invariant combination largely changes point, we get gulf on the result. Otherwise, we get better result and usually, the global invariant unchange, so the condition satisfies this.

# Tips around th. condorcet's jury.
If input accuracy is rough enough in fixed point float and recursively apply p0 with smaller arguments, they causes \#\{structure\} \> \#\{combination\}, so they gain some slice in numerically if Th. Condorcet's jury is being applied, but the theorem can be used with blending PRNG into original stream. And with extending such input accuracy doesn't causes prediction structure change in this.

# Tips on jamming
There always exists jammer to any of the predictor. So to avoid bared predictor jammer, we apply all of the invariants smaller dimension than arguments. However, it is harder but there also exists such jammer on this predictor, they behaves like smaller invariants are balanced in linear combination as to be near 0 vector not to be predicted by some of the invariants. So with this condition, we are only able to do enlarge input argument, otherwise, all of the inputs are full rank differed to input argument, so it's information-rich in the status bits meaning on the stream. So because of them, we can avoid them if input argument is larger than input not to be full rank.

-\> this isn't used anymore. -\> retry to use this. -\> re-absent this.

# Tips on dimension
We vanish prediction dimension by 2 of pairs invert/plain, complex/sedenion prediction.

-\> for speed and memory, we use only plain and complex.

# Tips on 0 invariant chain.
If there's non categorized clustered invariant chain, there should be 0 invariant on the categorize chain. So to avoid this, we can pass the argument into p0 as linearly weighted, this causes slips the 0 chain into some of the large number of clustered ones. With this, if there's a invariants balanced whole sum 0 vector, we can break the balance of them. Also multiplication inverse condition breaks them (reverse order on frequency).

# Tips on all 0 invariant.
If we cannot get any of invariant into 0 invariant chain, it's only a return to the average series in expectation value. This is because \<\[x0, ..., xn, x\+\], a\>==0 satisfies only a==0 vector in expectation value, so after we get \[x0, ..., xn, x\+\], \<\[x0, ..., xn\], a\> == 0 for any a and \<\[x0, ..., xn, x\+\], b\> == 0 for any b in non constant meaning. So with constant, \<\[x0, ..., xn\], 0\> == 0 causes \<\[x0, ..., xn\], \[(n+1), ..., 1\]\> == const (walk average.) for such series in expectation value, differ 2 of them, they concludes \<\[x0, ..., xn, x\+\], 1\> == 0 in general.

# Tips on possible jamming
Even with there's measurable condition, the jam to this prediction can exist. To fight with them, we produce 3 of predicts described in p2/readme.md. We choose triplet because original function has dimension 3 or more because (const., x, f(x), status) has them, and const. part can be described by x\_k - x\_(k - 1), so in literally literally 3 dimension.

# Tips on 3 of choices
This p0 depends on the decomposition low-frequency part (plain p0), high-frequency part (multiplication inverse), linear part (minimum square) and no information be gained part (average). So this is a max rank to categorize input stream, also be a enough dimension to predict.
And linear part might be included some plain or inverse part, so integrate them.

# Tips on hard jamming
Even with some of the counter measure on jamming, if there exists slow and nonstatistical illegal jammer (almost every condition have on any of one predictor), they fails in long span. So we select simple ones which output triplets.
Also, jammer can jam out any two of the predictors into 0 expectation value condition because the jammer can cheat (f(x), status). The status is literally 1 dimension but since they are possible to have more dimensions virtualy, the predictor triplet is able to be being jammed in short range.
N.B. in the long range, if the input is deterministic stream, the predictor triplet cannot be jammed because status can be treated as real one dimension. But with modern PRNGs, the status has sporadic sub groups they appears in the f(x) stream rarely, so the gulf will appear on the prediction triplet. Also, if the input stream has a noised ones, since duumy data will be inserted, so prediction on long range will be hard ones because the hypothesis invariant can be noised. In the worst case, noise itself is something biased causes categorized invariant to be biased ones, the prediction itself disturbed.

# Tips on pseudo Brown walk
If we make pseudo Brown walk by frequency space \[(random.uniform(- 1, 1) + 1.j * random.uniform(-1, 1)) / abs(...), ...\] series, the prediction might be harder than original stream from this. (This is the analogy of shuffle operation via DFT).
So from numerical test, we might need the initial (might be supreme of) status dimension in arithmetic operation meaning to status length. So the below extends them arbitrary length from original PRNG.

# Tips on dimension 2
(x, status) pair can be (counter, status) pair, so if the counter is cyclic, it's in status in general, so we only need first 2 predictors if \|counter\| \< status.
But in general, \|counter\| -\> infinite, and if we cut some range and slide them, it is different value not fixed from counter causes real 1 dimension, we need to bet the return to average that is treated as 0 invariant in whole range.
However, if we predict with finite range, \|counter\| \< infinite, so they can integrated to status in literally, so if we only have finite and upper bounded size stream and predict with large enough status length, we can ignore the last prediction.

# Tips on non measurable condition
When the original stream doesn't have Lebesgue measurable condition in discrete, the description f(x):=<a,x>/<b,x> has the
hyperbolic triangular function complemented region, and the series they have countable number elements has the description for n elements,
sqrt(n) sample blocks division in whole.
The division itself is not only exists but also there might be multiple of them, so if the first point is the same condition, they're virtually measurable.
So this predictor adds such measurable condition softly.
-&gt; disabled in the latest version.

# Tips on insurance run.
We input lineary complemented stream on behalf of insurace run the result. Because the result of the stack is nonlinear, we can earn some better result on them.

# Tips on statistical illegal result.
We normalize input delta with logscale/expscale chain. With this, better stable result gained.

# Tips on recursive.
We do raw prediction twice with and without the condition normalization.
So the context itself will change on recursive apply result, twice is too much but we suppose the benefit first to shaggy, second to oval.
-&gt; disabled in the latest version, we always take to normalize. -&gt; disabled, this isn't needed with clang -O0 condition.

# Tips on improvement on continuity
We can logscale input then expscale output. This causes improve continuity but this isn't implemented.

# Tips on non linear
We can log scale input (both x,y-axis). This causes some series of automorphic form to get better result on this.

# Tips on adding continuity.
We can add some of the continuity by recursive logscale input expscale output.

Also, we can add some of them by sliding average window to make some of the riemann measurable condition.

p0l.cc doing some of this however, this needs huge memory usage and this isn't improve enough.

# Another Download Sites (Closed)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/

# Real close
2023/02/28
2023/03/09 bug fix after close #1.
2023/03/13 integrate all files into lieonn.hh after close #2.
2023/03/18 offset pnext, eliminate recursive pnext. after close #3.
2023/03/24 code clean, after close #4.
2023/03/29 pnext needs 3 &lt;= r, after close #5.
2023/03/31 merge prand.
2023/04/02 merge catg fix.
2023/04/03 merge.
2023/04/04 update readme.
2023/04/21 only single side northPole method, we don't need them on stability in such case.
2023/06/24 fix to avoid observation matters.
2023/07/01 argv range minimum value change.
2023/07/07 update .cc comment.
2023/10/30 add p0l to add some continuity. update readme.
2024/05/06 integrate p0l into p0, p0's argv\[1\]\<0 is integrated into bitsofcotton/p1.
2024/05/31 compile jammer.
2024/06/01 fix JAM.
2024/06/02 fix _JAM replace with _JAM_, this caused pseudo correct results.
2024/06/02 JAM into p2/cr.py.
2024/06/13 add p0i.cc.
2024/06/15 punch p0i.cc, it's no use.
2024/06/15 add progression.
2024/06/16 add progression \<0 arg.
2024/06/17 fix progression.

