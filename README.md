# p2
A predictor formatter for p0, p1 either p0, p1 integrator compete with patternizable jammer.

# Usage
    cr.py ([Rm]|[HD]?) ...
    
     *** reformation part ***
    # take delta     on input stream
    p d
    # take summation on input stream
    p s
    # take skip      on input stream
    p k <interval>
    # take skip head on input stream
    p S <margin>
    # take reverse   on input stream
    p v
    # pick {0,1} str on input stream
    p x
    # take reform [-1,1] on input stream
    p X
    # take reform [-1,1] on input stream without offset
    p Z
    # take inverse   on input stream
    p i
    # take picked column      on input stream
    p l <col0index> ...
    # take duplicate toeplitz on input stream
    p z <column number>
    # take multiply each      on input stream
    p t <ratio>
    # take offset   each      on input stream
    p o <offset>
    # take absolute each      on input stream
    p a
    # take sign     each      on input stream
    p b
    # take sum columns each line on input stream
    p G
    # take walk condition each on input stream
    p w <range>
    # take column 0 horizontal cut output to each column
    p E <number>
    # take opposite type output string each on input stream
    p F <bit number>
    
     *** PRNG part ***
    # make [-1,1]   PRNG stream
    p [rR]  <proto>
    # make {-1,0,1} PRNG stream
    p [rR]b <proto>
    # make {-1,1}   PRNG stream
    p [rR]B <proto>
    # blend [-1,1]  PRNG stream
    p m<proto> <number of output columns>
    # flip or not   PRNG stream
    p M<proto> <number of output columns>
    
     *** predictor part ***
    # predict with Riemann measureable condition (c for difference output)
    p 0c? <arg>
    # predict with untangle combination condition (c for difference output)
    p 1c? <arg>
    # feed patternizable jammer input entropy (C for difference output)
    p [cC] <arg>
    # jammer to the jammer output (+ for short fixed range target)
    p j+?
    # jam out input column 0 by input column 1+
    p Q
    # trivial id. prediction (plain for flip last, + for return to average)
    p I+
    # ddpmopt compatible prediction
    p Ac?
    
     *** vector operation part ***
    # input serial stream to vector stream
    p f <dimension>
    # input vector stream to serial stream
    p h
    # input vector stream to pgm graphics output
    p P
    # input vector stream to harden PRNG part
    p e
    
     *** multi process call part ***
    # do double prediction on same input
    p D <command set 0> <command set 1>
    # do each of all column input prediction parallel, take output column 0 each.
    p H <command>
    # do each of all column input prediction parallel, take output column 0 as a prediction value, prediction average * input stream average output.
    p @ <command>
    
     *** other part ***
    # multiple file load into same line columns
    p L <file0> ...
    # show output statistics whether arg<|x - 1/2|*2<{1-arg,infty} (+ for infty)
    p T+ <arg>

    # N.B. there's many much of the combibnations, if we upload them on
    #      this README.md they'd slip from somehow in our experience some years
    #      ago. so we absent intentionally blank them.

# Description:
Predictor formatter and some of the toolset fighting with jammers for non usual input streams however, even this predictor can have the jammer to us.

# Another Download Sites (Leave)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/

# Real close
2023/02/28
2023/03/09 bug fix after close #1.
2023/03/13 valid only short range bug fix after close #2.
2023/03/13 integrate all files into lieonn.hh after close #3.
2023/03/18 merge latest p0, after close #4.
2023/03/22 update README from osdn.net, after close #5.
2023/03/24 code clean, after close #6.
2023/03/31 persistent resistance to integer calculation jamming, after close #7.
2023/04/01 source code comment remedy from osdn.jp generic-p2 summary, after close #8.
2023/04/02 balance between jammer and predictor.
2023/04/03 merge.
2023/04/05 update readme.
2023/04/14 additional closed .cc files.
2023/04/21 make/revert ProgramInvariant algorithm change.
2023/05/18 single multiple single fix, so input y-axis ratio is the matter. also fix rand.cc for MAX_RAND range for knuth_b.
2023/05/18 cr.py rewrote to per single command.
2023/06/07 delete pfork.py, add pchain.sh, update readme.
2023/06/18 update readme, p.cc comment (no logic change), cr.py p command.
2023/06/24 fix to avoid observation matters.
2023/07/01 op2.py m command relative change, fix indent on p.cc.
2023/07/07 update readme.md, update comment in .cc files.
2023/07/08 delete rand.cc, cr.py r is enough to use. invariant causes +1. update readme, eliminate rand.cc, p.cc, comments are moved into readme. op2.py integrated into cr.py.
2023/08/02 update cr.py, readme, lieonn.hh P012L::next.
2023/09/20 update cr.py e command.
2023/09/25 some prng improvement test and update readme.
2023/10/26 update cr.py v command use with m command.
2023/10/30 copy structure reliably with randtools meaning.
2024/04/25 add tips on simple controller.
2024/04/29 add pr4.cc .
2024/04/30 add tips on bothside controller condition. brush up readme.md.
2024/05/05 add pr4c.cc update readme, pr4.cc .
2024/05/06 add tips on jammer, fix some of the readme. eliminate unused commands in cr.py. fix pr4.cc, integrate pr4c.cc into pr4.cc with definition. add tips when being targetted. BRUSH UP README.MD, so some of the description had be DELETED.
2024/05/07 fix XXX, add tips on PRNGs. really close.
2024/05/31 compile jammer.
2024/06/01 fix JAM, also our system is infected because of the test.
2024/06/02 fix readme, it's from p0.cc:_JAM mis spell.
2024/06/02 add cr.py j command, elim _JAM.
2024/06/14 fix cr.py l command, add cr.py 0, - command.
2024/06/15 add progression and readme.
2024/06/16 add progression \<0 argv.
2024/06/17 add readme, something goes well with some of the PRNGs with static, absent repository with this. also fix progression.
2024/06/19 merge latest lieonn.
2024/06/21 merge latest lieonn.
2024/06/21 add p2prng.cc, update readme, merge latest lieonn. INCLUDES command line argument change.
2024/06/22 update readme, p01 fatal fix. make/revert program invariants change friendly to predictors.
2024/06/23 large change around class instance initializer, also have progression short range fix.
2024/06/23 fatal fix around last update.
2024/06/24 fix addp == true progression case.
2024/06/26 cr.py o command and so on.
2024/06/26 fix Ppersistent.
2024/07/07 Pprogression uses shorter range but enough internal states.
2024/09/09 add z cmd, change s cmd.
2024/09/22 add pp3dft.cc, elim p210.cc . update readme.
2024/09/22 integrate pp3dft.cc into p1/pp3.cc so elim them, releave.
2024/09/27 add usual prediction section readme.
2024/09/28 eliminate Pprogression, it's a waste of the accuracy.
2024/09/29 update readme. fix P012L::next avg.O().
2024/10/26 update readme.
2024/11/01 update readme, add N command to cr.py.
2024/11/14 add cr.py p command to continue whole input p0 with fft without accuracy.
2024/11/16 fix readme.
2024/12/05 fix readme. exchanged argv[1] and argv[2] meaning.
2024/12/08 update readme compat with latest p0, p1.
2024/12/13 should really leave, close readme.md.
2024/12/14 update readme.md, might still have some glitches on readme.md after update.
2024/01/27 update readme.md.
2025/02/05 add f command to cr.py.
2025/02/07 add persistent.cc .
2025/02/10 fix and refactor persistent.cc to pp8q.cc (pp8 - .25).
2025/02/11 add readme.md persistent.cc concerns.
2025/02/13 a little fix around Persistent.
2025/02/27 append readme.md, eliminate variable step length.
2025/03/01 update readme eliminate last section, append second last section.
2025/03/04 add _CHAIN_ compile option to catgp.cc, merge latest lieonn.
2025/03/08 add B command for cr.py they returns predicted result and original input pair slided causes some of the controlled condition original stream dump.
2025/03/09 merge latest lieonn.
2025/03/17 revert step param.
2025/03/19 add g, G commands to cr.py.
2025/03/22 close with this README.md.
2025/03/28 after burn README.md.
2025/03/30 adding test result glance around xeb scored ones.
2025/04/01 append p, q command to cr.py, add readme.md, real close.
2025/04/02 append 2 command and readme.md notation, real close.
2025/04/03 compat with simplefloat p, 2 commands also fixed fatal error on them.
2025/04/04 add y command for condorcet jury with original stream intensity. also simplify separate readme.md on same file.
2025/04/04 fix y command, add l command option, readme update, real close.
2025/04/05 fix y command median condition, update readme, really close.
2025/04/06 fix readme.md usage. ok really close.
2025/04/16 fix p2prng.cc GETENTROPY pointer and value mis interpretation, I don't know why compiler doesn't mattered this on our environment even the running result was so.
2025/04/17 merge latest dimension auto tuner from ddpmopt. fixed step argv works well.
2025/04/18 eliminate step parameter, they doesn't improve results.
2025/04/19 merge latest lieonn.
2025/05/02 add p2.cc predictor for catgp p1 p0 chain with different x-asis scaling twice with different step.
2025/05/03 Select the strategy to separate binary x-axis linear and nonlinear case. update readme close condition. code clean.
2025/05/04 fix last up arctanFeeder initial scaling.
2025/05/05 pseudo-harden prng command. close.
2025/05/10 add cr.py H cmd, change cr.py B cmd. update readme.
2025/05/11 change cr.py E cmd, elim M cmd, add actual sectional measurement. update readme.
2025/05/14 fix last of usage, add readme.md close section.
2025/05/15 change cr.py M cmd meaning, add usage.
2025/05/16 change cr.py M cmd meaning, add X, Y command, freeze close because we cannot find our targetted data stream isn't affected by the universal one. What we need is graphics predictor, we move to ddpmopt impl.
2025/05/17 fix typo on readme.md, we don't know why this mixed into readme.md.
2025/05/17 elim Y dependant commands.
2025/05/19 applying p0p can have slim cr.py result.
2025/05/27 update cr.py locking, D cmd output +1 original col, readme.md.
2025/05/31 cr.py Q command match the double of input stream. update readme.
2025/06/03 fix crash on p2... around zero division on _FLOAT_BITS_ specified case, not with cpu float they can crash in some condition. cr.py Q cmd comment change, G command chg.
2025/06/04 import Q command chain try into lieonn.
2025/06/05 fix p2qt as reasonable one but not with real pred-value. add B cmd to cr.py.
2025/06/06 ok p2qt for our machine which original stream seems better continuous on ours.
2025/06/06 exclude p0maxRank(idFeeder(3)) prediction from p2qt into p0.
2025/06/06 add p2qt argv causes ok result for us test stream on our machine. update readme.
2025/06/07 update readme. change reasonable implementation p2qt.cc and around them. elim p2.cc and concerns, it's not used deeply in any of the predictors.
2025/06/07 revert, Q command result sign doesn't depend on the place sign, so fixed.
2025/06/08 fix p2qt.cc param != 0 case but our machine is infected from the conclusion, either needs {-1, 1}-input stream.
2025/06/08 we're in infected condition, so we rework possible thin layered condition.
2025/06/08 we move p2qt.cc into p2sg.cc correct meaning from test result also update readme. also the concerned functions refactoring. fix catgp 0 crash.
2025/06/09 move into p2cr from p2prng.cc and cr.py concerned parts this makes python to exclude from prediction chain however this once improves enough but will slips. also fix p2cr Q cmd align. fix p2sg as invertible.
2025/06/09 merge catgp.cc, p2sg,cc, p2cr.cc into p2.cc one binaried. (should we integrate p0, p1, p2 into p as one binaried?). refactor readme.md.
2025/06/10 persistent uint32_t use in calculation option but they also use int32_t, also void*.
2025/06/11 compat compile option to gcc4.2.1.
2025/06/12 compat compile option with one variant of gcc2.95.3.
2025/06/17 merge latest ddpmopt fix.
2025/06/19 merge latest ddpmopt fix. add M command on p2.cc. add cr.py H[+-] command.
2025/06/20 cleanup some commands and source code. fix around persistent especially operator >> accuracy on simplefloat. Intentionally absent documentation but on the source codes.
2025/06/20 code cleaning, cr.py to p2.cc implement.
2025/06/22 change T command output to better reasonable one, PERSISTENT improve.
2025/06/23 code clean, flush.
2025/06/25 merge latest ddpmopt result. also need to implant readme.md into lieonn.hh.
2025/06/25 readme.md moved into implant lieonn.hh as a comment. add usage into p2.cc.
2025/06/25 fix usage also 'e' command crash, update lieonn.hh comment.
2025/06/28 refactor and fix around lieonn, re-compat with gcc2953.
2025/06/29-30 merge latest ddpmopt change, also add Ac? command also T command IMPORTANT change.
2025/07/01 merge latest ddpmopt change causes various bug fixes.
2025/07/02-03 merge latest ddpmopt causes new prediction algorithm debug ok.
2025/07/04 merge latest ddpmopt debug ok.

