# usage: awk -v excl=0 -v tol=1e-4 -f allp.awk out.txt
# ddpmopt k+ only (rankdata removed). 1 block = 56 lines:
#   v[0..2]  orig target     : med avg rand63
#   v[3]     orig target     : last
#   v[4..6]  replaced target : med avg rand63
#   v[7..54] 16 sign patterns t=0..15, each: med avg rand63
#   last line "x, 4, 5"  (x = sign mean of rand63; |x|=1 -> constant, excluded if excl=1)

function abs(x) { return x < 0 ? -x : x }

# mid-rank position of val among the 16 |c_t| (ties count half)
function midp(pred, val,    u, w, g, e) {
  g = 0; e = 0
  for (u = 0; u < 16; u++) {
    w = abs(ct[pred, u])
    if (w < val - tol) g++
    else if (w <= val + tol) e++
  }
  return (g + 0.5 * e) / 16
}

function block(    pred, p, t, cc, obs, rep, a, s, q, ea, es, va, vs, sm, sm2, mpt, vmp, mpo, mpr, ok, dd, d, ed, qd, vd, od, L) {
  nblk++
  for (pred = 0; pred < 3; pred++) {
    obs = v[pred]
    rep = v[4 + pred]
    a = 0; s = 0; q = 0
    for (t = 0; t < 16; t++) {
      cc = v[7 + 3 * t + pred]
      ct[pred, t] = cc
      a += abs(cc); s += cc; q += cc * cc
    }
    ea = a / 16; es = s / 16; q = q / 16
    va = q - ea * ea; vs = q - es * es
    sm = 0; sm2 = 0
    for (t = 0; t < 16; t++) {
      mpt = midp(pred, abs(ct[pred, t]))
      sm += mpt; sm2 += mpt * mpt
    }
    vmp = sm2 / 16 - (sm / 16) * (sm / 16)
    mpo = midp(pred, abs(obs))
    mpr = midp(pred, abs(rep))

    n[pred]++
    so[pred] += abs(obs)
    sE[pred] += ea
    da[pred] += abs(obs) - ea
    Va[pred] += va
    ds[pred] += obs - es
    Vs[pred] += vs
    dm[pred] += mpo - 0.5
    Vm[pred] += vmp
    dr[pred] += abs(rep) - ea
    dmr[pred] += mpr - 0.5
    dd = abs(obs) - abs(rep)
    d1[pred] += dd
    d2[pred] += dd * dd

    ok = 0
    for (t = 0; t < 16; t++) if (abs(ct[pred, t] - obs) <= tol) ok = 1
    mat[pred] += ok

    for (t = 0; t < 16; t++) {
      cc = ct[pred, t]
      if (n[pred] == 1) { mn[pred, t] = cc; mx[pred, t] = cc }
      else {
        if (cc < mn[pred, t]) mn[pred, t] = cc
        if (cc > mx[pred, t]) mx[pred, t] = cc
      }
    }
  }
  # paired differences: med - avg, med - rand63 (same target => enumerate jointly)
  for (p = 1; p < 3; p++) {
    ed = 0; qd = 0
    for (t = 0; t < 16; t++) {
      d = abs(ct[0, t]) - abs(ct[p, t])
      ed += d; qd += d * d
    }
    ed = ed / 16; vd = qd / 16 - ed * ed
    od = abs(v[0]) - abs(v[p])
    np[p]++
    pd[p] += od - ed
    pV[p] += vd
    pm[p] += od
  }
  L = v[3]
  ls += L; las += abs(L)
  if (abs(L) > 0.999) l1++
}

BEGIN {
  if (tol == "") tol = 1e-4
  if (excl == "") excl = 0
  pn[0] = "med"; pn[1] = "avg"; pn[2] = "rand63"
}

/^(med|avg|rand_pgm|last_rand)/ { v[c++] = $NF + 0; next }

/, 4, 5$/ {
  x = $1 + 0
  if (c != 55) { nbad++; printf "warn: skip block (%d values)\n", c > "/dev/stderr" }
  else if (excl && abs(x) > 0.999) nex++
  else block()
  c = 0
}

END {
  printf "blocks used=%d  excluded(rand63 const)=%d  bad=%d  tol=%g\n\n", nblk, nex, nbad, tol
  if (nblk == 0) exit
  print "[1] null-corrected |cos| per predictor, offset (k+) only (|z|<~2 = indistinguishable from null)"
  printf "%-6s %4s %8s %8s %8s %7s %7s %7s %7s %6s %7s %7s\n", "pred", "N", "mean|o|", "E|c|", "diff", "se", "z_abs", "z_sgn", "z_midp", "match", "cal_abs", "cal_mp"
  for (pred = 0; pred < 3; pred++) {
    N = n[pred]
    se = sqrt(Va[pred]) / N
    printf "%-6s %4d %8.4f %8.4f %+8.4f %7.4f %7s %7s %7s %5.0f%% %7s %7s\n", pn[pred], N, \
      so[pred] / N, sE[pred] / N, da[pred] / N, se, \
      (Va[pred] > 0 ? sprintf("%+.2f", da[pred] / sqrt(Va[pred])) : "NA"), \
      (Vs[pred] > 0 ? sprintf("%+.2f", ds[pred] / sqrt(Vs[pred])) : "NA"), \
      (Vm[pred] > 0 ? sprintf("%+.2f", dm[pred] / sqrt(Vm[pred])) : "NA"), \
      100 * mat[pred] / N, \
      (Va[pred] > 0 ? sprintf("%+.2f", dr[pred] / sqrt(Va[pred])) : "NA"), \
      (Vm[pred] > 0 ? sprintf("%+.2f", dmr[pred] / sqrt(Vm[pred])) : "NA")
  }
  print ""
  print "[2] med vs other, null-corrected (paired over the same 16 targets)"
  printf "%-12s %4s %9s %9s %7s %7s\n", "pair", "N", "raw diff", "adj diff", "se", "z"
  for (p = 1; p < 3; p++) {
    N = np[p]
    printf "%-12s %4d %+9.4f %+9.4f %7.4f %7s\n", "med-" pn[p], N, pm[p] / N, pd[p] / N, \
      sqrt(pV[p]) / N, (pV[p] > 0 ? sprintf("%+.2f", pd[p] / sqrt(pV[p])) : "NA")
  }
  print ""
  print "[3] orig target vs replaced target: mean(|obs|-|repl|) +- se (empirical, independent of the null model)"
  for (pred = 0; pred < 3; pred++) {
    N = n[pred]; m = d1[pred] / N
    v2 = (N > 1) ? (d2[pred] / N - m * m) / (N - 1) : 0
    printf "%-6s %+8.4f +- %.4f\n", pn[pred], m, (v2 > 0 ? sqrt(v2) : 0)
  }
  print ""
  print "[4] diagnostics"
  for (pred = 0; pred < 3; pred++) {
    sp = 0
    for (t = 0; t < 16; t++) { s = mx[pred, t] - mn[pred, t]; if (s > sp) sp = s }
    printf "%-6s max spread of c_t across blocks = %.6f%s\n", pn[pred], sp, \
      (sp < tol ? "   <-- CONSTANT across blocks (carries no per-block information)" : "")
  }
  printf "last: mean=%+.3f  mean|.|=%.3f  frac(|last|=1)=%.3f\n", ls / nblk, las / nblk, l1 / nblk
  print "match = share of blocks whose observed cos equals one of the 16 pattern cos (should be ~100%)"
  print "cal_*  = same z on the replaced (independent) target: should look like N(0,1)"
}
