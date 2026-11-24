# usage: awk -v excl=0 -v tol=1e-4 -f all.awk out.txt   (default excl=0: keep blocks where rand63 is constant)
#
# 1 block = 111 lines:
#   v[0..3]   rank   orig target : med avg rand63 last
#   v[4..7]   offset orig target : med avg rand63 last
#   v[8..10]  rank   replaced target : med avg rand63
#   v[11..13] offset replaced target : med avg rand63
#   v[14..109] 16 sign patterns t=0..15, each: rank(med avg rand63) offset(med avg rand63)
#   last line "x, 4, 5" (x = sign mean of rand63; |x|=1 -> constant, excluded if excl=1)

function abs(x) { return x < 0 ? -x : x }

# mid-rank position of val among the 16 |c_t| (share below val, ties count half); larger = |obs| larger than null
function midp(col, pred, val,    u, w, g, e) {
  g = 0; e = 0
  for (u = 0; u < 16; u++) {
    w = abs(ct[col, pred, u])
    if (w < val - tol) g++
    else if (w <= val + tol) e++
  }
  return (g + 0.5 * e) / 16
}

function block(    col, pred, p, t, cc, obs, rep, a, s, q, ea, es, va, vs, sm, sm2, mpt, vmp, mpo, mpr, ok, dd, d, ed, qd, vd, od, L) {
  nblk++
  for (col = 0; col < 2; col++) {
    for (pred = 0; pred < 3; pred++) {
      obs = v[4 * col + pred]
      rep = v[8 + 3 * col + pred]
      a = 0; s = 0; q = 0
      for (t = 0; t < 16; t++) {
        cc = v[14 + 6 * t + 3 * col + pred]
        ct[col, pred, t] = cc
        a += abs(cc); s += cc; q += cc * cc
      }
      ea = a / 16; es = s / 16; q = q / 16
      va = q - ea * ea; vs = q - es * es
      sm = 0; sm2 = 0
      for (t = 0; t < 16; t++) {
        mpt = midp(col, pred, abs(ct[col, pred, t]))
        sm += mpt; sm2 += mpt * mpt
      }
      vmp = sm2 / 16 - (sm / 16) * (sm / 16)
      mpo = midp(col, pred, abs(obs))
      mpr = midp(col, pred, abs(rep))

      n[col, pred]++
      so[col, pred] += abs(obs)
      sE[col, pred] += ea
      da[col, pred] += abs(obs) - ea
      Va[col, pred] += va
      ds[col, pred] += obs - es
      Vs[col, pred] += vs
      dm[col, pred] += mpo - 0.5
      Vm[col, pred] += vmp
      dr[col, pred] += abs(rep) - ea
      dmr[col, pred] += mpr - 0.5
      dd = abs(obs) - abs(rep)
      d1[col, pred] += dd
      d2[col, pred] += dd * dd

      ok = 0
      for (t = 0; t < 16; t++) if (abs(ct[col, pred, t] - obs) <= tol) ok = 1
      mat[col, pred] += ok

      for (t = 0; t < 16; t++) {
        cc = ct[col, pred, t]
        if (n[col, pred] == 1) { mn[col, pred, t] = cc; mx[col, pred, t] = cc }
        else {
          if (cc < mn[col, pred, t]) mn[col, pred, t] = cc
          if (cc > mx[col, pred, t]) mx[col, pred, t] = cc
        }
      }
    }
    # paired differences: med - avg, med - rand63 (same target => dependent, so enumerate jointly)
    for (p = 1; p < 3; p++) {
      ed = 0; qd = 0
      for (t = 0; t < 16; t++) {
        d = abs(ct[col, 0, t]) - abs(ct[col, p, t])
        ed += d; qd += d * d
      }
      ed = ed / 16; vd = qd / 16 - ed * ed
      od = abs(v[4 * col]) - abs(v[4 * col + p])
      np[col, p]++
      pd[col, p] += od - ed
      pV[col, p] += vd
      pm[col, p] += od
    }
    L = v[3 + 4 * col]
    ls[col] += L; las[col] += abs(L)
    if (abs(L) > 0.999) l1[col]++
  }
}

BEGIN {
  if (tol == "") tol = 1e-4
  if (excl == "") excl = 0
  cn[0] = "rank"; cn[1] = "offset"
  pn[0] = "med"; pn[1] = "avg"; pn[2] = "rand63"
}

/^(med|avg|rand_pgm|last_rand)/ { v[c++] = $NF + 0; next }

/, 4, 5$/ {
  x = $1 + 0
  if (c != 110) { nbad++; printf "warn: skip block (%d values)\n", c > "/dev/stderr" }
  else if (excl && abs(x) > 0.999) nex++
  else block()
  c = 0
}

END {
  printf "blocks used=%d  excluded(rand63 const)=%d  bad=%d  tol=%g\n\n", nblk, nex, nbad, tol
  if (nblk == 0) exit
  print "[1] null-corrected |cos| per predictor (z_* use the exact 16-pattern null; |z|<~2 = indistinguishable from null)"
  printf "%-6s %-6s %4s %8s %8s %8s %7s %7s %7s %7s %6s %7s %7s\n", "col", "pred", "N", "mean|o|", "E|c|", "diff", "se", "z_abs", "z_sgn", "z_midp", "match", "cal_abs", "cal_mp"
  for (col = 0; col < 2; col++) for (pred = 0; pred < 3; pred++) {
    N = n[col, pred]
    se = sqrt(Va[col, pred]) / N
    printf "%-6s %-6s %4d %8.4f %8.4f %+8.4f %7.4f %7s %7s %7s %5.0f%% %7s %7s\n", cn[col], pn[pred], N, \
      so[col, pred] / N, sE[col, pred] / N, da[col, pred] / N, se, \
      (Va[col, pred] > 0 ? sprintf("%+.2f", da[col, pred] / sqrt(Va[col, pred])) : "NA"), \
      (Vs[col, pred] > 0 ? sprintf("%+.2f", ds[col, pred] / sqrt(Vs[col, pred])) : "NA"), \
      (Vm[col, pred] > 0 ? sprintf("%+.2f", dm[col, pred] / sqrt(Vm[col, pred])) : "NA"), \
      100 * mat[col, pred] / N, \
      (Va[col, pred] > 0 ? sprintf("%+.2f", dr[col, pred] / sqrt(Va[col, pred])) : "NA"), \
      (Vm[col, pred] > 0 ? sprintf("%+.2f", dmr[col, pred] / sqrt(Vm[col, pred])) : "NA")
  }
  print ""
  print "[2] med vs other, null-corrected (paired over the same 16 targets)"
  printf "%-6s %-12s %4s %9s %9s %7s %7s\n", "col", "pair", "N", "raw diff", "adj diff", "se", "z"
  for (col = 0; col < 2; col++) for (p = 1; p < 3; p++) {
    N = np[col, p]
    printf "%-6s %-12s %4d %+9.4f %+9.4f %7.4f %7s\n", cn[col], "med-" pn[p], N, pm[col, p] / N, pd[col, p] / N, \
      sqrt(pV[col, p]) / N, (pV[col, p] > 0 ? sprintf("%+.2f", pd[col, p] / sqrt(pV[col, p])) : "NA")
  }
  print ""
  print "[3] orig target vs replaced target: mean(|obs|-|repl|) +- se (empirical, independent of the null model)"
  for (col = 0; col < 2; col++) for (pred = 0; pred < 3; pred++) {
    N = n[col, pred]; m = d1[col, pred] / N
    v2 = (N > 1) ? (d2[col, pred] / N - m * m) / (N - 1) : 0
    printf "%-6s %-6s %+8.4f +- %.4f\n", cn[col], pn[pred], m, (v2 > 0 ? sqrt(v2) : 0)
  }
  print ""
  print "[4] diagnostics"
  for (col = 0; col < 2; col++) for (pred = 0; pred < 3; pred++) {
    sp = 0
    for (t = 0; t < 16; t++) { s = mx[col, pred, t] - mn[col, pred, t]; if (s > sp) sp = s }
    printf "%-6s %-6s max spread of c_t across blocks = %.6f%s\n", cn[col], pn[pred], sp, \
      (sp < tol ? "   <-- CONSTANT across blocks (carries no per-block information)" : "")
  }
  for (col = 0; col < 2; col++)
    printf "%-6s last: mean=%+.3f  mean|.|=%.3f  frac(|last|=1)=%.3f\n", cn[col], ls[col] / nblk, las[col] / nblk, l1[col] / nblk
  print "match = share of blocks whose observed cos equals one of the 16 pattern cos (should be ~100%)"
  print "cal_*  = same z on the replaced (independent) target: should look like N(0,1)"
}

# Claude.ai (2026/11/24) generated.
