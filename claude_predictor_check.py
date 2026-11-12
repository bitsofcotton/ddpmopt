#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
predictor_check.py  --  予測器 p が元画像 x をどの程度予測できているかを検定する(1 ファイル版)

ここまでの議論で固定した指標
------------------------------------------------------------------------------
主指標(PRIMARY。これだけを「結果」として報告する)
  周囲 margin 画素(既定 4)を除いた領域で、avg と last を差し引いた p と x の
  偏順位相関(partial Spearman)。ヌルは x の位相ランダム化 surrogate。
  向きは「正」で事前固定。片側 p・両側 p・z・ブロックブートストラップ区間を出す。

補助指標(SECONDARY。主指標の読み方を助ける)
  raw     : 差し引きなしの素の Spearman
  sign    : 符号一致率(中央値で 2 分割、偶然水準 0.5)と kappa・バランス精度。ヌルは AAFT
  diag    : 偏相関をブロック寄与 u に分解した u_pos / u_neg / mean_rb2
            (u_neg の z が「正」なら負の寄与が偶然より抑えられている = 正の偏り)
  edge    : 周囲の除外幅 0 / 4 / 8 での偏相関(端の列への依存を見る)
  coarse  : S×S ブロック平均した画像どうしの偏相関(「全体の質感」のスケール。既定 S=8)
  border  : 偏相関への寄与の、外周(幅 4)と内部への分解。x・p・avg・last の外周の平均順位も出す

任意(--extras)  : コントラスト層別(hi/lo)と周波数帯分解。記述的な補足で、判定には使わない。
任意(--calibrate): x を大きくずらした疑似陰性で、検定の z が N(0,1) に近いかを確認。
任意(--selftest) : 合成データでの陽性・陰性対照(動作確認)。

重要な前提
  * avg は「x を含まない」平均でなければならない(以前の混入事故の再発防止)。
  * p・x・avg・last は同じ形の 2 次元配列(.npy / PNM(P2,P3,P5,P6))。
  * 画像 1 枚では結論は出ない。--manifest で複数画像を回して集計する。
    同じ x を共有する複数の予測は独立な証拠ではない。--ensemble-by-x で x ごとに束ね(p どうしの符号で
    向きを揃えて平均。x は使わない)、異なる x を単位に集計する。向きは任意になるので、束ねた結果は両側で読む。
    集計では、平均 rho(主)のほかに、符号を見ない補助(mean|rho|, sum z^2)と、
    sha256 による重複(同じ p,x の組・同じ x の共有)の検出も出す。補助は探索的で判定には使わない。
    同じ x に対する別々の予測は独立な証拠ではない(対象が別の画像を使うこと)。

フレーム番号による時間順の検査(任意。取り違え対策)
  manifest に x_idx, last_idx, avg_max_idx 列を足す(または --x-idx --last-idx --avg-max-idx)と、
  last / avg が x より前のフレームだけから作られているかを検査し、違えばエラーで止める。

使い方
  python3 predictor_check.py --p p.pgm --x x.pgm --avg avg.pgm --last last.pgm --out result
  python3 predictor_check.py --manifest manifest.csv --out result      # 列: name,p,x,avg,last
  python3 predictor_check.py --selftest --out selftest
出力は result.json と result.txt(リダイレクトに頼らず、各画像の終了ごとに書き出す)。
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time
import traceback

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import rankdata

# ----------------------------------------------------------------------------
# バージョン(確認用の実行の前に固定し、結果にスクリプトの sha256 を残す)
# ----------------------------------------------------------------------------
VERSION = "predictor_check v1.0 (2026-11-12)"

# ----------------------------------------------------------------------------
# 固定パラメータ(変更したら報告に明記すること)
# ----------------------------------------------------------------------------
DEFAULTS = dict(margin=4, block=4, n_null=2000, seed=0, n_boot=2000,
                margins=(0, 4, 8), p_thr="median", x_thr="median")


# ----------------------------------------------------------------------------
# 入出力
# ----------------------------------------------------------------------------
def read_pnm(path):
    with open(path, "rb") as f:
        data = f.read()
    pos = 0

    def token():
        nonlocal pos
        while True:
            while pos < len(data) and data[pos:pos + 1].isspace():
                pos += 1
            if pos < len(data) and data[pos:pos + 1] == b"#":
                while pos < len(data) and data[pos:pos + 1] not in (b"\n", b"\r"):
                    pos += 1
            else:
                break
        start = pos
        while pos < len(data) and not data[pos:pos + 1].isspace():
            pos += 1
        return data[start:pos]

    magic = token().decode()
    if magic not in ("P2", "P3", "P5", "P6"):
        raise ValueError("未対応の形式: %r" % magic)
    w, h, maxval = int(token()), int(token()), int(token())
    ch = 3 if magic in ("P3", "P6") else 1
    if magic in ("P5", "P6"):
        pos += 1
        dt = ">u2" if maxval > 255 else "u1"
        arr = np.frombuffer(data, dtype=dt, count=w * h * ch, offset=pos)
    else:
        arr = np.array(data[pos:].split(), dtype=float)[: w * h * ch]
    arr = arr.astype(float).reshape(h, w, ch)
    if ch == 3:
        return arr @ np.array([0.299, 0.587, 0.114])
    return arr[:, :, 0]


def load_array(path):
    if path.endswith(".npy"):
        a = np.load(path)
    else:
        try:
            a = read_pnm(path)
        except Exception:
            from PIL import Image  # 任意(入っていれば使う)
            a = np.asarray(Image.open(path).convert("L"), dtype=float)
    a = np.asarray(a, dtype=float)
    if a.ndim != 2:
        raise ValueError("%s: 2 次元配列ではありません shape=%s" % (path, a.shape))
    return a


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.round(6).tolist()
    if isinstance(o, np.generic):
        return o.item()
    return o


# ----------------------------------------------------------------------------
# 基本部品
# ----------------------------------------------------------------------------
def rank01(a):
    a = np.asarray(a, float).ravel()
    return (rankdata(a, method="average") - 0.5) / a.size


def spearman(a, b):
    return float(np.corrcoef(rank01(a), rank01(b))[0, 1])


def phase_surrogate(x, rng):
    """位相ランダム化(振幅スペクトル = 自己相関を保つ)。順位化して使う前提。"""
    F = np.fft.fft2(x - x.mean())
    ph = np.angle(np.fft.fft2(rng.standard_normal(x.shape)))
    return np.real(np.fft.ifft2(np.abs(F) * np.exp(1j * ph)))


def aaft_surrogate(x, rng):
    """周辺分布(値の集合・平均・符号の割合)を厳密に保ち、スペクトルをほぼ保つ。"""
    flat = x.ravel()
    order = np.argsort(np.argsort(flat))
    g = np.sort(rng.standard_normal(flat.size))[order]
    G = np.fft.fft2(g.reshape(x.shape))
    ph = np.angle(np.fft.fft2(rng.standard_normal(x.shape)))
    g2 = np.real(np.fft.ifft2(np.abs(G) * np.exp(1j * ph))).ravel()
    out = np.empty_like(flat)
    out[np.argsort(g2)] = np.sort(flat)
    return out.reshape(x.shape)


def prep(a, margin, B):
    """周囲 margin 画素を除き、さらにブロック幅 B の倍数に切り詰める。"""
    if margin:
        a = a[margin:-margin, margin:-margin]
    H, W = a.shape
    return a[: H // B * B, : W // B * B]


def make_resid(controls):
    """定数 + 各 control の順位への射影を引く(SVD で共線性にも頑健)。"""
    n = controls[0].size if controls else None
    cols = [np.ones(n)] if controls else None
    if not controls:
        return lambda y: y - y.mean()
    cols += [rank01(c) for c in controls]
    U, s, _ = np.linalg.svd(np.column_stack(cols), full_matrices=False)
    Q = U[:, s > 1e-10 * s[0]]
    return lambda y: y - Q @ (Q.T @ y)


def pvals(obs, null, two_sided=False):
    m = null.mean()
    up = (1 + (null >= obs).sum()) / (1 + len(null))
    if two_sided:
        return up, (1 + (np.abs(null - m) >= abs(obs - m)).sum()) / (1 + len(null))
    return up


# ----------------------------------------------------------------------------
# 中核: 偏順位相関 + ブロック寄与の診断(すべて同じ surrogate から出す)
# ----------------------------------------------------------------------------
def core_test(p, x, controls, B, n_null, rng, n_boot=2000):
    """p, x, controls: 周囲除外・ブロック整列済みの同形 2 次元配列。"""
    H, W = x.shape
    hb, wb = H // B, W // B
    resid = make_resid(controls)
    rP = resid(rank01(p))
    nP = np.linalg.norm(rP)

    def bsum(flat):
        return flat.reshape(H, W).reshape(hb, B, wb, B).sum((1, 3))

    syy = bsum(rP * rP)

    def stats(xv, full=False):
        rX = resid(rank01(xv))
        nX = np.linalg.norm(rX)
        c = bsum(rX * rP)
        sxx = bsum(rX * rX)
        u = c / (nX * nP)
        ok = (sxx * syy) > 0
        rb = np.zeros_like(c)
        rb[ok] = c[ok] / np.sqrt(sxx[ok] * syy[ok])
        v = np.array([u.sum(), u[u > 0].sum(), u[u < 0].sum(),
                      (rb[ok] ** 2).mean() if ok.any() else 0.0])
        return (v, c, sxx, u) if full else v

    obs, c, sxx, u = stats(x.ravel(), full=True)
    null = np.array([stats(phase_surrogate(x, rng).ravel()) for _ in range(n_null)])
    names = ["rho", "u_pos", "u_neg", "mean_rb2"]
    out = {}
    for k, nm in enumerate(names):
        m, s = null[:, k].mean(), null[:, k].std()
        d = dict(obs=float(obs[k]), null_mean=float(m), null_sd=float(s),
                 z=float((obs[k] - m) / s) if s > 0 else float("nan"))
        if nm == "rho":
            d["p_one_sided"], d["p_two_sided"] = [float(v) for v in
                                                   pvals(obs[k], null[:, k], two_sided=True)]
        else:
            d["p_upper"] = float(pvals(obs[k], null[:, k]))
        out[nm] = d
    nb = c.size
    idx = rng.integers(0, nb, (n_boot, nb))
    boots = c.ravel()[idx].sum(1) / np.sqrt(sxx.ravel()[idx].sum(1) * syy.ravel()[idx].sum(1))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    out["ci95"] = (float(lo), float(hi))
    out["n_blocks"] = int(nb)
    out["shape"] = (int(H), int(W))
    out["u_blocks"] = u
    return out


def sign_test(p, x, p_thr, x_thr, n_null, rng):
    """符号一致率。x 側のヌルは AAFT(周辺分布を保つ)。thr は値(または 'median')。"""
    pt = np.median(p) if p_thr == "median" else float(p_thr)
    xt = np.median(x) if x_thr == "median" else float(x_thr)

    def stat(xv):
        pv, xs = p.ravel(), xv.ravel()
        m = (pv != pt) & (xs != xt)
        sp, sx = pv[m] > pt, xs[m] > xt
        a = (sp == sx).mean()
        fp, fx = sp.mean(), sx.mean()
        e = fp * fx + (1 - fp) * (1 - fx)
        kappa = (a - e) / (1 - e) if e < 1 else np.nan
        tpr = (sp & sx).sum() / max(sx.sum(), 1)
        tnr = (~sp & ~sx).sum() / max((~sx).sum(), 1)
        return np.array([a, kappa, (tpr + tnr) / 2, fx, fp, m.mean()])

    obs = stat(x)
    null = np.array([stat(aaft_surrogate(x, rng)) for _ in range(n_null)])
    out = {}
    for k, nm in enumerate(["agree", "kappa", "balanced_acc", "fx", "fp", "frac_used"]):
        m, s = null[:, k].mean(), null[:, k].std()
        out[nm] = dict(obs=float(obs[k]), null_mean=float(m),
                       z=float((obs[k] - m) / s) if s > 0 else float("nan"),
                       p_upper=float(pvals(obs[k], null[:, k])))
    out["thresholds"] = dict(p=float(pt), x=float(xt))
    return out


# ----------------------------------------------------------------------------
# 任意の補足(記述的。判定には使わない)
# ----------------------------------------------------------------------------
def contrast_strata(p, x, controls, q, n_null, rng):
    """x の勾配の大きさが上位 (1-q) の画素(hi)とそれ以外(lo)での偏相関。マスクもヌル側で作り直す。"""
    resid = make_resid(controls)
    rP = resid(rank01(p))

    def stat(xv):
        gy, gx = np.gradient(xv.astype(float))
        g = np.hypot(gx, gy).ravel()
        hi = g >= np.quantile(g, q)
        rX = resid(rank01(xv))
        f = lambda m: float(np.corrcoef(rX[m], rP[m])[0, 1])
        return np.array([f(hi), f(~hi)])

    obs = stat(x)
    null = np.array([stat(phase_surrogate(x, rng)) for _ in range(n_null)])
    z = lambda o, nl: float((o - nl.mean()) / nl.std())
    d_obs, d_null = obs[0] - obs[1], null[:, 0] - null[:, 1]
    return dict(hi=(float(obs[0]), z(obs[0], null[:, 0])),
                lo=(float(obs[1]), z(obs[1], null[:, 1])),
                hi_minus_lo=(float(d_obs), z(d_obs, d_null)))


def coarse_test(p, x, controls, S, n_null, rng):
    """S×S ブロック平均した画像どうしの偏順位相関(「全体の質感」のスケール)。
    ヌルは全解像度の x の位相ランダム化 → 同じブロック平均。p・x・controls は周囲除外済みの 2 次元配列。"""
    def blk(a):
        H, W = a.shape
        hb, wb = H // S, W // S
        return a[: hb * S, : wb * S].reshape(hb, S, wb, S).mean((1, 3))
    P, X, C = blk(p), blk(x), [blk(c) for c in controls]
    resid = make_resid(C)
    rP = resid(rank01(P))
    nP = np.linalg.norm(rP)

    def stat(xv):
        rX = resid(rank01(blk(xv)))
        return float(rX @ rP / (np.linalg.norm(rX) * nP))

    obs = stat(x)
    null = np.array([stat(phase_surrogate(x, rng)) for _ in range(n_null)])
    m, sd = null.mean(), null.std()
    up, two = pvals(obs, null, two_sided=True)
    return dict(S=int(S), n_blocks=int(P.size), rho=float(obs),
                z=float((obs - m) / sd) if sd > 0 else float("nan"),
                null_sd=float(sd), p_one_sided=float(up), p_two_sided=float(two))


def border_diag(P, X, controls, w):
    """偏相関への寄与を、外周(幅 w)と内部に分解する(合計 = margin 0 の偏相関)。
    あわせて、x・p・controls の外周の平均順位(0.5 なら偏りなし。外周が暗い・低いと 0.5 より小さい)を出す。"""
    H, W = X.shape
    mask = np.ones((H, W), bool)
    mask[w:-w, w:-w] = False
    resid = make_resid(controls)
    rP, rX = resid(rank01(P)), resid(rank01(X))
    tot = np.linalg.norm(rP) * np.linalg.norm(rX)
    prod = (rX * rP).reshape(H, W)

    def bm(a):
        return float(rank01(a).reshape(a.shape)[mask].mean())

    ranks = dict(x=bm(X), p=bm(P))
    for nm, c in zip(("avg", "last"), controls):
        ranks[nm] = bm(c)
    return dict(width=int(w), area_share=float(mask.mean()),
                contrib_border=float(prod[mask].sum() / tot),
                contrib_interior=float(prod[~mask].sum() / tot),
                mean_rank_border=ranks)


def band_test(p, x, controls, margin, B, sigma, which, n_null, rng):
    def f(a):
        low = gaussian_filter(a.astype(float), sigma, mode="reflect")
        return prep(low if which == "low" else a - low, margin, B)
    r = core_test(f(p), f(x), [f(c) for c in controls], B, n_null, rng, n_boot=200)
    return dict(rho=r["rho"]["obs"], z=r["rho"]["z"])


def calibration_check(p, x, controls, margin, B, trials, n_null, rng, obs_rho=None):
    """x を大きく巡回シフトして p と無関係にし(疑似陰性)、主指標の z が N(0,1) に近いかを見る。
    あわせて、疑似陰性の rho の分布(ずらしヌル)に対する観測 rho の両側 p も出す。"""
    H, W = x.shape
    zs, rhos, shifts = [], [], set()
    for _ in range(trials):
        dy = int(rng.integers(H // 4, 3 * H // 4 + 1))
        dx = int(rng.integers(W // 4, 3 * W // 4 + 1))
        shifts.add((dy, dx))
        xs = np.roll(x, (dy, dx), (0, 1))
        r = core_test(prep(p, margin, B), prep(xs, margin, B),
                      [prep(c, margin, B) for c in controls], B, n_null, rng, n_boot=10)
        zs.append(r["rho"]["z"])
        rhos.append(r["rho"]["obs"])
    zs, rhos = np.array(zs), np.array(rhos)
    out = dict(trials=int(trials), z_mean=float(zs.mean()), z_sd=float(zs.std(ddof=1)),
               frac_z_gt_1_645=float((zs > 1.645).mean()),
               frac_abs_z_gt_1_96=float((np.abs(zs) > 1.96).mean()),
               rho_mean=float(rhos.mean()), rho_sd=float(rhos.std(ddof=1)),
               n_distinct_shifts=len(shifts),
               note="理想: z_mean≈0, z_sd≈1, frac(z>1.645)≈0.05, frac(|z|>1.96)≈0.05。z_sd が 1 を大きく超えるとヌルの分散が過小。"
                    "z_mean が 0 から大きくずれると、ずらしでも p と相関が出る構造がある")
    if obs_rho is not None:
        out["shift_null_p_two_sided"] = float((1 + (np.abs(rhos) >= abs(obs_rho)).sum()) / (1 + len(rhos)))
        out["shift_null_note"] = "ずらしヌルの p の実質的な下限は 1/(1+ずらしの種類数)。種類が少ないと粗い目安"
    return out


# ----------------------------------------------------------------------------
# 健全性チェック(取り違え・リークの検出)
# ----------------------------------------------------------------------------
def sanity(p, x, avg, last):
    errs, warns = [], []
    arrs = dict(p=p, x=x, avg=avg, last=last)
    for k, a in arrs.items():
        if a is not None and a.shape != x.shape:
            errs.append("形が違う: %s %s vs x %s" % (k, a.shape, x.shape))
    if errs:
        return errs, warns
    for k, a in arrs.items():
        if k != "x" and a is not None and np.array_equal(a, x):
            errs.append("%s が x と完全に一致(test と prediction の取り違え・リークの疑い)" % k)
    if np.array_equal(p, last if last is not None else p + 1):
        warns.append("p が last と完全に一致")
    if np.array_equal(p, avg):
        warns.append("p が avg と完全に一致")
    for k, a in (("p", p), ("avg", avg), ("last", last)):
        if a is not None and np.ptp(a) > 0:
            r = spearman(a, x)
            if r > 0.9:
                warns.append("rho(%s, x)=%.3f が非常に高い(x の混入・リークの疑い)" % (k, r))
    return errs, warns


# ----------------------------------------------------------------------------
# 1 画像の解析
# ----------------------------------------------------------------------------
def check_order(idx):
    """フレーム番号が与えられていれば、入力が対象より前であることを検査する(1 枚後ろの取り違え対策)。"""
    errs, warns = [], []
    xi, li, am = idx.get("x_idx"), idx.get("last_idx"), idx.get("avg_max_idx")
    if xi is not None and li is not None:
        if li >= xi:
            errs.append("last_idx(%d) が x_idx(%d) 以降(対象より後ろの画像を入力に使っている)" % (li, xi))
        elif li != xi - 1:
            warns.append("last_idx(%d) が x_idx-1(%d) ではない" % (li, xi - 1))
    if xi is not None and am is not None and am >= xi:
        errs.append("avg_max_idx(%d) が x_idx(%d) 以降(avg に対象以降の画像が混入)" % (am, xi))
    return errs, warns


def analyze(name, p, x, avg, last, cfg, seed_seq, paths=None, idx=None, log=print):
    rng = np.random.default_rng(seed_seq)
    B, m0 = cfg["block"], cfg["margin"]
    controls_full = [avg] + ([last] if last is not None else [])
    res = dict(name=name, shape=tuple(x.shape),
               controls="avg+last" if last is not None else "avg のみ")
    if paths:
        res["sha256"] = {k: sha256_of(v) for k, v in paths.items()}
    errs, warns = sanity(p, x, avg, last)
    if idx:
        e2, w2 = check_order(idx)
        errs, warns = errs + e2, warns + w2
        res["frame_idx"] = idx
    res["sanity_errors"], res["sanity_warnings"] = errs, warns
    if errs:
        return res

    P, X = prep(p, m0, B), prep(x, m0, B)
    C = [prep(c, m0, B) for c in controls_full]
    if min(X.shape) < 3 * B:
        res["sanity_errors"] = ["周囲除外後の領域が小さすぎる %s" % (X.shape,)]
        return res

    ctl = dict(rho_avg_x=spearman(C[0], X))
    if last is not None:
        rng_x = float(np.ptp(X)) or 1.0
        ctl["rho_last_x"] = spearman(C[1], X)
        ctl["frac_last_identical"] = float(np.mean(C[1] == X))
        ctl["frac_last_near_identical"] = float(np.mean(np.abs(C[1] - X) <= 0.01 * rng_x))
        if ctl["frac_last_near_identical"] > 0.3:
            res["sanity_warnings"].append(
                "last が x とほぼ同一の画素が %.0f%%。差し引きで残差が差の部分に集中し、偏相関が不安定・過大になりうる(例外的な入力を外して再計算を推奨)"
                % (100 * ctl["frac_last_near_identical"]))
    res["control_strength"] = ctl
    thr_c = cfg.get("max_control_rho", 0.85)
    reasons = []
    if abs(ctl["rho_avg_x"]) > thr_c:
        reasons.append("|rho(avg,x)|=%.2f > %.2f" % (abs(ctl["rho_avg_x"]), thr_c))
    if "rho_last_x" in ctl and abs(ctl["rho_last_x"]) > thr_c:
        reasons.append("|rho(last,x)|=%.2f > %.2f" % (abs(ctl["rho_last_x"]), thr_c))
    if reasons:
        res["excluded_reason"] = ("controls が x をほとんど取り尽くしている(偏相関が残差のノイズになり解釈できない): "
                                  + "; ".join(reasons))
    t0 = time.time()
    res["PRIMARY"] = core_test(P, X, C, B, cfg["n_null"], rng, cfg["n_boot"])
    wb = cfg.get("border_width", 4)
    P0, X0 = prep(p, 0, B), prep(x, 0, B)
    if wb and min(X0.shape) >= 2 * wb + 3:
        res["border_diag"] = border_diag(P0, X0, [prep(c, 0, B) for c in controls_full], wb)
    res["raw_spearman"] = {k: v for k, v in
                           core_test(P, X, [], B, min(cfg["n_null"], 1000), rng, 200).items()
                           if k in ("rho", "ci95")}
    res["sign"] = sign_test(P, X, cfg["p_thr"], cfg["x_thr"], min(cfg["n_null"], 1000), rng)
    S = cfg.get("coarse", 8)
    if S and min(X.shape) >= 3 * S:
        res["coarse"] = coarse_test(P, X, C, S, min(cfg["n_null"], 1000), rng)

    edge = {}
    for m in cfg["margins"]:
        Pm, Xm = prep(p, m, B), prep(x, m, B)
        if min(Xm.shape) < 3 * B:
            continue
        r = core_test(Pm, Xm, [prep(c, m, B) for c in controls_full], B,
                      min(cfg["n_null"], 1000), rng, 500)
        edge["margin=%d" % m] = dict(rho=r["rho"]["obs"], z=r["rho"]["z"], ci95=r["ci95"])
    res["edge_sensitivity"] = edge

    if cfg.get("extras"):
        n_ex = min(cfg["n_null"], 500)
        ex = {"strata_q50": contrast_strata(P, X, C, 0.5, n_ex, rng),
              "strata_q80": contrast_strata(P, X, C, 0.8, n_ex, rng)}
        for s in (2, 4, 8):
            for w in ("low", "high"):
                ex["band_sigma%d_%s" % (s, w)] = band_test(
                    p, x, controls_full, m0, B, s, w, n_ex, rng)
        res["extras"] = ex
    if cfg.get("calibrate"):
        res["calibration"] = calibration_check(p, x, controls_full, m0, B,
                                               cfg["calibrate"], 200, rng,
                                               obs_rho=res["PRIMARY"]["rho"]["obs"])
    cal = res.get("calibration")
    if cal and (cal["frac_z_gt_1_645"] > 0.15 or cal["frac_abs_z_gt_1_96"] > 0.15):
        why = "校正不良(疑似陰性で偽陽性率が高い: z>1.645 %.0f%%, |z|>1.96 %.0f%%)" % (
            100 * cal["frac_z_gt_1_645"], 100 * cal["frac_abs_z_gt_1_96"])
        res["excluded_reason"] = (res["excluded_reason"] + "; " + why) if res.get("excluded_reason") else why
    res["seconds"] = round(time.time() - t0, 1)
    return res


def verdict(pr, two_sided=False):
    if two_sided:
        z = abs(pr["rho"]["z"])
        sg = "正" if pr["rho"]["obs"] >= 0 else "負"
        if z >= 2.0:
            return "有意(|z|≥2, 向きは %s。向きは任意なので x を見ずに固定した規則がない限り符号は解釈しない)" % sg
        if z >= 1.645:
            return "境界的(|z|≥1.645, 向き %s)。別の x で再現するかを見る" % sg
        return "有意でない(|z|<1.645)"
    r, z, lo = pr["rho"]["obs"], pr["rho"]["z"], pr["ci95"][0]
    if z >= 2.0 and lo > 0:
        return "有意(z≥2 かつ区間が 0 を含まない)。ただし 1 画像。別の対象画像での再現が必要"
    if z >= 1.645:
        return "境界的(片側 5% 水準)。区間が 0 を含むかを確認し、画像を増やす"
    if r > 0:
        return "点推定は正だが有意でない"
    return "正の関連は見えない"


def fmt_item(res):
    L = []
    L.append("=" * 78)
    L.append("[%s] shape=%s controls=%s" % (res["name"], res["shape"], res.get("controls")))
    for e in res.get("sanity_errors", []):
        L.append("  !! エラー: " + e)
    for w in res.get("sanity_warnings", []):
        L.append("  ?? 警告: " + w)
    if res.get("excluded_reason"):
        L.append("  !! 集計から除外: " + res["excluded_reason"])
    if "PRIMARY" not in res:
        return "\n".join(L)
    pr = res["PRIMARY"]
    r = pr["rho"]
    L.append("  PRIMARY 偏相関(avg,last 差し引き, 周囲除外): rho=%.4f  z=%.2f  p_片側=%.4f  p_両側=%.4f"
             % (r["obs"], r["z"], r["p_one_sided"], r["p_two_sided"]))
    L.append("          null_sd=%.3f → この画像で検出できる最小の効果(検出力80%%・片側2.5%%)≈ %.2f  ※これより小さい効果は見えない"
             % (r["null_sd"], 2.8 * r["null_sd"]))
    L.append("          区間95%%(ブロックBS, %d ブロック; 滑らかな画像では狭すぎる傾向があり参考値)=[%.3f, %.3f]"
             % (pr["n_blocks"], pr["ci95"][0], pr["ci95"][1]))
    L.append("          判定: " + verdict(pr, two_sided=bool(res.get("ensemble"))))
    if res.get("ensemble"):
        en = res["ensemble"]
        L.append("          ensemble: 同じ x に対する %d 個の予測を、p どうしの符号(x は使わない)で揃えて平均。"
                 "負に揃えた予測=%d  p-p 相関の平均 |ρ|=%.2f  実効的に独立な予測数≈%.1f"
                 % (en["n_members"], en["n_flipped"], en["mean_abs_pp_corr"], en["n_eff_members"]))
        L.append("          p-p 相関行列の固有値の寄与率(上位): " +
                 " ".join("%.2f" % v for v in en["eig_share_top5"]) +
                 "  ※実効数は参加比 (Σλ)²/Σλ² で、階数(次元)ではない。第1成分が大きいと次元が大きくても小さくなる")
        L.append("          基準にした予測=%s  向きを反転した予測=%s"
                 % (en.get("reference_member"), ",".join(en.get("flipped_members", [])) or "なし"))
        if en.get("warn"):
            L.append("          ?? " + en["warn"])
    bd = res.get("border_diag")
    if bd:
        mr = bd["mean_rank_border"]
        L.append("  border(幅 %d, 面積 %.0f%%): 偏相関への寄与 外周 %+.3f / 内部 %+.3f   外周の平均順位(0.5=偏りなし) x=%.2f p=%.2f%s"
                 % (bd["width"], 100 * bd["area_share"], bd["contrib_border"], bd["contrib_interior"],
                    mr["x"], mr["p"],
                    "".join(" %s=%.2f" % (k, mr[k]) for k in ("avg", "last") if k in mr)))
    ctl = res.get("control_strength", {})
    if ctl:
        L.append("  controls と x の関係: rho(avg,x)=%.2f  rho(last,x)=%s  last が x と(ほぼ)同一の画素の割合=%s"
                 % (ctl["rho_avg_x"],
                    "%.2f" % ctl["rho_last_x"] if "rho_last_x" in ctl else "-",
                    "%.0f%%(完全一致 %.0f%%)" % (100 * ctl["frac_last_near_identical"], 100 * ctl["frac_last_identical"])
                    if "frac_last_near_identical" in ctl else "-"))
    rs = res["raw_spearman"]
    L.append("  raw Spearman(差し引きなし): rho=%.4f  区間=[%.3f, %.3f]"
             % (rs["rho"]["obs"], rs["ci95"][0], rs["ci95"][1]))
    sg = res["sign"]
    L.append("  sign(中央値で2分割): agree=%.4f(偶然 0.5 付近: null %.4f) kappa=%.3f bal_acc=%.3f  z=%.2f  p_upper=%.4f  fx=%.3f fp=%.3f"
             % (sg["agree"]["obs"], sg["agree"]["null_mean"], sg["kappa"]["obs"],
                sg["balanced_acc"]["obs"], sg["agree"]["z"], sg["agree"]["p_upper"],
                sg["fx"]["obs"], sg["fp"]["obs"]))
    if res.get("coarse"):
        cs = res["coarse"]
        L.append("  coarse(%d×%d ブロック平均=「全体の質感」のスケール, %d ブロック): rho=%.4f  z=%.2f  p_両側=%.4f  null_sd=%.3f"
                 % (cs["S"], cs["S"], cs["n_blocks"], cs["rho"], cs["z"], cs["p_two_sided"], cs["null_sd"]))
    up, un, mr = pr["u_pos"], pr["u_neg"], pr["mean_rb2"]
    ratio = up["obs"] / -un["obs"] if un["obs"] < 0 else float("inf")
    L.append("  diag: u_pos=%.3f(z %.2f)  u_neg=%.3f(z %.2f; 正なら負の寄与が抑制)  pos/|neg|=%.2f  mean_rb2=%.3f(z %.2f)"
             % (up["obs"], up["z"], un["obs"], un["z"], ratio, mr["obs"], mr["z"]))
    if mr["z"] > 5 and un["z"] < -2:
        L.append("        → 局所的に強い依存があるが符号が領域ごとに混在し、全体で相殺している型")
    if res.get("edge_sensitivity"):
        L.append("  edge: " + " | ".join("%s rho=%.3f z=%.2f" % (k, v["rho"], v["z"])
                                          for k, v in res["edge_sensitivity"].items()))
    if "extras" in res:
        ex = res["extras"]
        for k in ("strata_q50", "strata_q80"):
            s = ex[k]
            L.append("  %s: hi=%.3f(z %.2f) lo=%.3f(z %.2f) hi-lo=%.3f(z %.2f)"
                     % (k, s["hi"][0], s["hi"][1], s["lo"][0], s["lo"][1],
                        s["hi_minus_lo"][0], s["hi_minus_lo"][1]))
        L.append("  bands: " + " ".join("%s=%.3f(z %.2f)" % (k.replace("band_", ""), v["rho"], v["z"])
                                         for k, v in ex.items() if k.startswith("band_")))
    if "diversity_curve" in res:
        dc = res["diversity_curve"]
        L.append("  diversity(同じ x の中でメンバーの部分集合を取った記述的な曲線; k=メンバー数):")
        for k, v in dc["by_k"].items():
            L.append("     k=%-2d  平均 n_eff=%.2f  平均|z|=%.2f  平均|rho|=%.3f" % (k, v["mean_n_eff"], v["mean_abs_z"], v["mean_abs_rho"]))
        L.append("     n_eff と |z| の順位相関: 全体 %+.2f / k を固定した中の平均 %+.2f (permutation p=%.3f, 部分集合は非独立なので甘い目安)"
                 % (dc["spearman_neff_absz_overall"], dc["within_k_spearman"], dc["within_k_perm_p_upper"]))
    if "calibration" in res:
        c = res["calibration"]
        L.append("  calibration(疑似陰性 %d 回): z_mean=%.2f z_sd=%.2f frac(z>1.645)=%.3f frac(|z|>1.96)=%.3f  (理想 0 / 1 / 0.05 / 0.05)"
                 % (c["trials"], c["z_mean"], c["z_sd"], c["frac_z_gt_1_645"], c["frac_abs_z_gt_1_96"]))
        if "shift_null_p_two_sided" in c:
            L.append("              疑似陰性の rho: 平均 %+.3f sd %.3f(ずらしの種類 %d)  観測 rho に対するずらしヌルの両側 p=%.3f(下限は粗い)"
                     % (c["rho_mean"], c["rho_sd"], c["n_distinct_shifts"], c["shift_null_p_two_sided"]))
    L.append("  (%s 秒)" % res.get("seconds"))
    return "\n".join(L)


# ----------------------------------------------------------------------------
# 複数画像の集計
# ----------------------------------------------------------------------------
def aggregate(results, seed=0, n_boot=5000):
    ok = [r for r in results if "PRIMARY" in r and not r.get("excluded_reason")]
    if not ok:
        return None
    rho = np.array([r["PRIMARY"]["rho"]["obs"] for r in ok])
    z = np.array([r["PRIMARY"]["rho"]["z"] for r in ok])
    rng = np.random.default_rng(seed)
    boots = [rng.choice(rho, len(rho)).mean() for _ in range(n_boot)]
    nsd = np.array([r["PRIMARY"]["rho"]["null_sd"] for r in ok])
    pooled_sd = float(np.sqrt((nsd ** 2).sum()) / len(nsd))       # 画像が独立なら
    from scipy.stats import chi2
    sz2 = float((z ** 2).sum())
    ne = [(r["ensemble"]["n_eff_members"], abs(r["PRIMARY"]["rho"]["z"]))
          for r in ok if r.get("ensemble")]
    neff_corr = None
    if len(ne) >= 3:
        a_, b_ = np.array(ne).T
        rs = spearman(a_, b_)
        perm = [spearman(a_, np.random.default_rng(i).permutation(b_)) for i in range(2000)]
        neff_corr = dict(n=len(ne), spearman=rs,
                         perm_p_two_sided=float((1 + sum(abs(t) >= abs(rs) for t in perm)) / 2001),
                         points=[(float(u), float(v)) for u, v in ne],
                         note="x の数が少ないと有意になり得ない(n=3 の Spearman の両側 p の下限は 1/3)")
    groups = {}
    for r in ok:
        h = r.get("sha256", {})
        if h.get("x"):
            groups.setdefault(("px", h.get("p"), h.get("x")), []).append(r["name"])
            groups.setdefault(("x", h.get("x")), []).append(r["name"])
    dup_px = [v for k, v in groups.items() if k[0] == "px" and len(v) > 1]
    same_x = [v for k, v in groups.items() if k[0] == "x" and len(v) > 1]
    n_distinct_x = len({k[1] for k in groups if k[0] == "x"}) if groups else len(ok)
    k_eff, p_adj = None, None
    U = [np.asarray(r["PRIMARY"]["u_blocks"]).ravel() for r in ok]
    if len(U) > 1 and len({u.size for u in U}) == 1:
        Rm = np.corrcoef(np.stack(U))
        if not np.isnan(Rm).any():
            evm = np.clip(np.linalg.eigvalsh(Rm), 0, None)
            k_eff = float(evm.sum() ** 2 / (evm ** 2).sum())
            p_adj = float(chi2.sf(sz2 * k_eff / len(z), k_eff))
    return dict(n_images=len(ok), n_distinct_x=int(n_distinct_x), mean_rho=float(rho.mean()),
                n_negative=int((rho < 0).sum()), mean_abs_rho=float(np.abs(rho).mean()),
                sum_z2=sz2, chi2_p_if_independent=float(chi2.sf(sz2, len(z))),
                duplicate_p_x=dup_px, same_x_groups=same_x,
                n_eff_targets_from_u_maps=k_eff, chi2_p_adjusted_by_keff=p_adj,
                neff_members_vs_absz=neff_corr,
                approx_mde_80pct_if_independent=float(2.8 * pooled_sd),
                ci95_over_images=[float(v) for v in np.percentile(boots, [2.5, 97.5])],
                n_positive=int((rho > 0).sum()),
                stouffer_z=float(z.sum() / np.sqrt(len(z))),
                note="Stouffer は画像が独立(対象が別・入力が重ならない)な場合のみ有効。"
                     "同じ x の別予測や重なる入力では過大になるので区間(ci95_over_images)を主に見る")


# ----------------------------------------------------------------------------
# 同じ x に対する複数の予測を束ねる(x は使わずに p どうしの符号で揃える)
# ----------------------------------------------------------------------------
def build_ensembles(items):
    """items: [(name, paths, idx), ...] → [(name, p_ens, x, avg, last, paths, idx), ...]
    x の sha256 でグループ化し、各グループで最初の p を基準に p どうしの Spearman の符号で向きを揃え、
    順位を平均する。向きの決定に x を使わないのでリークしない(ただし絶対的な向きは任意)。"""
    groups = {}
    for name, paths, idx in items:
        groups.setdefault(sha256_of(paths["x"]), []).append((name, paths, idx))
    out, info, data = [], {}, {}
    for h, grp in groups.items():
        arrs = [{k: load_array(v) for k, v in paths.items()} for _, paths, _ in grp]
        ref = arrs[0]["p"]
        aligned, flipped, flipped_names = [], 0, []
        for (mname, _, _), a_ in zip(grp, arrs):
            rk = rank01(a_["p"]).reshape(a_["p"].shape)
            if spearman(a_["p"], ref) < 0:
                rk, flipped = 1.0 - rk, flipped + 1
                flipped_names.append(mname)
            aligned.append(rk)
        p_ens = np.mean(aligned, axis=0)
        k = len(aligned)
        if k > 1:
            R = np.corrcoef(np.stack([a.ravel() for a in aligned]))
            off = R[np.triu_indices(k, 1)]
            ev = np.clip(np.linalg.eigvalsh(R), 0, None)
            n_eff = float(ev.sum() ** 2 / (ev ** 2).sum())
            mean_abs = float(np.abs(off).mean())
            share = (np.sort(ev)[::-1] / ev.sum())[:5]
        else:
            n_eff, mean_abs, share = 1.0, float("nan"), np.array([1.0])
        warn = None
        base = arrs[0]
        for a_ in arrs[1:]:
            if not (np.array_equal(a_["avg"], base["avg"]) and
                    np.array_equal(a_.get("last", base.get("last")), base.get("last"))):
                warn = "グループ内で avg / last が一致しない。最初の項目のものを使用"
        name = "ens[x=%s, n=%d]" % (h[:8], k)
        meta = {kk: grp[0][1][kk] for kk in ("x", "avg", "last") if kk in grp[0][1]}
        out.append((name, p_ens, base["x"], base["avg"], base.get("last"), meta, grp[0][2]))
        data[name] = dict(aligned=np.stack(aligned), x=base["x"], avg=base["avg"], last=base.get("last"))
        info[name] = dict(n_members=k, n_flipped=int(flipped), flipped_members=flipped_names,
                          reference_member=grp[0][0], mean_abs_pp_corr=mean_abs,
                          n_eff_members=n_eff, eig_share_top5=[float(v) for v in share],
                          members=[g[0] for g in grp], warn=warn)
    return out, info, data


# ----------------------------------------------------------------------------
# 予測の多様性(実効的に独立な予測数 n_eff)と |z| の関係(記述的)
# ----------------------------------------------------------------------------
def effective_number(stack):
    """stack: (k, H, W) の向き揃え済み順位画像。相関行列の参加比 (Σλ)²/Σλ²。"""
    k = stack.shape[0]
    if k < 2:
        return 1.0
    ev = np.clip(np.linalg.eigvalsh(np.corrcoef(stack.reshape(k, -1))), 0, None)
    return float(ev.sum() ** 2 / (ev ** 2).sum())


def diversity_curve(d, cfg, rng, sizes, n_subsets, n_null):
    """同じ x の中で、メンバーの部分集合をランダムに取り、(部分集合の n_eff, 束ねた予測の |z|) を集める。
    k(メンバー数)を固定した中での n_eff と |z| の順位相関(k の効果を除く)も出す。
    部分集合どうしは同じメンバー・同じ x を共有するので、p は記述的な目安(甘め)。"""
    B, m0 = cfg["block"], cfg["margin"]
    A, x, avg, last = d["aligned"], d["x"], d["avg"], d["last"]
    controls = [avg] + ([last] if last is not None else [])
    K = A.shape[0]
    X = prep(x, m0, B)
    C = [prep(c, m0, B) for c in controls]
    rows = []
    for k in sizes:
        if k > K:
            continue
        for _ in range(n_subsets):
            idx = rng.choice(K, k, replace=False)
            sub = A[idx]
            r = core_test(prep(sub.mean(0), m0, B), X, C, B, n_null, rng, n_boot=10)
            rows.append((k, effective_number(sub), r["rho"]["obs"], abs(r["rho"]["z"])))
    rows = np.array(rows)
    ks = sorted({int(v) for v in rows[:, 0]})
    by_k = {int(k): dict(mean_n_eff=float(rows[rows[:, 0] == k, 1].mean()),
                         mean_abs_z=float(rows[rows[:, 0] == k, 3].mean()),
                         mean_abs_rho=float(np.abs(rows[rows[:, 0] == k, 2]).mean())) for k in ks}

    def within_k(z):
        cs = []
        for k in ks:
            m = rows[:, 0] == k
            if k > 1 and m.sum() >= 5 and np.ptp(rows[m, 1]) > 0 and np.ptp(z[m]) > 0:
                cs.append(spearman(rows[m, 1], z[m]))
        return float(np.mean(cs)) if cs else float("nan")

    z = rows[:, 3]
    obs = within_k(z)
    null = []
    for _ in range(1000):
        zp = z.copy()
        for k in ks:
            m = np.where(rows[:, 0] == k)[0]
            zp[m] = rng.permutation(z[m])
        null.append(within_k(zp))
    null = np.array([v for v in null if not np.isnan(v)])
    p_perm = float((1 + (null >= obs).sum()) / (1 + len(null))) if len(null) and not np.isnan(obs) else float("nan")
    return dict(by_k=by_k, spearman_neff_absz_overall=spearman(rows[:, 1], z),
                within_k_spearman=obs, within_k_perm_p_upper=p_perm,
                n_points=int(len(rows)),
                note="記述的。部分集合は同じメンバー・同じ x を共有するため p は甘い。k を固定した相関が本命")


# ----------------------------------------------------------------------------
# 合成データによる動作確認(陽性・陰性対照)
# ----------------------------------------------------------------------------
def synth_items(rng, H=48, W=32):
    def field(s):
        f = gaussian_filter(rng.standard_normal((H, W)), s, mode="wrap")
        return (f - f.mean()) / f.std()
    f1, f2, f3, f4, f5 = (field(2.5) for _ in range(5))
    x = 128 + 40 * f1
    last = 128 + 40 * (0.3 * f1 + 0.95 * f2)
    avg = 128 + 40 * (0.2 * f1 + 0.98 * f3)
    p_pos = 84 + 5 * (0.7 * f1 + 0.71 * f4)       # x と相関 0.7 の予測(陽性対照)
    p_neg = 84 + 5 * f5                            # x と独立(陰性対照)
    return x, avg, last, p_pos, p_neg


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--p", "--x", "--avg", "--last"):
        ap.add_argument(k)
    ap.add_argument("--manifest", help="CSV(列: name,p,x,avg,last)")
    ap.add_argument("--out", help="出力の接頭辞(.json と .txt を書く)")
    ap.add_argument("--margin", type=int, default=DEFAULTS["margin"])
    ap.add_argument("--block", type=int, default=DEFAULTS["block"])
    ap.add_argument("--n-null", type=int, default=DEFAULTS["n_null"])
    ap.add_argument("--n-boot", type=int, default=DEFAULTS["n_boot"])
    ap.add_argument("--seed", type=int, default=DEFAULTS["seed"])
    ap.add_argument("--p-thr", default=DEFAULTS["p_thr"], help="p の符号の閾値(数値 or median)")
    ap.add_argument("--x-thr", default=DEFAULTS["x_thr"], help="x の符号の閾値(数値 or median)")
    ap.add_argument("--extras", action="store_true", help="コントラスト層別・周波数帯(記述的)")
    ap.add_argument("--calibrate", type=int, default=0, help="疑似陰性の回数(例 100)")
    ap.add_argument("--selftest", action="store_true", help="合成データでの動作確認")
    ap.add_argument("--coarse", type=int, default=8, help="coarse 指標のブロック幅(0 で無効)")
    ap.add_argument("--max-control-rho", type=float, default=0.85,
                    help="|rho(avg,x)| または |rho(last,x)| がこれを超えたら集計から除外(controls が x を取り尽くしている)")
    ap.add_argument("--border-width", type=int, default=4,
                    help="border 診断の外周の幅(寄与の分解。0 で無効)")
    ap.add_argument("--diversity-curve", action="store_true",
                    help="(--ensemble-by-x と併用) 同じ x 内でメンバーの部分集合を取り、実効予測数 n_eff と |z| の関係を見る")
    ap.add_argument("--subsets", type=int, default=15, help="diversity-curve: 各サイズの部分集合の数")
    ap.add_argument("--ensemble-by-x", action="store_true",
                    help="同じ x に対する複数の予測を、p どうしの符号(x は使わない)で揃えて平均し、x ごとに 1 つにする")
    ap.add_argument("--x-idx", type=int, help="対象 x のフレーム番号(任意。指定すると時間順を検査)")
    ap.add_argument("--last-idx", type=int, help="last のフレーム番号")
    ap.add_argument("--avg-max-idx", type=int, help="avg に使った最大のフレーム番号")
    a = ap.parse_args()

    cfg = dict(DEFAULTS, margin=a.margin, block=a.block, n_null=a.n_null, n_boot=a.n_boot,
               seed=a.seed, p_thr=a.p_thr, x_thr=a.x_thr, extras=a.extras, calibrate=a.calibrate,
               coarse=a.coarse, border_width=a.border_width,
               max_control_rho=a.max_control_rho)
    txt = open(a.out + ".txt", "w", encoding="utf-8") if a.out else None

    def log(s=""):
        print(s, flush=True)
        if txt:
            txt.write(s + "\n")
            txt.flush()

    meta = dict(version=VERSION, script_sha256=sha256_of(os.path.abspath(__file__)),
                argv=sys.argv, python=platform.python_version(), numpy=np.__version__,
                platform=platform.platform(), params={k: v for k, v in cfg.items()},
                time=time.strftime("%Y-%m-%d %H:%M:%S"))
    try:
        import scipy
        meta["scipy"] = scipy.__version__
    except Exception:
        pass
    items = []
    if a.selftest:
        rng = np.random.default_rng(12345)
        x, avg, last, p_pos, p_neg = synth_items(rng)
        cfg["n_null"] = min(cfg["n_null"], 500)
        cfg["calibrate"] = cfg["calibrate"] or 40
        items = [("selftest_positive(期待: 有意)", p_pos, x, avg, last, None),
                 ("selftest_negative(期待: z≈0)", p_neg, x, avg, last, None)]
    elif a.manifest:
        with open(a.manifest, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                paths = {k: row[k] for k in ("p", "x", "avg", "last") if row.get(k)}
                idx = {k: int(row[k]) for k in ("x_idx", "last_idx", "avg_max_idx")
                       if row.get(k) not in (None, "")}
                items.append((row.get("name") or paths["p"], paths, idx))
    elif a.p and a.x and a.avg:
        paths = {k: getattr(a, k) for k in ("p", "x", "avg", "last") if getattr(a, k)}
        idx = {k: getattr(a, k) for k in ("x_idx", "last_idx", "avg_max_idx")
               if getattr(a, k) is not None}
        items.append((os.path.basename(a.p), paths, idx))
    else:
        ap.error("--p --x --avg(任意 --last)、--manifest、--selftest のいずれかを指定してください")

    meta["params"] = {k: v for k, v in cfg.items()}
    log("%s  script_sha256=%s" % (VERSION, meta["script_sha256"][:16]))
    log("predictor_check  params=%s" % json.dumps(jsonable(meta["params"]), ensure_ascii=False))
    ens_info, ens_data = {}, {}
    if a.ensemble_by_x and a.manifest and not a.selftest:
        items, ens_info, ens_data = build_ensembles(items)
        log("ensemble-by-x: %d 個の x に束ねた" % len(items))
    results = []
    for i, it in enumerate(items):
        try:
            if len(it) == 7:
                name, p, x, avg, last, paths, idx = it
            elif a.selftest:
                name, p, x, avg, last, paths = it
            else:
                name, paths, idx = it
                arr = {k: load_array(v) for k, v in paths.items()}
                p, x, avg, last = arr["p"], arr["x"], arr["avg"], arr.get("last")
            res = analyze(name, p, x, avg, last, cfg, [cfg["seed"], i], paths=paths,
                          idx=None if a.selftest else idx)
            if name in ens_info:
                res["ensemble"] = ens_info[name]
                if a.diversity_curve and "PRIMARY" in res:
                    K = ens_data[name]["aligned"].shape[0]
                    sizes = [k for k in (1, 2, 3, 5, 8, 12) if k <= K]
                    res["diversity_curve"] = diversity_curve(
                        ens_data[name], cfg, np.random.default_rng([cfg["seed"], i, 99]),
                        sizes, a.subsets, min(cfg["n_null"], 200))
        except Exception:
            res = dict(name=str(it[0]), sanity_errors=["例外: " + traceback.format_exc()])
        results.append(res)
        log(fmt_item(res))

    agg = aggregate(results, cfg["seed"])
    excl = [(r["name"], r["excluded_reason"]) for r in results if r.get("excluded_reason")]
    if excl:
        log("=" * 78)
        log("集計から除外した項目(%d 個): %s" % (len(excl), ", ".join(n for n, _ in excl)))
    if agg and agg["n_images"] > 1 and agg["n_distinct_x"] < agg["n_images"]:
        log("=" * 78)
        log("!! 異なる x は %d 種類だけ(%d 個の予測)。同じ x の別予測は独立な証拠ではないので、"
            "Stouffer・sum z^2・検出できる最小効果の数字は無効です。" % (agg["n_distinct_x"], agg["n_images"]))
        log("   x ごとに予測を束ねる(--ensemble-by-x)か、別の x を増やしてから集計してください。")
        log("   参考: mean_rho=%.4f  mean|rho|=%.3f  負の予測=%d / %d"
            % (agg["mean_rho"], agg["mean_abs_rho"], agg["n_negative"], agg["n_images"]))
    elif agg and agg["n_images"] > 1:
        log("=" * 78)
        log("集計(主指標): n=%d  mean_rho=%.4f  区間(画像BS)=[%.3f, %.3f]  正の画像数=%d  Stouffer z=%.2f"
            % (agg["n_images"], agg["mean_rho"], agg["ci95_over_images"][0],
               agg["ci95_over_images"][1], agg["n_positive"], agg["stouffer_z"]))
        log("  この枚数で検出できる最小の平均効果(独立と仮定、検出力80%%)≈ %.3f" % agg["approx_mde_80pct_if_independent"])
        log("  符号を見ない補助(探索的・事後導入): mean|rho|=%.3f  負の画像数=%d  sum z^2=%.1f (df=%d) p=%.2g(独立と仮定)"
            % (agg["mean_abs_rho"], agg["n_negative"], agg["sum_z2"], agg["n_images"],
               agg["chi2_p_if_independent"]))
        if agg["n_eff_targets_from_u_maps"] is not None:
            log("  寄与マップの相関から見積もった、実効的に独立な x の数≈%.2f(概算)。Σz² に k_eff/n を掛けて χ²(k_eff) と比べた p≈%.3f"
                % (agg["n_eff_targets_from_u_maps"], agg["chi2_p_adjusted_by_keff"]))
        if agg.get("neff_members_vs_absz"):
            c = agg["neff_members_vs_absz"]
            log("  x をまたいだ『束ねた予測の実効数 n_eff』と『|z|』の順位相関: %+.2f (n=%d, permutation 両側 p=%.3f)  点=%s"
                % (c["spearman"], c["n"], c["perm_p_two_sided"],
                   " ".join("(%.2f, %.2f)" % pt for pt in c["points"])))
        if agg["duplicate_p_x"]:
            log("  !! 同一の (p, x) の組が複数ある: %s" % agg["duplicate_p_x"])
        if agg["same_x_groups"]:
            log("  ?? 同じ x を共有する組(独立な画像ではない): %s" % agg["same_x_groups"])
        log("  注: " + agg["note"])
    if a.out:
        with open(a.out + ".json", "w", encoding="utf-8") as f:
            json.dump(jsonable(dict(meta=meta, items=results, aggregate=agg)), f,
                      ensure_ascii=False, indent=1)
        log("wrote %s.json / %s.txt" % (a.out, a.out))
    if txt:
        txt.close()


if __name__ == "__main__":
    main()
