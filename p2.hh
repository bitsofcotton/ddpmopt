/*
 BSD 3-Clause License

Copyright (c) 2022-2023, bitsofcotton (kazunobu watatsu)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#if !defined(_P2_)

using std::pair;
using std::make_pair;
using std::vector;
using std::move;
using std::cerr;
using std::endl;

// N.B. prediction on some range do persistent.
template <typename T> class P {
public:
  inline P() { ; }
  inline P(const int& status) {
    assert(0 < status);
    const int var1(max(T(int(2)), pow(T(status), T(int(1)) / T(int(3)))));
    const int var2(max(T(int(2)), pow(T(status), T(int(1)) / T(int(4)))));
    M.resize(7, Mx0 = MM = T(int(0)));;
    Mx.resize(M.size(), MM);
    p0 = P0maxRank<T>(status);
    p1 = shrinkMatrix<T, P1I<T, idFeeder<T> > >(P1I<T, idFeeder<T> >(status, var1, var1), var1);
    p2 = shrinkMatrix<T, P012L<T, idFeeder<T> > >(P012L<T, idFeeder<T> >(status, var2, var2), var2);
    q  = idFeeder<T>(status);
    q0 = SimpleVector<T>(status + 1).O();
    f  = idFeeder<T>(status);
  }
  inline ~P() { ; }
  inline const T& next(T d) {
    Mx0   = max(Mx0, abs(d));
    if(Mx0 == T(int(0)) || (d /= Mx0) == T(int(0)) || ! isfinite(d))
      return MM;
    M[1]  = p0.next(d);
    M[2]  = p1.next(d);
    M[3]  = p2.next(d);
    if((Mx[1] = max(Mx[1], abs(M[1]))) != T(int(0))) M[1] /= Mx[1];
    if((Mx[2] = max(Mx[2], abs(M[2]))) != T(int(0))) M[2] /= Mx[2];
    if((Mx[3] = max(Mx[3], abs(M[3]))) != T(int(0))) M[3] /= Mx[3];
    {
      auto qm(makeProgramInvariant<T>(q.next(d)));
      q0 += move(qm.first) * pow(qm.second, ceil(- log(SimpleMatrix<T>().epsilon())));
      auto qq(q);
      auto qqm(makeProgramInvariant<T>(qq.next(d)));
      M[4] = revertProgramInvariant<T>(make_pair(
        - (q0.dot(qqm.first) - q0[q0.size() - 2] *
             qqm.first[qqm.first.size() - 2]) / q0[q0.size() - 2] /
           T(int(q0.size())), qqm.second));
    }
    M[4] /= (Mx[4] = max(Mx[4], abs(M[4])));
    {
      const auto& ff(f.next(d));
      Mx[6] = max(abs(d), Mx[6]);
      if(f.full) {
        M[6] = ff[0];
        for(int i = 1; i < ff.size(); i ++) M[6] += ff[i];
        M[6] = (- M[6]) / (Mx[6] = max(Mx[6], abs(M[6])));
      }
    }
    MM = T(int(0));
    for(int i = 0; i < M.size(); i ++) if(isfinite(M[i])) MM += M[i];
    return MM *= Mx0 / T(int(M.size() - 2));
  }
  P0maxRank<T> p0;
  shrinkMatrix<T, P1I<T, idFeeder<T> > > p1;
  shrinkMatrix<T, P012L<T, idFeeder<T> > > p2;
  idFeeder<T> q;
  SimpleVector<T> q0;
  idFeeder<T> f;
  T MM;
  T Mx0;
  vector<T> M;
  vector<T> Mx;
};

// N.B. subtract minimum square.
template <typename T, typename P> class Pmss {
public:
  inline Pmss() { ; }
  inline Pmss(P&& p, const int& len = 2) {
    this->p = p;
    r = q = idFeeder<T>(len);
  }
  inline ~Pmss() { ; }
  inline T next(const T& d) {
    const auto& msv(q.next(d));
    const auto  ms(msv.dot(mscache<T>(msv.size())));
    if(! q.full) return T(int(0));
    const auto& mssv(r.next(ms));
    return r.full ? p.next(d - ms) + pnextcacher<T>(mssv.size(), 1, 2).dot(mssv) / T(int(2)) : T(int(0));
  }
  P p;
  idFeeder<T> q;
  idFeeder<T> r;
};

// In laurent series, we treat a_-1 as both side a_1 and a_-1 arithmetric avg.
template <typename T, typename P> class PWalkBoth {
public:
  inline PWalkBoth() { ; }
  inline PWalkBoth(P&& p) { S = T(int(0)); this->p = q = p; }
  inline ~PWalkBoth() { ; }
  inline T next(const T& d) {
    static const T zero(int(0));
    static const T one(int(1));
    static const T two(int(2));
    const auto bS(S);
    if((S += d) != zero && bS != zero) {
      const auto dd(one / S - one / bS);
            auto qd(q.next(dd));
      qd = qd == zero ? qd : one / (one / S + one / qd) - S;
      return (p.next(d) + (isfinite(qd) ? qd : zero)) / two;
    }
    return p.next(d);
  }
  T S;
  P p;
  P q;
};

// N.B. [- 1, 1] prediction with certain range.
template <typename T> class Prange {
public:
  inline Prange() { ; }
  inline Prange(const int& status) {
    assert(0 < status);
    const int var0(max(T(int(1)), T(int(exp(sqrt(log(T(status)))))) ) );
    const int var1(max(T(int(2)), pow(T(status), T(int(1)) / T(int(3)))));
    const int var2(max(T(int(2)), pow(T(status), T(int(1)) / T(int(4)))));
    p0 = P0maxRank<T>(status - var0 - 3);
    p1 = shrinkMatrix<T, P1I<T, idFeeder<T> > >(P1I<T, idFeeder<T> >(status - var1 * 2, var1, var1), var1);
    p2 = shrinkMatrix<T, P012L<T, idFeeder<T> > >(P012L<T, idFeeder<T> >(status - var2 * 2, var2, var2), var2);
    const int qstatus(sqrt(T(status)));
    q  = idFeeder<T>(qstatus);
    q0 = SimpleVector<T>(qstatus + 1).O();
  }
  inline ~Prange() { ; }
  inline T next(T d) {
    static const T one(int(1));
    auto M(max(- one, min(one, p0.next(d))) );
    M += max(- one, min(one, p1.next(d)));
    M += max(- one, min(one, p2.next(d)));
    {
      auto qm(makeProgramInvariant<T>(q.next(d)));
      q0 += qm.first * pow(qm.second, ceil(- log(SimpleMatrix<T>().epsilon())));
      auto qq(q);
      auto qqm(makeProgramInvariant<T>(qq.next(d)));
      M += max(- one, min(one, revertProgramInvariant<T>(make_pair(
        - (q0.dot(qqm.first) - q0[q0.size() - 2] *
             qqm.first[qqm.first.size() - 2]) / q0[q0.size() - 2] /
           T(int(q0.size())), qqm.second)) ));
    }
    return max(- T(int(1)), min(T(int(1)), (M += d) *= T(int(2)) / T(int(5)) ));
  }
  P0maxRank<T> p0;
  shrinkMatrix<T, P1I<T, idFeeder<T> > > p1;
  shrinkMatrix<T, P012L<T, idFeeder<T> > > p2;
  idFeeder<T> q;
  SimpleVector<T> q0;
};

// N.B. x -> x^(2/3) replacement. This causes (cbrt(x))^2 series or rich of
//                   information on sparse input matrix.
template <typename T, typename P> class Ppad {
public:
  inline Ppad() { ; }
  inline Ppad(P&& p, const int& pad = 1, const int& absent = - 1) {
    this->p = p0 = p;
    nxt = int(pow(T(this->pad = pad0 = pad), T(int(2)) / T(int(3))));
    assert(0 < pad);
    d = res = T(tt = t ^= t);
    tt ++;
    this->absent = absent < 0 ? pad : absent;
  }
  inline ~Ppad() { ; }
  inline const T& next(const T& in) {
    d += in;
    if(t ++ < nxt) return res;
    const auto nxt2(pow(T(pad += pad0), T(int(2)) / T(int(3))));
    if(nxt2 < T(nxt + 1) + T(int(1)) / T(int(4))) {
      nxt = int(pow(T(pad = pad0), T(int(2)) / T(int(3))));
      p   = p0;
      return res = d = T(tt = t ^= t);
    }
    nxt = int(nxt2);
    res = p.next(d);
    d   = T(int(0));
    return tt ++ < absent ? res = T(int(0)) : res;
  }
  int t;
  int tt;
  int absent;
  int pad0;
  int pad;
  int nxt;
  P p0;
  P p;
  T d;
  T res;
};

// N.B. auto retry Ppad, we need this because x -> large, delta x^(2/3) -> 0.
template <typename T, typename P> class Ppretry {
public:
  inline Ppretry() { ; }
  inline Ppretry(P&& p) { p0 = p; flag = true; }
  inline ~Ppretry() { ; }
  inline T next(const T& in) {
    if(flag) {
      if(! p.size() || p[0].tt)
        p.emplace_back(p0);
      else if(p.size()) flag = ! flag;
    }
    T res(int(0));
    for(int i = 0; i < p.size(); i ++)
      res += p[i].next(in);
    return res;
  }
  P p0;
  vector<P> p;
  bool flag;
};

template <typename T> pair<vector<SimpleVector<T> >, vector<SimpleVector<T> > > predv(const vector<SimpleVector<T> >& in) {
  vector<Prange<T> > p0;
  for(int ext = 0; ext < in.size() / 2; ext ++) {
    const int status(in.size() / (ext + 1) - 2);
    const int var0(max(T(int(1)), T(int(exp(sqrt(log(T(status)))))) ) );
    if(status < var0 + 3 * 2) break;
    p0.emplace_back(Prange<T>(status));
    auto pp(p0[ext]);
    for(int i = 0; i < status * 2 + 4; i ++)
      pp.next(T(i + 1) / T(status * 2 + 5) - T(int(1)) / T(int(2)));
    cerr << "(volatile dummy:)" << pp.next(T(int(0))) << endl;
  }
  vector<SimpleVector<T> > invariant;
  invariant.resize(in.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int i = 0; i < in.size(); i ++) {
    auto inv(makeProgramInvariant<T>(in[i]));
    invariant[i]  = move(inv.first);
    invariant[i] *=
      pow(inv.second, ceil(- log(SimpleMatrix<T>().epsilon()) ));
  }
  vector<SimpleVector<T> > p;
  p.resize(p0.size());
  auto q(p);
  for(int i = 0; i < p0.size(); i ++) {
    cerr << i << " / " << p0.size() << endl;
    p[i].resize(invariant[0].size());
    q[i].resize(invariant[0].size());
    p[i].O();
    q[i].O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int j = 0; j < p[i].size(); j ++) {
      auto pf(p0[i]);
      auto pb(p0[i]);
      T    qmax(int(0));
      T    pmax(int(0));
      for(int k = 0; k < invariant.size() / (i + 1); k ++)
        qmax += invariant[(invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] *
          invariant[(invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j];
      for(int k = 0; k < invariant.size() / (i + 1); k ++)
        pmax += invariant[invariant.size() - 1 -
            (invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] *
                invariant[invariant.size() - 1 -
            (invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j];
      qmax = sqrt(qmax * T(int(2)));
      pmax = sqrt(pmax * T(int(2)));
      try {
        for(int k = 0; k < invariant.size() / (i + 1); k ++) {
          assert(0 <= (invariant.size() / (i + 1) - (k + 1)) * (i + 1));
          assert((invariant.size() / (i + 1) - (k + 1)) * (i + 1) < invariant.size());
          assert(- T(int(1)) <= invariant[(invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] / qmax);
          assert(invariant[(invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] / qmax <= T(int(1)));
          q[i][j] = pb.next(invariant[(invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] / qmax);
        }
      } catch(const char* e) {
        q[i][j] = T(int(0));
      }
      try {
        for(int k = 0; k < invariant.size() / (i + 1); k ++) {
          assert(0 <= (invariant.size() / (i + 1) - (k + 1)) * (i + 1));
          assert((invariant.size() / (i + 1) - (k + 1)) * (i + 1) < invariant.size());
          assert(- T(int(1)) <= invariant[invariant.size() - 1 -
            (invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] / pmax);
          assert(invariant[invariant.size() - 1 -
            (invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] / pmax <= T(int(1)));
          p[i][j] = pf.next(invariant[invariant.size() - 1 -
              (invariant.size() / (i + 1) - (k + 1)) * (i + 1)][j] / pmax);
        }
      } catch(const char* e) {
        p[i][j] = T(int(0));
      }
      q[i][j] *= qmax;
      p[i][j] *= pmax;
    }
  }
  return make_pair(move(p), move(q));
}

// N.B. invariant gathers some of the group on the input pattern.
template <typename T> SimpleMatrix<T> concat(const SimpleMatrix<T>& m0, const SimpleMatrix<T>& m1) {
  // det diag result = det diag m0 + det diag m1
  // [1 x x^reverse 1]
  // if we met rank shrink, assert exit then.
  // we can handle this with compiling operation adding period-depend values.
  assert(m0.rows() == m1.rows() && m0.cols() == m1.cols());
  SimpleMatrix<T> work0(m0);
  SimpleMatrix<T> work1(m1);
  auto res(m0);
  for(int i = 0; i < m0.rows(); i ++) {
    auto qw1(work1.transpose().QR());
    auto rw1(qw1 * work1.transpose());
    // XXX : assert exit here.
    work0 = (rw1.inverse() * qw1 * work0.transpose()).transpose();
    assert(work0.rows() == work0.cols());
    SimpleMatrix<T> lwork(work0.rows() * 2, work0.cols() * 2);
    const auto ii(SimpleMatrix<T>(work0.rows(), work0.cols()).I());
    lwork.setMatrix(0, 0, ii).setMatrix(0, work0.cols(), work0 - ii).setMatrix(work0.rows(), work0.cols(), ii);
    for(int j = 0; j < work0.rows(); j ++)
      for(int k = 0; k < work0.cols(); k ++)
        lwork(lwork.rows() - j, work0.cols() - k) = work0(j, k);
    for(int j = 0; j < lwork.rows() / 2; j ++)
      for(int k = 0; k < lwork.cols(); k ++)
        std::swap(lwork(j, k), lwork(lwork.rows() - j, k));
    for(int j = 0; j < lwork.rows(); j ++)
      for(int k = 0; k < lwork.cols() / 2; k ++)
        std::swap(lwork(j, k), lwork(j, lwork.cols() - k));
    // LDLt:
/*
    auto L();
    auto D();
    auto Linv();
    for(int j = 0; j < L.rows() / 2; j ++)
      for(int k = 0; k < L.cols(); k ++)
        std::swap(L(j, k), L(L.rows() - j, k));
    for(int j = 0; j < Linv.rows(); j ++)
      for(int k = 0; k < Linv.cols() / 2; k ++)
        std::swap(Linv(j, k), Linv(j, Linv.cols() - k));
*/
    // factor apply res, work0, work1:
  }
  return res;
}


template <typename T> SimpleMatrix<T> diff(const SimpleMatrix<T>& m, const int& idx) {
  SimpleMatrix<T> res(m.rows() - 1, m.cols());
  res.O();
  for(int i = 0; i < m.rows(); i ++) {
    auto lres(res);
    lres.O();
    for(int j = 0; j < i; j ++) lres.row(j) = m.row(j);
    for(int j = i + 1; j < m.rows(); j ++) lres.row(j - 1) = m.row(j);
    if(m(i, i) == T(int(0))) continue;
    else if(m(i, i) < T(int(0))) lres.row(0) = - lres.row(0);
    res = concat(res, lres /= pow(abs(m(i, i)), T(int(1)) / T(int(lres.rows()))));
  }
}

template <typename T> SimpleMatrix<T> integrate(const SimpleMatrix<T>& m, const int& idx, const int& stage = 0) {
  // N.B. S^x det diag Ax dx (= S u'v) =
  //  (S^x dx) * det diag Ax (= S(uv)') -
  //  S^x(x(det diag Ax)')dx (= S uv')
  //      S^x(x(det diag Ax)')dx   (= S u'v) =
  //  (S^x x dx) * (det diag Ax)'  (= S(uv)') -
  //  S^x(x^2/2 (det diag Ax)'')dx (= S uv')
  //    ...
  SimpleMatrix<T> factorial(m.rows(), m.cols());
  factorial.O();
  for(int i = 0; i < factorial.rows(); i ++)
    factorial(i, idx) = T(int(1)) / T(int(i + 1));
  if(stage == m.rows() - 1) return factorial * m(m.rows() - 1, idx);
  return concat(m, factorial.setMatrix(stage + 1, 0, integrate(diff(m, idx), idx, stage + 1)), true);
}

// N.B. we need huge computing power depends on m at least O(m.rows()^4).
template <typename T> SimpleVector<T> reduce(const SimpleMatrix<T> m) {
  SimpleMatrix<T> work(m);
  for(int i = 0; i < m.rows() - 1; i ++)
    work = integrate(work, i);
  for(int i = 0; i < m.rows() - 1; i ++)
    work = diff(work, i);
  return work.row(0);
}

#define _P2_
#endif

