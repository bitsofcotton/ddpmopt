#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#include <random>
#include <assert.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

//#define int int64_t
#define int int32_t
#include "lieonn.hh"
typedef myfloat num_t;
#include "p0.hh"
#include "p1.hh"
#include "catg.hh"
#include "decompose.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::to_string;
using std::vector;
using std::sort;
using std::binary_search;
using std::make_pair;
using std::istringstream;

#include <stdlib.h>

template <typename T> class P {
public:
  inline P() { ; }
  inline P(const int& status) {
    assert(0 < status);
    const int var0(max(T(int(1)), T(int(exp(sqrt(log(T(status)))))) ) );
    const int var1(max(T(int(2)), pow(T(status), T(int(1)) / T(int(3)))));
    const int var2(max(T(int(2)), pow(T(status), T(int(1)) / T(int(4)))));
    p0 = P0maxRank<T>(status - var0 - 3);
    p1 = shrinkMatrix<T, P1I<T, idFeeder<T> > >(P1I<T, idFeeder<T> >(status - var1 * 2, var1, var1), var1);
    p2 = shrinkMatrix<T, P012L<T, idFeeder<T> > >(P012L<T, idFeeder<T> >(status - var2 * 2, var2, var2), var2);
    const int qstatus(sqrt(num_t(status)));
    q  = idFeeder<T>(qstatus);
    q0 = SimpleVector<T>(qstatus + 1).O();
  }
  inline ~P() { ; }
  inline T next(T d) {
    static const T one(int(1));
    auto M(max(- one, min(one, p0.next(d))) );
    M += max(- one, min(one, p1.next(d)));
    M += max(- one, min(one, p2.next(d)));
    {
      auto qm(makeProgramInvariant<T>(q.next(d)));
      q0 += std::move(qm.first) * pow(qm.second, ceil(- log(SimpleMatrix<T>().epsilon())));
      auto qq(q);
      auto qqm(makeProgramInvariant<T>(qq.next(d)));
      M += max(- one, min(one, revertProgramInvariant<T>(make_pair(
        - (q0.dot(qqm.first) - q0[q0.size() - 2] *
             qqm.first[qqm.first.size() - 2]) / q0[q0.size() - 2] /
           T(int(q0.size())), qqm.second)) ));
    }
    return max(- num_t(int(1)), min(num_t(int(1)), (M += d) *= num_t(int(2)) / num_t(int(5)) ));
  }
  P0maxRank<T> p0;
  shrinkMatrix<T, P1I<T, idFeeder<T> > > p1;
  shrinkMatrix<T, P012L<T, idFeeder<T> > > p2;
  idFeeder<T> q;
  SimpleVector<T> q0;
};

static inline bool whiteline(const string& s) {
  for(auto ss(s.begin()); ss < s.end(); ++ ss)
    if(! std::isspace(* ss) && *ss != '\n')
      return false;
  return true;
}

template <typename T> static inline bool loadstub(ifstream& input, const int& nmax, const int& ncolor, vector<SimpleMatrix<T> >& datas) {
  int i = 0, j = 0, k = 0;
  char buf;
  int  work = 0;
  bool mode = false;
  while(input.get(buf) && j < datas[0].rows()) {
    if('0' <= buf && buf <= '9') {
      work *= 10;
      work += buf - '0';
      mode  = true;
      continue;
    } else if(mode) {
      mode = false;
      datas[k](j, i) = T(work) / T(nmax);
      work = 0;
      if(++ k >= ncolor) {
        if(++ i >= datas[0].cols()) {
          if(++ j >= datas[0].rows())
            return true;
          i = 0;
        }
        k = 0;
      }
    }
  }
  return true;
}

template <typename T> bool loadp2or3(vector<SimpleMatrix<T> >& data, const char* filename) {
  string line;
  string line2;
  string line3;
  ifstream input;
  input.open(filename);
  if(input.is_open()) {
    try {
      data.resize(3, SimpleMatrix<T>());
      getline(input, line);
      while((whiteline(line) || line[0] == '#') && getline(input, line) && !input.eof() && !input.bad()) ;
      getline(input, line2);
      while((whiteline(line2) || line2[0] == '#') && getline(input, line2) && !input.eof() && !input.bad()) ;
      getline(input, line3);
      while((whiteline(line3) || line3[0] == '#') && getline(input, line3) && !input.eof() && !input.bad()) ;
      istringstream iline2(line2);
      int w, h;
      iline2 >> w;
      iline2 >> h;
      if(line.size() < 2 || w <= 0 || h <= 0) {
        cerr << "unknown size." << endl;
        input.close();
        return false;
      } 
      istringstream iline3(line3);
      int nmax;
      iline3 >> nmax;
      if(line[0] == 'P') {
        if(line[1] == '2') {
          data[0] = SimpleMatrix<T>(h, w).O();
          loadstub<T>(input, nmax, 1, data);
          data[1] = data[2] = data[0];
        } else if(line[1] == '3') {
          for(int i = 0; i < 3; i ++)
            data[i] = SimpleMatrix<T>(h, w).O();
          loadstub<T>(input, nmax, 3, data);
        } else {
          cerr << "unknown file type." << endl;
          input.close();
          return false;
        }
      } else {
        cerr << "unknown file type." << endl;
        input.close();
        return false;
      }
    } catch (...) {
      cerr << "Exception while reading." << endl;
    }
    input.close();
  } else {
    cerr << "Unable to open file for read: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> bool savep2or3(const char* filename, const vector<SimpleMatrix<T> >& data, const bool& gray, const int& depth = 255) { ofstream output;
  output.open(filename);
  if(output.is_open()) {
    try {
      if(gray)
        output << "P2" << "\n";
      else
        output << "P3" << "\n";
      output << data[0].cols() << " " << data[0].rows() << "\n";
      output << depth << "\n";
      for(int i = 0; i < data[0].rows(); i ++)
        for(int j = 0; j < data[0].cols(); j ++)
          if(gray)
            output << int(data[0](i, j) * T(depth)) << "\n";
          else
            for(int k = 0; k < 3; k ++)
              output << int(data[k](i, j) * T(depth)) << "\n";
    } catch (...) {
      cerr << "An error has occured while writing file." << endl;
    }
    output.close();
    output.close();
  } else {
    cerr << "Unable to open file for write: " << filename << endl;
    return false;
  }
  return true;
}

template <typename T> vector<SimpleMatrix<T> > normalize(const vector<SimpleMatrix<T> >& data, const T& upper = T(1)) {
  T MM(0), mm(0);
  bool fixed(false);
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        if(! fixed || (isfinite(data[k](i, j)) && ! isinf(data[k](i, j)) && ! isnan(data[k](i, j)))) {
          if(! fixed)
            MM = mm = data[k](i, j);
          else {
            MM = max(MM, data[k](i, j));
            mm = min(mm, data[k](i, j));
          }
          fixed = true;
        }
  if(MM == mm || ! fixed)
    return data;
  auto result(data);
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++) {
        if(isfinite(result[k](i, j)) && ! isinf(data[k](i, j)) && ! isnan(result[k](i, j)))
          result[k](i, j) -= mm;
        else
          result[k](i, j)  = T(0);
        assert(T(0) <= result[k](i, j) && result[k](i, j) <= MM - mm);
        result[k](i, j) *= upper / (MM - mm);
      }
  return result;
}

static inline num_t rng() {
  myuint res(0);
  // XXX: we don't trust system or compiler PRNG.
  // static std::random_device rd;
  // XXX: we want natural, deterministic, better PRNG, however,
  //      we don't search deepinside of this PRNG.
  //      (might be predecessor exists.)
  static uint64_t t(1);
  assert(t && "rng() should not be periodical.");
#if defined(_FLOAT_BITS_)
  for(int i = 0; i < _FLOAT_BITS_ / sizeof(uint32_t) / 8; i ++) {
#else
  for(int i = 0; i < 2; i ++) {
#endif
    res <<= sizeof(uint32_t) * 8;
#if defined(_FLOAT_BITS_)
typedef SimpleFloat<uint64_t, unsigned __int128, 64, int64_t> thisfl;
#else
typedef long double thisfl;
#endif
    auto buf(sin(thisfl(t ++)) * pow(thisfl(int(2)), thisfl(int(32))));
    buf  -= floor(buf);
    res  |= uint32_t(int(buf * pow(thisfl(int(2)), thisfl(int(32)) )));
#undef thisfl
    // res  |= uint32_t(rd());
  }
  return max(num_t(int(0)), min(num_t(int(1)), num_t(res) / num_t(~ myuint(0)) ));
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

template <typename T> SimpleMatrix<T> shrinken(const SimpleMatrix<T>& in, const int& h = 2, const int& w = 2) {
  SimpleMatrix<T> shrink(h, w);
  shrink.O();
  for(int ii = 0; ii < h; ii ++)
    for(int jj = 0; jj < w; jj ++) {
      int cnt(0);
      for(int iik = 0;
          iik < min(in.rows() / h,
            in.rows() - ii * (in.rows() / h)); iik ++)
        for(int jjk = 0; jjk < min(in.cols() / w,
              in.cols() - jj * (in.cols() / w));
            jjk ++, cnt ++)
          shrink(ii, jj) +=
            in(ii * (in.rows() / h) + iik,
                     jj * (in.cols() / w) + jjk);
      shrink(ii, jj) /= num_t(cnt);
    }
  return shrink;
}

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  const auto ext(std::atoi(argv[1]));
  for(int i0 = 2; i0 < argc; i0 ++) {
    vector<SimpleMatrix<num_t> > in;
    if(! loadp2or3<num_t>(in, argv[i0])) continue;
    const num_t diag(in[0].rows() * in[0].cols());
    const int   sz(sqrt(diag));
    const int   num(diag * log(diag));
    vector<SimpleVector<num_t> > noisep;
    vector<SimpleVector<num_t> > noiseq;
    noisep.resize(num, SimpleVector<num_t>(sz));
    noiseq.resize(num, SimpleVector<num_t>(sz));
    for(int j = 0; j < noisep.size(); j ++) {
      for(int n = 0; n < noisep[j].size(); n ++)
        noisep[j][n] = rng();
      noisep[j] = (dft<num_t>(- noisep[j].size()) * noisep[j].template cast<complex<num_t> >()).template real<num_t>();
    }
    for(int j = 0; j < noiseq.size(); j ++) {
      for(int n = 0; n < noiseq[j].size(); n ++)
        noiseq[j][n] = rng();
      noiseq[j] = (dft<num_t>(- noiseq[j].size()) * noiseq[j].template cast<complex<num_t> >()).template real<num_t>();
    }
    auto out(in);
    for(int i = 0; i < out.size(); i ++) {
      out[i].resize(in[i].rows() + ext * 2, in[i].cols() + ext * 2);
      out[i].O().setMatrix(ext, ext, in[i]);
    }
    SimpleMatrix<num_t> riny(ext, out[0].cols());
    SimpleMatrix<num_t> rinx(out[0].rows(), ext);
    for(int n = 0; n < riny.rows(); n ++)
      for(int nn = 0; nn < riny.cols(); nn ++)
        riny(n, nn) = rng();
    for(int n = 0; n < rinx.rows(); n ++)
      for(int nn = 0; nn < rinx.cols(); nn ++)
        rinx(n, nn) = rng();
    for(int i = 0; i < ext; i ++) {
      P<num_t> p0(in[0].rows() + i * 2 - 2);
      vector<SimpleVector<num_t> > yp;
      yp.resize(out.size(), SimpleVector<num_t>(in[0].cols() + i * 2).O());
      auto ym(yp);
      cerr << "Step 1: " << i << " / " << ext << " over " << i0 << " / " << argc - 1 << std::endl;
      for(int k = 0; k < out.size(); k ++) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int ii = i; ii < out[k].cols() - i; ii ++) {
          auto pf(p0);
          auto pb(p0);
          try {
            for(int jj = i; jj < out[k].rows() - i; jj ++)
              ym[k][ii - i] = pf.next(out[k](out[k].rows() - i - 1 - jj, ii));
          } catch(const char* e) {
            ym[k][ii - i] = num_t(int(0));
          }
          try {
            for(int jj = i; jj < out[k].rows() - i; jj ++)
              yp[k][ii - i] = pb.next(out[k](jj, ii));
          } catch(const char* e) {
            yp[k][ii - i] = num_t(int(0));
          }
        }
      }
      P<num_t> q0(in[0].cols() + i * 2 - 2);
      vector<SimpleVector<num_t> > xp;
      xp.resize(out.size(), SimpleVector<num_t>(in[0].rows() + i * 2).O());
      auto xm(xp);
      cerr << "Step 2: " << i << " / " << ext << " over " << i0 << " / " << argc - 1 << std::endl;
      for(int k = 0; k < out.size(); k ++) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int ii = i; ii < out[k].rows() - i; ii ++) {
          auto pf(q0);
          auto pb(q0);
          try {
            for(int jj = i; jj < out[k].cols() - i; jj ++)
              ym[k][ii - i] = pf.next(out[k](ii, out[k].cols() - i - 1 - jj));
          } catch(const char* e) {
            ym[k][ii - i] = num_t(int(0));
          }
          try {
            for(int jj = i; jj < out[k].cols() - i; jj ++)
              yp[k][ii - i] = pb.next(out[k](ii, jj));
          } catch(const char* e) {
            yp[k][ii - i] = num_t(int(0));
          }
        }
      }
      cerr << "Step 3: " << i << " / " << ext << " over " << i0 << " / " << argc - 1 << std::endl;
      vector<SimpleMatrix<num_t> > Lp;
      Lp.resize(out.size());
      for(int ii = 0; ii < Lp.size(); ii ++) {
        auto shrinky(shrinken<num_t>(in[ii].subMatrix(i, i, in[ii].rows() - i * 2, in[ii].cols() - i * 2), in[ii].rows() - i * 2, sz));
        Lp[ii].resize(in[0].cols(), sz + 2);
        Lp[ii].O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int m = 0; m < in[0].cols(); m ++) {
          SimpleMatrix<num_t> work(num, sz + 2);
          for(int jj = 0; jj < num; jj ++) {
            SimpleVector<num_t> vwork(shrinky.cols() + 1);
            for(int n = 0; n < shrinky.cols(); n ++)
              vwork[n] = shrinky(jj * shrinky.rows() / num, n) * noisep[jj][n];
            vwork[vwork.size() - 1] = in[ii](jj * shrinky.rows() / num, m);
            auto mpi(makeProgramInvariant<num_t>(vwork));
            work.row(jj)  = move(mpi.first);
            work.row(jj) *=
              pow(mpi.second, ceil(- log(in[0].epsilon()) ));
          }
          Lp[ii].row(m) = linearInvariant(work);
        }
      }
      vector<SimpleMatrix<num_t> > Lq;
      Lq.resize(0);
      Lq.resize(out.size());
      for(int ii = 0; ii < Lq.size(); ii ++) {
        auto shrinkx(shrinken<num_t>(in[ii].subMatrix(i, i, in[ii].rows() - i * 2, in[ii].cols() - i * 2), sz, in[ii].cols() - i * 2));
        Lq[ii].resize(in[0].rows(), sz + 2);
        Lq[ii].O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int m = 0; m < in[0].rows(); m ++) {
          SimpleMatrix<num_t> work(num, sz + 2);
          for(int jj = 0; jj < num; jj ++) {
            SimpleVector<num_t> vwork(shrinkx.rows() + 1);
            for(int n = 0; n < shrinkx.rows(); n ++)
              vwork[n] = shrinkx(n, jj * shrinkx.cols() / num) * noiseq[jj][n];
            vwork[vwork.size() - 1] = in[ii](m, jj * shrinkx.cols() / num);
            auto mpi(makeProgramInvariant<num_t>(vwork));
            work.row(jj)  = move(mpi.first);
            work.row(jj) *=
              pow(mpi.second, ceil(- log(in[0].epsilon()) ));
          }
          Lq[ii].row(m) = linearInvariant(work);
        }
      }
      for(int k = 0; k < Lp.size(); k ++) {
        auto pL(Lp[k]);
        auto qL(Lq[k]);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int ii = 0; ii < Lp[k].rows(); ii ++) {
          for(int jj = 0; jj < Lp[k].cols(); jj ++) {
            auto pf(p0);
            auto pb(p0);
            num_t qLn(int(0));
            num_t pLn(int(0));
            for(int kk = 0; kk < Lp.size() / (i + 1); kk ++)
              qLn += Lp[(Lp.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) *
                     Lp[(Lp.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj);
            for(int kk = 0; kk < Lp.size() / (i + 1); kk ++)
              pLn += Lp[Lp.size() - 1 -
                (Lp.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) *
                     Lp[Lp.size() - 1 -
                (Lp.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj);
            qLn = sqrt(qLn) * num_t(int(2));
            pLn = sqrt(pLn) * num_t(int(2));
            try {
              for(int kk = 0; kk < Lp.size() / (i + 1); kk ++)
                qL(ii, jj) = pb.next(
                  Lp[(Lp.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) / qLn);
            } catch(const char* e) {
              qL(ii, jj) = num_t(int(0));
            }
            try {
              for(int kk = 0; kk < Lp.size() / (i + 1); kk ++)
                pL(ii, jj) = pf.next(
                  Lp[Lp.size() - 1 -
                  (Lp.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) / pLn);
            } catch(const char* e) {
              pL(ii, jj) = num_t(int(0));
            }
            qL(ii, jj) *= qLn;
            pL(ii, jj) *= pLn;
          }
          pL.row(ii) /= num_t(pL(ii, pL.cols() - 2));
          pL(ii, pL.cols() - 2) = num_t(int(0));
          qL.row(ii) /= num_t(qL(ii, qL.cols() - 2));
          qL(ii, qL.cols() - 2) = num_t(int(0));
        }
        out[k].row(i).setVector(ext - i - 1, qL * makeProgramInvariant<num_t>(
          SimpleVector<num_t>(ym[k].size() + 1).O().setVector(0, ym[k])).first);
        out[k].row(out[k].rows() - i - 1).setVector(ext - i - 1,
          pL * makeProgramInvariant<num_t>(
            SimpleVector<num_t>(yp[k].size() + 1).O().setVector(0, yp[k])).first);
        pL = Lp[k];
        qL = Lq[k];
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int ii = 0; ii < Lq[k].rows(); ii ++) {
          for(int jj = 0; jj < Lq[k].cols(); jj ++) {
            auto pf(q0);
            auto pb(q0);
            num_t qLn(int(0));
            num_t pLn(int(0));
            for(int kk = 0; kk < Lq.size() / (i + 1); kk ++)
              qLn += Lq[(Lq.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) *
                     Lq[(Lq.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj);
            for(int kk = 0; kk < Lq.size() / (i + 1); kk ++)
              pLn += Lq[Lq.size() - 1 -
                (Lq.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) *
                     Lq[Lq.size() - 1 -
                (Lq.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj);
            qLn = sqrt(qLn) * num_t(int(2));
            pLn = sqrt(pLn) * num_t(int(2));
            try {
              for(int kk = 0; kk < Lq.size() / (i + 1); kk ++)
                qL(ii, jj) = pb.next(
                  Lq[(Lq.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) / qLn);
            } catch(const char* e) {
              qL(ii, jj) = num_t(int(0));
            }
            try {
              for(int kk = 0; kk < Lq.size() / (i + 1); kk ++)
                pL(ii, jj) = pf.next(
                  Lq[Lq.size() - 1 -
                  (Lq.size() / (i + 1) - (kk + 1)) * (i + 1)](ii, jj) / pLn);
            } catch(const char* e) {
              pL(ii, jj) = num_t(int(0));
            }
            qL(ii, jj) *= qLn;
            pL(ii, jj) *= pLn;
          }
          pL.row(ii) /= num_t(pL(ii, pL.cols() - 2));
          pL(ii, pL.cols() - 2) = num_t(int(0));
          qL.row(ii) /= num_t(qL(ii, qL.cols() - 2));
          qL(ii, qL.cols() - 2) = num_t(int(0));
        }
        out[k].setCol(i, SimpleVector<num_t>(out[k].col(i)).setVector(
          ext - i - 1,
          qL * makeProgramInvariant<num_t>(
            SimpleVector<num_t>(xm[k].size() + 1).O().setVector(0, xm[k])).first));
        out[k].setCol(out[k].cols() - i - 1, SimpleVector<num_t>(out[k].col(
          out[k].cols() - i - 1)).setVector(ext - i - 1,
          pL * makeProgramInvariant<num_t>(
            SimpleVector<num_t>(xp[k].size() + 1).O().setVector(0, xp[k])).first));
      }
    }
    if(! savep2or3<num_t>(argv[i0], normalize<num_t>(out), false, 65535) )
      std::cerr << "failed to save." << std::endl;
  }
  return 0;
}

