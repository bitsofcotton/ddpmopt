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
      q0 += qm.first * pow(qm.second, ceil(- log(SimpleMatrix<T>().epsilon())));
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

template <typename T> bool savep2or3(const char* filename, const vector<SimpleMatrix<T> >& data, const bool& gray, const int& depth = 255) {
  ofstream output;
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

template <typename T> vector<vector<SimpleMatrix<T> > > shrinken(const vector<vector<SimpleMatrix<T> > >& in, const int& sz = 2) {
  vector<vector<SimpleMatrix<T> > > shrink;
  shrink.resize(in.size());
  for(int i = 0; i < shrink.size(); i ++) {
    shrink[i].resize(in[i].size(), SimpleMatrix<T>(sz, sz).O());
    for(int j = 0; j < shrink[i].size(); j ++)
      for(int ii = 0; ii < sz; ii ++)
        for(int jj = 0; jj < sz; jj ++) {
          int cnt(0);
          for(int iik = 0;
              iik < min(in[i][j].rows() / sz,
                in[i][j].rows() - ii * (in[i][j].rows() / sz)); iik ++)
            for(int jjk = 0; jjk < min(in[i][j].cols() / sz,
                  in[i][j].cols() - jj * (in[i][j].cols() / sz));
                jjk ++, cnt ++)
              shrink[i][j](ii, jj) +=
                in[i][j](ii * (in[i][j].rows() / sz) + iik,
                         jj * (in[i][j].cols() / sz) + jjk);
          shrink[i][j](ii, jj) /= num_t(cnt);
        }
  }
  return shrink;
}

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  vector<vector<SimpleMatrix<num_t> > > in;
  vector<vector<SimpleMatrix<num_t> > > noise;
  vector<vector<SimpleMatrix<num_t> > > p;
  vector<vector<SimpleMatrix<num_t> > > L;
  vector<vector<SimpleMatrix<num_t> > > pL;
  vector<P<num_t> > p0;
  int sz(0);
  int num(0);
  in.resize(argc - 1);
  noise.resize(in.size());
  for(int ext = 0; ext < in.size() / 2; ext ++) {
    const int status(in.size() / (ext + 1) - 2);
    const int var0(max(num_t(int(1)), num_t(int(exp(sqrt(log(num_t(status)))))) ) );
    if(status < var0 + 3 * 2) break;
    p0.emplace_back(P<num_t>(status));
    auto pp(p0[ext]);
    for(int i = 0; i < status * 2 + 4; i ++)
      pp.next(num_t(i + 1) / num_t(status * 2 + 5) - num_t(int(1)) / num_t(int(2)));
    cerr << "(volatile dummy:)" << pp.next(num_t(int(0))) << std::endl;
  }
  p.resize(p0.size());
  auto q(p);
  for(int i = 1; i < argc; i ++) {
    if(! loadp2or3<num_t>(in[i - 1], argv[i])) continue;
    assert(in[0][0].rows() == in[i - 1][0].rows() &&
           in[0][0].cols() == in[i - 1][0].cols());
    if(i == 1) {
      sz  = int(sqrt(num_t(min(int(in.size()), int(in[i - 1][0].rows())))));
      // XXX:
      num = int(num_t(in.size()) * log(num_t(in.size())));
    }
    noise[i - 1].resize(num, SimpleMatrix<num_t>(sz, sz));
    for(int j = 0; j < num; j ++) {
      for(int n = 0; n < sz; n ++)
        for(int nn = 0; nn < sz; nn ++)
          noise[i - 1][j](n, nn) = rng();
      noise[i - 1][j] = (dft<num_t>(- noise[i - 1][j].rows()) * noise[i - 1][j].template cast<complex<num_t> >() * dft<num_t>(- noise[i - 1][j].cols())).template real<num_t>();
    }
  }
  for(int i = 0; i < p.size(); i ++) {
    p[i].resize(in[i].size(), SimpleMatrix<num_t>(in[i][0].rows(), in[i][0].cols()).O());
    q[i].resize(in[i].size(), SimpleMatrix<num_t>(in[i][0].rows(), in[i][0].cols()).O());
    for(int k = 0; k < p[i].size(); k ++) {
      cerr << "Step 1: " << i * p[i].size() + k << " / " << p.size() * p[i].size() << std::endl;
      for(int ii = 0; ii < p[i][k].rows(); ii ++) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int jj = 0; jj < p[i][k].cols(); jj ++) {
          auto pf(p0[i]);
          auto pb(p0[i]);
          try {
            for(int kk = 0; kk < in.size() / (i + 1); kk ++) {
              assert(0 <= (in.size() / (i + 1) - (kk + 1)) * (i + 1));
              assert((in.size() / (i + 1) - (kk + 1)) * (i + 1) < in.size());
              assert(- num_t(int(1)) <= in[(in.size() / (i + 1) - (kk + 1)) * (i + 1)][k](ii, jj));
              assert(in[(in.size() / (i + 1) - (kk + 1)) * (i + 1)][k](ii, jj) <= num_t(int(1)));
              q[i][k](ii, jj) = pb.next(in[(in.size() / (i + 1) - (kk + 1)) * (i + 1)][k](ii, jj));
            }
          } catch(const char* e) {
            q[i][k](ii, jj) = num_t(int(0));
          }
          try {
            for(int kk = 0; kk < in.size() / (i + 1); kk ++) {
              assert(0 <= (in.size() / (i + 1) - (kk + 1)) * (i + 1));
              assert((in.size() / (i + 1) - (kk + 1)) * (i + 1) < in.size());
              assert(- num_t(int(1)) <= in[in.size() - 1 -
                (in.size() / (i + 1) - (kk + 1)) * (i + 1)][k](ii, jj));
              assert(in[in.size() - 1 -
                (in.size() / (i + 1) - (kk + 1)) * (i + 1)][k](ii, jj) <= num_t(int(1)));
              p[i][k](ii, jj) = pf.next(in[in.size() - 1 -
                (in.size() / (i + 1) - (kk + 1)) * (i + 1)][k](ii, jj));
            }
          } catch(const char* e) {
            p[i][k](ii, jj) = num_t(int(0));
          }
        }
      }
    }
  }
  auto shrink(shrinken<num_t>(in, sz));
  L.resize(in.size());
  for(int i = 0; i < in.size(); i ++) {
    L[i].resize(in[0].size());
    for(int j = 0; j < in[0].size(); j ++) {
      L[i][j].resize(in[0][0].rows() * in[0][0].cols(), shrink[0][0].rows() * shrink[0][0].cols() + 2);
      cerr << "Step 2: " << i * in[0].size() + j << " / " << in.size() * in[0].size() << std::endl;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int m = 0; m < in[0][0].rows() * in[0][0].cols(); m ++) {
        SimpleMatrix<num_t> work(num, shrink[0][0].rows() * shrink[0][0].cols() + 2);
        for(int jj = 0; jj < num; jj ++) {
          SimpleVector<num_t> vwork(shrink[i][j].rows() * shrink[i][j].cols() + 1);
          for(int n = 0; n < shrink[i][j].rows(); n ++)
            for(int nn = 0; nn < shrink[i][j].cols(); nn ++)
              vwork[n * shrink[i][j].cols() + nn] = shrink[i][j](n, nn) * noise[i][jj](n, nn);
          vwork[vwork.size() - 1] = in[i][j](m / in[0][0].cols(), m % in[0][0].cols());
          auto mpi(makeProgramInvariant<num_t>(vwork));
          work.row(jj)  = move(mpi.first);
          work.row(jj) *=
            pow(mpi.second, ceil(- log(in[0][0].epsilon()) ));
        }
        L[i][j].row(m) = linearInvariant(work);
      }
    }
  }
  pL.resize(p.size());
  auto qL(pL);
  for(int i = 0; i < p.size(); i ++) {
    pL[i].resize(p[i].size(), SimpleMatrix<num_t>(L[0][0].rows(), L[0][0].cols()).O());
    qL[i].resize(p[i].size(), SimpleMatrix<num_t>(L[0][0].rows(), L[0][0].cols()).O());
    for(int j = 0; j < pL[i].size(); j ++) {
      cerr << "Step 3: " << i * pL[i].size() + j << " / " << pL.size() * pL[i].size() << std::endl;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int ii = 0; ii < pL[i][j].rows(); ii ++) {
        for(int jj = 0; jj < pL[i][j].cols(); jj ++) {
          auto pf(p0[i]);
          auto pb(p0[i]);
          num_t qLn(int(0));
          num_t pLn(int(0));
          for(int kk = 0; kk < L.size() / (i + 1); kk ++)
            qLn += L[(L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) *
                   L[(L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj);
          for(int kk = 0; kk < L.size() / (i + 1); kk ++)
            pLn += L[L.size() - 1 -
              (L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) *
                   L[L.size() - 1 -
              (L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj);
          // N.B. multiply by 2 for division accuracy.
          //      we need this because some of the implementation int32_t
          //      have 64bit integer type causes SimpleFloat glitch.
          qLn = sqrt(qLn) * num_t(int(2));
          pLn = sqrt(pLn) * num_t(int(2));
          try {
            for(int kk = 0; kk < L.size() / (i + 1); kk ++) {
              assert(0 <= (L.size() / (i + 1) - (kk + 1)) * (i + 1));
              assert((L.size() / (i + 1) - (kk + 1)) * (i + 1) < L.size());
              assert(L[(L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / qLn <= num_t(int(1)));
              assert(- num_t(int(1)) <= L[(L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / qLn);
              qL[i][j](ii, jj) = pb.next(
                L[(L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / qLn);
            }
          } catch(const char* e) {
            qL[i][j](ii, jj) = num_t(int(0));
            for(int kk = 0; kk < L.size() / (i + 1); kk ++)
              qL[i][j](ii, jj) += L[(L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / qLn;
            qL[i][j](ii, jj) /= num_t(L.size() / (i + 1));
          }
          try {
            for(int kk = 0; kk < L.size() / (i + 1); kk ++) {
              assert(0 <= (L.size() / (i + 1) - (kk + 1)) * (i + 1));
              assert((L.size() / (i + 1) - (kk + 1)) * (i + 1) < L.size());
              assert(L[L.size() - 1 - 
                (L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / pLn <=
                num_t(int(1)));
              assert(- num_t(int(1)) <= L[L.size() - 1 - 
                (L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / pLn);
              pL[i][j](ii, jj) = pf.next(L[L.size() - 1 - 
                (L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / pLn);
            }
          } catch(const char* e) {
            pL[i][j](ii, jj) = num_t(int(0));
            for(int kk = 0; kk < L.size() / (i + 1); kk ++)
              pL[i][j](ii, jj) += L[L.size() - 1 - 
                (L.size() / (i + 1) - (kk + 1)) * (i + 1)][j](ii, jj) / pLn;
            pL[i][j](ii, jj) /= num_t(L.size() / (i + 1));
          }
          qL[i][j](ii, jj) *= qLn;
          pL[i][j](ii, jj) *= pLn;
        }
        pL[i][j].row(ii) /= num_t(pL[i][j](ii, pL[i][j].cols() - 2));
        pL[i][j](ii, pL[i][j].cols() - 2) = num_t(int(0));
        qL[i][j].row(ii) /= num_t(qL[i][j](ii, qL[i][j].cols() - 2));
        qL[i][j](ii, qL[i][j].cols() - 2) = num_t(int(0));
      }
    }
  }
  cerr << "Step 4 " << std::flush;
  p.insert(p.end(), q.begin(), q.end());
  pL.insert(pL.end(), qL.begin(), qL.end());
  auto sp(shrinken<num_t>(p, sz));
  for(int i = 0; i < p.size(); i ++) {
    auto& out(sp[i]);
    auto& outs(p[i]);
    auto  rin(out[0]);
    for(int n = 0; n < rin.rows(); n ++)
      for(int nn = 0; nn < rin.cols(); nn ++)
        rin(n, nn) = rng();
    rin = (dft<num_t>(- rin.rows()) * rin.template cast<complex<num_t> >() * dft<num_t>(- rin.cols())).template real<num_t>();
    for(int j = 0; j < out.size(); j ++) {
      SimpleVector<num_t> vwork0(out[j].rows() * out[j].cols() + 1);
      for(int n = 0; n < out[j].rows(); n ++)
        for(int nn = 0; nn < out[j].cols(); nn ++)
          vwork0[n * out[j].cols() + nn] = out[j](n, nn) * rin(n, nn);
      vwork0[vwork0.size() - 1] = num_t(int(0));
      auto outwork(pL[i][j] * makeProgramInvariant<num_t>(vwork0).first);
      for(int n = 0; n < outwork.size(); n ++)
        outs[j](n / outs[j].cols(), n % outs[j].cols()) =
          revertProgramInvariant<num_t>(make_pair(outwork[n], num_t(int(1)) ));
    }
    if(! savep2or3<num_t>((std::string("ddpmoptp-") + std::to_string(i) + std::string(".ppm")).c_str(), normalize<num_t>(outs), false, 65535) )
      std::cerr << "failed to save." << std::endl;
  }
  cerr << " Done" << std::endl;
  return 0;
}

