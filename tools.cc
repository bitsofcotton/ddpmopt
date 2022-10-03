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

template <typename T> vector<SimpleMatrix<T> > autoLevel(const vector<SimpleMatrix<T> >& data, const int& count = 0) {
  vector<T> res;
  res.reserve(data[0].rows() * data[0].cols() * data.size());
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        res.emplace_back(data[k](i, j));
  sort(res.begin(), res.end());
  auto result(data);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for(int k = 0; k < data.size(); k ++)
    for(int i = 0; i < data[k].rows(); i ++)
      for(int j = 0; j < data[k].cols(); j ++)
        result[k](i, j) = max(min(data[k](i, j), res[res.size() - count - 1]), res[count]);
  return result;
}

num_t edge(0);
static inline num_t rng() {
  myuint res(0);
#if defined(_FLOAT_BITS_)
  for(int i = 0; i < _FLOAT_BITS_ / sizeof(uint32_t) / 8; i ++) {
#else
  for(int i = 0; i < 2; i ++) {
#endif
    res <<= sizeof(uint32_t) * 8;
    res  |= uint32_t(arc4random());
  }
  return num_t(res) / num_t(~ myuint(0)) * edge;
}

// XXX: Invariant integration meets gulf on finding regularity class and
//      applying them into data class.
//      This is avoidable if we can calculate each of 2 possible operands as
//      the same operation. Otherwise, from somehow, it is learning itself on
//      possible existance causes re-calculate has glitches to the first
//      calculation, 3 or more possible calculation methods divided in some
//      of the operand groups needs non-flat data integrity itself.
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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  vector<SimpleMatrix<num_t> > L;
  const auto step(std::atoi(argv[1]));
  const auto size0(std::atoi(argv[2]));
  const auto recur0(std::atoi(argv[3]));
  const auto size(abs(size0));
  const auto recur(abs(recur0));
  vector<SimpleMatrix<num_t> > out;
  if(! loadp2or3<num_t>(out, argv[4])) return - 1;
  std::random_device rd;
  std::mt19937_64 rde(rd());
  edge = num_t(8) / num_t(4 * int(abs(step) + 1));
  if(step < 0) {
    L.reserve(- step);
    for(int i = 0; i < - step; i ++) {
      SimpleMatrix<num_t> wL;
      std::cin >> wL;
      assert(wL.rows() == size * size && wL.cols() == size * size + 2);
      L.emplace_back(wL);
    }
    assert(L.size() == - step);
  } else {
    L.reserve(step);
    for(int i = 0; i < step; i ++)
      L.emplace_back(SimpleMatrix<num_t>(size * size, size * size + 2).O());
    vector<vector<vector<SimpleVector<num_t> > > > cache;
    {
      vector<vector<SimpleVector<num_t> > > work;
      work.resize(L[0].rows());
      cache.resize(L.size(), work);
    }
    auto cache2(cache);
    for(int i = 5; i < argc; i ++) {
      for(int k = 0; k < L.size(); k ++) cerr << L[k] << std::endl;
      cerr << "remains: ";
      for(int k = i; k < argc; k ++) cerr << "\"" << argv[k] << "\" ";
      cerr << std::endl;
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      for(int kkk = 0; kkk < recur; kkk ++)
        for(int j = 0; j < work.size(); j ++) {
          cerr << kkk * work.size() + j << " / " << recur * work.size() << " over " << i - 5 << " / " << argc - 5 << std::endl;
          auto rin(work[j]);
          rin.O();
          for(int mm = 0; mm < L.size(); mm ++) {
            auto rin0(rin);
            for(int n = 0; n < rin.rows(); n ++)
              for(int nn = 0; nn < rin.cols(); nn ++)
                rin(n, nn) += rng();
            for(int k = 0; k < work[j].rows() - size; k ++)
              for(int kk = 0; kk < work[j].cols() - size; kk ++) {
                auto orig(work[j].subMatrix(k, kk, size, size) / num_t(int(4)) +
                          rin.subMatrix(k, kk, size, size));
                for(int mmm = 0; mmm < L[mm].rows(); mmm ++) {
                  SimpleVector<num_t> vwork(size * size + 1);
                  for(int n = 0; n < orig.rows(); n ++)
                    vwork.setVector(n * orig.cols(), orig.row(n));
                  vwork[size * size] =
                    - (work[j](k + (mmm / size), kk + (mmm % size)) +
                       rin0(k + (mmm / size), kk + (mmm % size)) );
                  for(int nnn = 0; nnn < vwork.size(); nnn ++)
                    vwork[nnn] = isfinite(vwork[nnn])
                      ? max(- num_t(int(1)), min(num_t(int(1)), vwork[nnn]))
                      : num_t(int(1)) / num_t(int(8));
                  auto mpi(makeProgramInvariant<num_t>(vwork));
  // XXX: Invariant summation causes average invariant.
  //      We need p1 or catg for linear ones,
  //      Otherwise we need multiplication and reduce method
  //      described in randtools.
  //      It is recommended to use latter one because of accuracy
  //      but they needs huge calculation time.
  //      Also, p1 and extreme accuracy condition should work but this needs
  //      huge memory usage.
                  cache[mm][mmm].emplace_back(move(mpi.first) *
                    pow(mpi.second, ceil(- log(orig.epsilon()) )) );
                  // L[mm].row(mmm) += move(mpi.first) *
                  //  pow(mpi.second, ceil(- log(orig.epsilon()) ));
                }
              }
          }
          for(int n = 0; n < L.size(); n ++) {
            num_t normL(int(0));
            for(int nn = 0; nn < L[n].rows(); nn ++) {
              SimpleMatrix<num_t> work(cache[n][nn].size(), cache[n][nn][0].size());
              for(int nnn = 0; nnn < work.rows(); nnn ++)
                work.row(nnn) = move(cache[n][nn][nnn]);
              cache2[n][nn].emplace_back(linearInvariant(work));
              cache[n][nn].resize(0);
            }
          }
        }
    }
    for(int n = 0; n < L.size(); n ++) {
      num_t normL(int(0));
      for(int nn = 0; nn < L[n].rows(); nn ++) {
        SimpleMatrix<num_t> work(cache2[n][nn].size(), cache2[n][nn][0].size());
        for(int nnn = 0; nnn < work.rows(); nnn ++)
          work.row(nnn) = move(cache2[n][nn][nnn]);
        L[n].row(nn)  = linearInvariant(work);
        L[n].row(nn) /= - num_t(L[n](nn, L[n].cols() - 2));
        L[n](nn, L[n].cols() - 2) = num_t(int(0));
        normL += L[n].row(nn).dot(L[n].row(nn));
        cache2[n][nn].resize(0);
      }
      // XXX: don't know why, but *= 2 scales well on revert.
      //L[n] /= sqrt(normL /= num_t(int(L[n].rows()))) * num_t(int(2));
      L[n] /= sqrt(normL /= num_t(int(L[n].rows())));
      std::cout << L[n] << std::endl;
    }
  }
  auto outc(out);
  for(int i = 0; i < outc.size(); i ++) outc[i].O();
  auto outr(outc);
  SimpleMatrix<num_t> one(size, size);
  one.O(num_t(int(1)));
  if(recur0 < 0) for(int i = 0; i < out.size(); i ++) out[i].O();
  for(int rc = 0; rc < (size0 < 0 ? 1 : recur); rc ++) {
    cerr << rc << " / " << recur << std::endl;
    auto rin(out[0]);
    rin.O();
    if(0 < size0)
      for(int n = 0; n < rin.rows(); n ++)
        for(int nn = 0; nn < rin.cols(); nn ++)
          for(int nnn = 0; nnn < L.size() - 1; nnn ++)
            rin(n, nn) += rng();
    if(! rc) {
      vector<SimpleMatrix<num_t> > buf;
      for(int i = 0; i < out.size(); i ++) buf.emplace_back(out[i] / num_t(int(4)) + rin);
      savep2or3<num_t>((std::string(argv[4]) + std::string("-oneshot.ppm")).c_str(), buf, false, 65535);
    }
    for(int j = 0; j < out.size(); j ++)
      for(int k = 0; k < out[j].rows() - size; k ++)
        for(int kk = 0; kk < out[j].cols() - size; kk ++) {
          auto orig(out[j].subMatrix(k, kk, size, size) / num_t(int(4)) +
                    rin.subMatrix(k, kk, size, size));
          SimpleVector<num_t> vwork(size * size + 1);
          vwork.O();
          for(int n = 0; n < orig.rows(); n ++)
            vwork.setVector(n * orig.cols(), orig.row(n));
          for(int nn = L.size() - 1; 0 <= nn; nn --) {
            for(int nnn = 0; nnn < vwork.size(); nnn ++)
              vwork[nnn] = isfinite(vwork[nnn])
                ? max(- num_t(int(1)), min(num_t(int(1)), vwork[nnn]))
                : num_t(int(1)) / num_t(int(8));
            vwork[vwork.size() - 1] = - num_t(int(1));
            auto mpi(makeProgramInvariant<num_t>(vwork));
            vwork = L[nn] * move(mpi.first);
            for(int nnn = 0; nnn < vwork.size(); nnn ++) 
              // XXX: don't know why this needs + 1,
              //      but pair.second must be 1 because of scaling.
              // vwork[nnn] = revertProgramInvariant<num_t>(make_pair(vwork[nnn] + num_t(int(1)), num_t(int(1)) ));
              vwork[nnn] = revertProgramInvariant<num_t>(make_pair(vwork[nnn], num_t(int(1)) ));
            vwork = SimpleVector<num_t>(size * size + 1).O().setVector(0, vwork);
          }
          for(int nnn = 0; nnn < vwork.size(); nnn ++)
            vwork[nnn] = isfinite(vwork[nnn])
              ? max(num_t(int(0)), min(num_t(int(1)), abs(vwork[nnn])))
              : num_t(int(1)) / num_t(int(2));
          SimpleMatrix<num_t> temp(size, size);
          for(int n = 0; n < temp.rows(); n ++)
            temp.row(n) = vwork.subVector(n * size, size);
          outr[j].setMatrix(k, kk, outr[j].subMatrix(k, kk, size, size) + temp);
          outc[j].setMatrix(k, kk, outc[j].subMatrix(k, kk, size, size) + one);
        }
  }
  for(int i = 0; i < outr.size(); i ++)
    for(int k = 0; k < outr[i].rows(); k ++)
      for(int kk = 0; kk < outr[i].cols(); kk ++)
        if(outc[i](k, kk) != num_t(int(0))) outr[i](k, kk) /= outc[i](k, kk);
  if(! savep2or3<num_t>(argv[4], normalize<num_t>(autoLevel<num_t>(outr, (outr[0].rows() + outr[0].cols()) * 9)), false, 65535)) return - 2;
  return 0;
}

