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
  return max(- num_t(int(1)), min(num_t(int(1)), (num_t(res) / num_t(~ myuint(0)) - num_t(int(1)) / num_t(int(2))) * num_t(int(2)) ));
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
  const auto epoch0(std::atoi(argv[1]));
        auto epoch(abs(epoch0));
  if(epoch0 < 0) {
    vector<SimpleMatrix<num_t> > L;
    L.reserve(3);
    int rows(0);
    std::cin >> rows;
    assert(0 < rows);
    for(int j = 0; j < 3; j ++) {
      SimpleMatrix<num_t> wL(rows, rows + 2);
      for(int i = 0; i < wL.rows(); i ++)
        // XXX: clang++13 causes SimpleFloat ... operator >> while loop
        //      illegaly omitted causes zero vector for each.
        //      however, this can be only the clang binary set I got had be
        //      infected condition.
        std::cin >> wL.row(i);
      // L.emplace_back(move(wL));
      L.emplace_back(wL);
      assert(L[0].rows() == L[j].rows() && L[0].cols() == L[j].cols());
      assert(L[j].rows() + 2 == L[j].cols());
    }
    epoch *= L[0].rows();
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > out;
      if(! loadp2or3<num_t>(out, argv[i])) return - 1;
      assert(out[0].rows() * out[0].cols() == L[0].rows());
      auto outs(out);
      for(int j = 0; j < outs.size(); j ++) outs[j].O();
      for(int k = 0; k < epoch; k ++) {
        auto rin(out[0]);
        for(int n = 0; n < rin.rows(); n ++)
          for(int nn = 0; nn < rin.cols(); nn ++)
            rin(n, nn) = rng();
        for(int j = 0; j < out.size(); j ++) {
          cerr << k * out.size() + j << " / " << epoch * out.size() << " over " << i - 2 << " / " << argc - 2 << std::endl;
          SimpleVector<num_t> vwork(out[j].rows() * out[j].cols() + 1);
          for(int n = 0; n < out[j].rows(); n ++)
            vwork.setVector(n * out[j].cols(), (out[j].row(n) + rin.row(n)) / num_t(int(2)) );
          vwork[vwork.size() - 1] = num_t(int(0));
          auto outwork(L[j] * makeProgramInvariant<num_t>(vwork).first);
          for(int n = 0; n < outwork.size(); n ++)
            outwork[n] = revertProgramInvariant<num_t>(make_pair(outwork[n], num_t(int(1)) ));
          for(int n = 0; n < out[j].rows(); n ++)
            outs[j].row(n) += outwork.subVector(n * out[j].cols(), out[j].cols());
        }
      }
      if(! savep2or3<num_t>(argv[i], normalize<num_t>(outs), false, 65535) )
        std::cerr << "failed to save." << std::endl;
    }
  } else {
    vector<vector<SimpleMatrix<num_t> > > in;
    vector<vector<SimpleMatrix<num_t> > > noise;
    in.resize(argc - 2);
    noise.resize(in.size());
    for(int i = 2; i < argc; i ++) {
      if(! loadp2or3<num_t>(in[i - 2], argv[i])) continue;
      assert(in[0][0].rows() == in[i - 2][0].rows() &&
             in[0][0].cols() == in[i - 2][0].cols());
      assert(in[i - 2][0].rows() == in[i - 2][0].cols());
      if(i == 2) epoch *= max(1, in[0][0].rows() * in[0][0].cols() / (argc - 2));
      noise[i - 2].resize(epoch, SimpleMatrix<num_t>(in[i - 2][0].rows(), in[i - 2][0].cols()));
      for(int j = 0; j < epoch; j ++)
        for(int n = 0; n < noise[i - 2][j].rows(); n ++)
          for(int nn = 0; nn < noise[i - 2][j].cols(); nn ++)
            noise[i - 2][j](n, nn) = rng();
    }
    std::cout << in[0][0].rows() * in[0][0].cols() << std::endl;
    for(int j = 0; j < in[0].size(); j ++)
      for(int m = 0; m < in[0][0].rows() * in[0][0].cols(); m ++) {
        cerr << j * in[0][0].rows() * in[0][0].cols() + m << " / " << in[0][0].rows() * in[0][0].cols() * in[0].size() << std::endl;
        SimpleMatrix<num_t> work(epoch * in.size(), in[0][0].rows() * in[0][0].cols() + 2);
        for(int i = 0; i < in.size(); i ++)
          for(int jj = 0; jj < epoch; jj ++) {
            SimpleVector<num_t> vwork(in[i][j].rows() * in[i][j].cols() + 1);
            for(int n = 0; n < in[i][j].rows(); n ++)
              vwork.setVector(n * in[i][j].cols(), (in[i][j].row(n) + noise[i][jj].row(n)) / num_t(int(2)) );
            vwork[vwork.size() - 1] = in[i][j](m / in[i][j].cols(), m % in[i][j].cols());
  // XXX: Invariant summation causes average invariant.
  //      We need p1 or catg for linear ones,
  //      Otherwise we need multiplication and reduce method
  //      described in randtools.
  //      It is recommended to use latter one because of accuracy
  //      but they needs huge calculation time.
  //      Also, p1 and extreme accuracy condition should work but this needs
  //      huge memory usage.
  // XXX: in the both case, we need to treat different corresponding as
  //      different vectors. Otherwise, we get spoiled result.
            auto mpi(makeProgramInvariant<num_t>(vwork));
            work.row(i * epoch + jj)  = move(mpi.first);
            work.row(i * epoch + jj) *=
              pow(mpi.second, ceil(- log(in[0][0].epsilon()) ));
          }
        auto vwork(linearInvariant(work));
        vwork /= - num_t(vwork[vwork.size() - 2]);
        vwork[vwork.size() - 2] = num_t(int(0));
        std::cout << vwork;
      }
  }
  return 0;
}

