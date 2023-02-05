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
#include "p2.hh"

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
          data.resize(1);
          data[0] = SimpleMatrix<T>(h, w).O();
          loadstub<T>(input, nmax, 1, data);
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
              iik <= min(in[i][j].rows() / sz - 1,
                in[i][j].rows() - ii * (in[i][j].rows() / sz)); iik ++)
            for(int jjk = 0; jjk <= min(in[i][j].cols() / sz - 1,
                  in[i][j].cols() - jj * (in[i][j].cols() / sz));
                jjk ++, cnt ++)
              shrink[i][j](ii, jj) +=
                in[i][j](ii * (in[i][j].rows() / sz) + iik,
                         jj * (in[i][j].cols() / sz) + jjk);
          if(cnt) shrink[i][j](ii, jj) /= num_t(cnt);
        }
  }
  return shrink;
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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  vector<vector<SimpleVector<num_t> > > in;
  int rows(0);
  int cols(0);
  in.reserve(argc - 1);
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    in.emplace_back(vector<SimpleVector<num_t> >());
    in[i - 1].resize(work.size());
    for(int j = 0; j < work.size(); j ++) {
      in[i - 1][j].resize(work[j].rows() * work[j].cols());
      for(int k = 0; k < work[j].rows(); k ++)
        in[i - 1][j].setVector(k * work[j].cols(), work[j].row(k));
      rows = work[j].rows();
      cols = work[j].cols();
      assert(in[0][0].size() == in[i - 1][j].size());
    }
    assert(in[0].size() == work.size());
  }
  vector<vector<SimpleVector<num_t> > > out;
  out.resize(in[0].size());
  for(int i = 0; i < out.size(); i ++) {
    vector<SimpleVector<num_t> > d;
    d.reserve(in.size());
    for(int k = 0; k < in.size(); k ++)
      d.emplace_back(move(in[k][i]));
    auto p(predv(d));
    out[i] = std::move(p.first);
    out[i].insert(out[i].end(), p.second.begin(), p.second.end());
  }
  vector<SimpleMatrix<num_t> > outs;
  outs.resize(out.size());
  for(int i = 0; i < out[0].size(); i ++) {
    for(int j = 0; j < out.size(); j ++) {
      outs[j].resize(rows, cols);
      for(int k = 0; k < outs[j].rows() * outs[j].cols(); k ++)
        outs[j](k / outs[j].cols(), k % outs[j].cols()) =
          revertProgramInvariant<num_t>(make_pair(
            out[j][i][k] / out[j][i][out[j][i].size() - 1], num_t(int(1)) ));
    }
    if(! savep2or3<num_t>((std::string("predg-") + std::to_string(i) + std::string(".ppm")).c_str(), normalize<num_t>(outs), outs.size() == 1, 65535) )
      std::cerr << "failed to save." << std::endl;
  }
  cerr << " Done" << std::endl;
  return 0;
}

