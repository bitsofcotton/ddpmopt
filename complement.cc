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

#if defined(_FLOAT_BITS_)
num_t edge(0);
static inline num_t rng(std::ranlux48& rde) {
  myuint res(0);
  for(int i = 0; i < _FLOAT_BITS_ / sizeof(uint32_t) / 8; i ++) {
    res <<= sizeof(uint32_t) * 8;
    res  |= uint32_t(rde());
  }
  return num_t(res) / num_t(~ myuint(0)) * edge;
}
#endif

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  vector<SimpleMatrix<num_t> > L;
  const auto step(std::atoi(argv[1]));
  int size(0);
  const auto recur(std::atoi(argv[2]));
  vector<SimpleMatrix<num_t> > out;
  SimpleMatrix<num_t> mask;
  if(! loadp2or3<num_t>(out, argv[4])) return - 1;
  mask = out[0];
  const auto mask0(mask);
  out = vector<SimpleMatrix<num_t> >();
  if(! loadp2or3<num_t>(out, argv[3])) return - 1;
  L.reserve(step);
  for(int i = 0; i < step; i ++) {
    SimpleMatrix<num_t> wL;
    std::cin >> wL;
    if(! i) size = int(sqrt(num_t(int(wL.rows()))));
    assert(wL.rows() == size * size && wL.cols() == size * size + 2);
    L.emplace_back(wL);
  }
  assert(L.size() == step);
  std::random_device rd;
  std::ranlux48 rde(rd());
#if defined(_FLOAT_BITS_)
  edge = num_t(2) / num_t(4 * int(abs(step) + 1));
#else
  std::uniform_real_distribution<num_t> rng(- num_t(1) / num_t(4 * int(abs(step) + 1)), num_t(1) / num_t(4 * int(abs(step) + 1)) );
#endif
  auto outc(out);
  for(int i = 0; i < outc.size(); i ++) outc[i].O();
  auto outr(outc);
  SimpleMatrix<num_t> one(size, size);
  one.O(num_t(int(1)));
  out = normalize<num_t>(autoLevel<num_t>(out, (out[0].rows() + out[0].cols()) * 3));
  for(int i = 0; i < mask.rows(); i ++)
    for(int j = 0; j < mask.cols(); j ++)
      if(mask0(i, j) < num_t(int(1)) / num_t(int(2)) &&
                       num_t(int(1)) / num_t(int(2)) <= mask(i, j))
        for(int k = 0; k < out.size(); k ++) out[k](i, j) = num_t(int(1)) / num_t(int(2));
  while(1) {
    bool cont(true);
    for(int i = 0; i < mask.rows(); i ++)
      for(int j = 0; j < mask.cols(); j ++)
        if(mask0(i, j) < num_t(int(1)) / num_t(int(2)) &&
                         num_t(int(1)) / num_t(int(2)) <= mask(i, j))
          for(int k = 0; k < out.size(); k ++)out[k](i, j) = outr[k](i, j);
        else if(mask0(i, j) < num_t(int(1)) / num_t(int(2)) )
          cont = false;
    if(cont) break;
    for(int i = 0; i < outc.size(); i ++) outc[i].O();
    for(int i = 0; i < outr.size(); i ++) outr[i].O();
    for(int rc = 0; rc < recur; rc ++) {
      cerr << rc << " / " << recur << std::endl;
      auto rin(out[0]);
      rin.O();
      for(int n = 0; n < rin.rows(); n ++)
        for(int nn = 0; nn < rin.cols(); nn ++)
          for(int nnn = 0; nnn < L.size() - 1; nnn ++)
            rin(n, nn) += rng(rde);
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
                  : vwork[nnn] = num_t(int(1)) / num_t(int(8));
              vwork[vwork.size() - 1] = - num_t(int(1));
              auto mpi(makeProgramInvariant<num_t>(vwork));
              vwork = L[nn] * move(mpi.first);
              for(int nnn = 0; nnn < vwork.size(); nnn ++)
                vwork[nnn] = revertProgramInvariant<num_t>(make_pair(vwork[nnn], mpi.second)) / pow(mpi.second, ceil(- log(orig.epsilon()) ));
              vwork = SimpleVector<num_t>(size * size + 1).O().setVector(0, vwork);
            }
            SimpleMatrix<num_t> temp(size, size);
            for(int n = 0; n < temp.rows(); n ++)
              temp.row(n) = vwork.subVector(n * size, size);
            outr[j].setMatrix(k, kk, outr[j].subMatrix(k, kk, size, size) + temp);
            outc[j].setMatrix(k, kk, outc[j].subMatrix(k, kk, size, size) + one);
          }
    }
    for(int i = 0; i < outr.size(); i ++)
      for(int k = 0; k < outr[i].rows(); k ++)
        for(int kk = 0; kk < outr[i].cols(); kk ++) {
          if(outc[i](k, kk) != num_t(int(0))) outr[i](k, kk) /= outc[i](k, kk);
          outr[i](k, kk) = max(num_t(int(0)), min(num_t(int(1)), outr[i](k, kk)));
        }
    outr = normalize<num_t>(autoLevel<num_t>(outr, (outr[0].rows() + outr[0].cols()) * 3));
    auto mask2(mask);
    for(int i = 0; i < mask.rows(); i ++)
      for(int j = 0; j < mask.cols(); j ++)
        if(num_t(int(1)) / num_t(int(2)) <= mask(i, j))
          mask2(max(0, i - 1), j) =
            mask2(i, max(0, j - 1)) =
            mask2(min(i, mask2.rows() - 1), j) =
            mask2(i, min(j, mask2.cols() - 1)) = num_t(int(1));
    mask = mask2;
  }
  if(! savep2or3<num_t>(argv[4], outr, false, 65535)) return - 2;
  return 0;
}

