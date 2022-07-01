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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  vector<SimpleMatrix<num_t> > L;
  const auto step(std::atoi(argv[1]));
  const auto size(std::atoi(argv[2]));
  const auto recur(std::atoi(argv[3]));
  vector<SimpleMatrix<num_t> > out;
  if(! loadp2or3<num_t>(out, argv[4]))
    return - 1;
  L.reserve(step);
  for(int i = 0; i < step; i ++)
    L.emplace_back(SimpleMatrix<num_t>(size * size, size * size * 2 + 1).O());
  std::random_device rd;
  std::ranlux48 rde(rd());
  std::normal_distribution<num_t> rng(num_t(0), num_t(1));
  for(int i = 5; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    for(int j = 0; j < work.size(); j ++)
      for(int k = 0; k < work[j].rows() - size; k ++) {
        std::cerr << k << " / " << work[j].rows() - size << " over " << i - 5 << " / " << argc - 5 << std::endl;
        for(int kk = 0; kk < work[j].cols() - size; kk ++) {
          auto orig(work[j].subMatrix(k, kk, size, size));
          for(int m = 0; m < L[0].rows(); m ++)
            for(int kkk = 0; kkk < recur; kkk ++) {
              SimpleVector<num_t> vwork0(size * size);
              for(int n = 0; n < orig.rows(); n ++)
                vwork0.setVector(n * orig.cols(), orig.row(n));
              for(int mm = 0; mm < L.size(); mm ++) {
                vwork0 /= sqrt(vwork0.dot(vwork0));
                SimpleVector<num_t> vwork(size * size * 2);
                vwork.setVector(0, vwork0);
                for(int n = 0; n < vwork0.size(); n ++)
                  vwork0[n] += rng(rde) / num_t(int(L.size() + 1));
                vwork.setVector(vwork0.size(), vwork0);
                vwork /= sqrt(vwork.dot(vwork));
                auto mpi(makeProgramInvariant<num_t>(vwork));
                L[mm].row(m) += move(mpi.first) * pow(mpi.second, ceil(- log(orig.epsilon()) ));
              }
            }
        }
     }
  }
  auto outc(out);
  auto outr(out);
  for(int i = 0; i < outc.size(); i ++) outc[i].O();
  for(int i = 0; i < outr.size(); i ++) outr[i].O();
  SimpleMatrix<num_t> one(size, size);
  one.O(num_t(int(1)));
  for(int j = 0; j < out.size(); j ++)
    for(int k = 0; k < out[j].rows() - size; k ++) {
      std::cerr << k << " / " << out[j].rows() - size << std::endl;
      for(int kk = 0; kk < out[j].cols() - size; kk ++) {
        auto orig(out[j].subMatrix(k, kk, size, size));
        SimpleVector<num_t> vwork(size * size * 2);
        vwork.O();
        for(int n = 0; n < orig.rows(); n ++)
          vwork.setVector(size * size + n * orig.cols(), orig.row(n));
        for(int nn = L.size() - 1; 0 <= nn; nn --) {
          auto mpi(makeProgramInvariant<num_t>(vwork));
          vwork = move(mpi.first) * pow(mpi.second, ceil(- log(orig.epsilon()) ));
          vwork = L[nn].subMatrix(0, 0, size * size, size * size).solve(- L[nn].subMatrix(0, size * size, size * size, size * size + 1) * vwork.subVector(size * size, size * size + 1));
          SimpleVector<num_t> work(size * size * 2);
          work.O();
          for(int m = 0; m < size * size; m ++)
            work[size * size + m] = max(num_t(0), min(num_t(1), revertProgramInvariant<num_t>(make_pair(vwork[m], mpi.second)) / pow(mpi.second, ceil(- log(orig.epsilon()) )) ));
          vwork = work;
        }
        SimpleMatrix<num_t> temp(size, size);
        for(int n = 0; n < temp.rows(); n ++)
          temp.row(n) = vwork.subVector(n * size + size * size, size);
        outr[j].setMatrix(k, kk, outr[j].subMatrix(k, kk, size, size) + temp);
        outc[j].setMatrix(k, kk, outc[j].subMatrix(k, kk, size, size) + one);
      }
    }
  for(int i = 0; i < outr.size(); i ++)
    for(int k = 0; k < outr[i].rows(); k ++)
      for(int kk = 0; kk < outr[i].cols(); kk ++)
        outr[i](k, kk) /= outc[i](k, kk);
  if(! savep2or3<num_t>(argv[4], outr, true, 65535))
    return - 2;
  return 0;
}

