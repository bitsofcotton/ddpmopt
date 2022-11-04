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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  const auto m(argv[1][0]);
  if(m == '-') {
    vector<SimpleMatrix<num_t> > L;
    int sz0(0);
    std::cin >> sz0;
    assert(0 < sz0);
    L.reserve(3 * sz0 * sz0);
    for(int j = 0; j < 3; j ++)
      for(int kk = 0; kk < sz0 * sz0; kk ++) {
        SimpleMatrix<num_t> wL(sz0 * sz0, sz0 * sz0 + 2);
        for(int i = 0; i < wL.rows(); i ++)
          std::cin >> wL.row(i);
        // L.emplace_back(move(wL));
        L.emplace_back(wL);
        assert(L[0].rows() == L[j].rows() && L[0].cols() == L[j].cols());
        assert(L[j].rows() + 2 == L[j].cols());
      }
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > out;
      if(! loadp2or3<num_t>(out, argv[i])) return - 1;
      assert(out[0].rows() * out[0].cols() == L[0].rows());
      auto outs(out);
      for(int n = 0; n < outs.size(); n ++)
        outs[n] = SimpleMatrix<num_t>(sz0 * sz0, sz0 * sz0);
      for(int j = 0; j < out.size(); j ++) {
        cerr << j << " / " << out.size() << " over " << i - 2 << " / " << argc - 2 << std::endl;
        SimpleVector<num_t> vwork0(out[j].rows() * out[j].cols() + 1);
        for(int n = 0; n < out[j].rows(); n ++)
          vwork0.setVector(n * out[j].cols(), out[j].row(n));
        vwork0[vwork0.size() - 1] = num_t(int(0));
        const auto vwork(makeProgramInvariant<num_t>(vwork0).first);
        for(int k = 0; k < sz0 * sz0; k ++) {
          auto outwork(L[j * sz0 * sz0 + k] * vwork);
          for(int n = 0; n < outwork.size(); n ++)
            outwork[n] = revertProgramInvariant<num_t>(make_pair(outwork[n], num_t(int(1)) ));
          for(int n = 0; n < out[j].rows(); n ++)
            out[j].row(n) = outwork.subVector(n * out[j].cols(), out[j].cols());
          for(int n = 0; n < out[j].rows(); n ++)
            outs[j].setMatrix((k / sz0) * sz0, (k % sz0) * sz0, out[j]);
        }
      }
      if(! savep2or3<num_t>(argv[i], normalize<num_t>(outs), false, 65535) )
        std::cerr << "failed to save." << std::endl;
    }
  } else if(m == '+') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.resize(argc - 2);
    for(int i = 2; i < argc; i ++) {
      if(! loadp2or3<num_t>(in[i - 2], argv[i])) continue;
      assert(in[0][0].rows() == in[i - 2][0].rows() &&
             in[0][0].cols() == in[i - 2][0].cols());
      assert(in[i - 2][0].rows() == in[i - 2][0].cols());
    }
    const int sz0(pow(num_t(in[0][0].rows()), num_t(int(1)) / num_t(int(3))));
    assert(sz0 * sz0 == in[0][0].rows() && sz0 * sz0 == in[0][0].cols());
    assert(sz0 * sz0 + 2 < in.size() &&
           "Input number of graphics is too small differed from input graphics"
           "size.");
    std::cout << sz0 << std::endl;
    auto shrink(in);
    for(int i = 0; i < shrink.size(); i ++)
      for(int j = 0; j < shrink[i].size(); j ++) {
        shrink[i][j] = SimpleMatrix<num_t>(sz0, sz0).O();
        for(int ii = 0; ii < sz0; ii ++)
          for(int jj = 0; jj < sz0; jj ++)
            for(int iik = 0; iik < sz0; iik ++)
              for(int jjk = 0; jjk < sz0; jjk ++)
                shrink[i][j](ii, jj) += in[i][j](ii * sz0 + iik, jj * sz0 + jjk);
        shrink[i][j] /= num_t(sz0 * sz0);
      }
    for(int j = 0; j < in[0].size(); j ++)
      for(int m = 0; m < in[0][0].rows() * in[0][0].cols(); m ++) {
        cerr << j * in[0][0].rows() * in[0][0].cols() + m << " / " << in[0][0].rows() * in[0][0].cols() * in[0].size() << std::endl;
        SimpleMatrix<num_t> work(in.size(), shrink[0][0].rows() * shrink[0][0].cols() + 2);
        for(int i = 0; i < in.size(); i ++) {
          SimpleVector<num_t> vwork(shrink[i][j].rows() * shrink[i][j].cols() + 1);
          for(int n = 0; n < shrink[i][j].rows(); n ++)
            vwork.setVector(n * shrink[i][j].cols(), shrink[i][j].row(n));
          vwork[vwork.size() - 1] = in[i][j](m / in[i][j].cols(), m % in[i][j].cols());
          auto mpi(makeProgramInvariant<num_t>(vwork));
          work.row(i)  = move(mpi.first);
          work.row(i) *=
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

