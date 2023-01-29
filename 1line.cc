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
  const auto m(std::atoi(argv[1]));
  if(m) {
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i])) {
        cerr << "could not open " << argv[i] << std::endl;
        continue;
      }
      if(in[0].rows() != 1) {
        cerr << argv[i] << " : isn't 1-lined" << std::endl;
        continue;
      }
      if(in[0].cols() % abs(m) != 0) {
        cerr << argv[i] << " : argv[1] isn't proper" << std::endl;
        continue;
      }
      vector<SimpleMatrix<num_t> > out;
      out.reserve(m < 0 ? 1 : 3);
      for(int ii = 0; ii < in.size(); ii ++) {
        out.emplace_back(SimpleMatrix<num_t>(abs(m), in[ii].cols() / abs(m) / (m < 0 ? 1 : 3)));
        for(int j = 0; j < abs(m); j ++)
          out[ii].row(j) = in[ii].subMatrix(0, out[ii].cols() * j +
            out[0].rows() * out[0].cols() * ii, i, out[ii].cols()).row(0);
      }
      if(! savep2or3<num_t>(argv[i], out, ! (m < 0), 65535) )
        cerr << "could not save " << argv[i] << std::endl;
    }
  } else {
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i])) {
        cerr << "could not open " << argv[i] << std::endl;
        continue;
      }
      vector<SimpleMatrix<num_t> > out;
      out.emplace_back(SimpleMatrix<num_t>(1, in[0].rows() * in[0].cols() * in.size()));
      for(int ii = 0; ii < in.size(); ii ++)
        for(int j = 0; j < in[ii].rows(); j ++)
          out[0].row(0).setVector(j * in[ii].cols() + in[0].cols() * in[0].rows() * ii, in[ii].row(j));
      if(! savep2or3<num_t>(argv[i], out, true, 65535) )
        cerr << "could not save " << argv[i] << std::endl;
    }
  }
  return 0;
}

