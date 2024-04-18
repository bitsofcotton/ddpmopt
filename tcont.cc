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

#define int int64_t
//#define int int32_t
#include "lieonn.hh"
typedef myfloat num_t;

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::move;

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
//#define int int32_t
  assert(1 < argc);
  vector<num_t> score;
  score.resize(argc + 1, num_t(int(0)));
  switch(argv[1][0]) {
  case 'x':
  case 'i':
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      P0maxRank<num_t> p;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      for(int i = 0; i < work.size(); i ++)
        for(int ii = 0; ii < work[i].rows(); ii ++) {
          idFeeder<num_t> w(3);
          for(int jj = 0; jj < work[i].cols(); jj ++) {
            if(w.full) {
              const auto pp(p.next(w.res));
              score[i0] += (pp - work[i](ii, jj)) * (pp - work[i](ii, jj));
            } 
            w.next(work[i](ii, jj));
          }
        }
      if(argv[1][0] == 'x')
        score[i0] /= num_t(work[0].rows() * work[0].cols() * work.size());
    }
    if(argv[1][0] == 'x') break;
  case 'y':
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      P0maxRank<num_t> p;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      for(int i = 0; i < work.size(); i ++)
        for(int jj = 0; jj < work[i].cols(); jj ++) {
          idFeeder<num_t> w(3);
          for(int ii = 0; ii < work[i].rows(); ii ++) {
            if(w.full) {
              const auto pp(p.next(w.res));
              score[i0] += (pp - work[i](ii, jj)) * (pp - work[i](ii, jj));
            } 
            w.next(work[i](ii, jj));
          }
        }
      score[i0] /= num_t(work[0].rows() * work[0].cols() * work.size());
    }
    break;
  case 't':
    {
      vector<vector<SimpleMatrix<num_t> > > b;
      b.resize(argc + 1);
      for(int i = 2; i < argc; i ++) {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i])) continue;
        b[i] = move(work);
        assert(b[i].size() == b[2].size() &&
          b[i][0].rows() == b[2][0].rows() &&
          b[i][0].cols() == b[2][0].cols());
      }
      P0maxRank<num_t> p;
      int cnt(0);
      for(int i = 0; i < b[2][0].rows(); i ++)
        for(int j = 0; j < b[2][0].cols(); j ++)
          for(int k = 0; k < b[2].size(); k ++) {
            idFeeder<num_t> w(3);
            for(int ii = 2; ii < b.size(); ii ++, cnt ++) {
              if(! b[ii].size()) continue;
              if(w.full) {
                const auto pp(p.next(w.res));
                score[0] += (pp - b[ii][k](i, j)) * (pp - b[ii][k](i, j));
              }
              w.next(b[ii][k](i, j));
            }
          }
      score[0] /= num_t(cnt);
    }
    break;
  default:
    cerr << "Invalid argv[1]" << endl;
  }
  if(argv[1][0] == 'x' || argv[1][0] == 'y' || argv[1][0] == 'i')
    for(int i = 2; i < argc; i ++)
      cout << sqrt(score[i]) << ", " << argv[i] << endl;
  else
    cout << sqrt(score[0]) << ", whole image index" << endl;
  return 0;
}

