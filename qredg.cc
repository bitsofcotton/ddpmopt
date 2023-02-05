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
#include "util.hh"
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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    for(int j = 0; j < work.size(); j ++) {
      vector<SimpleVector<num_t> > pwork;
      pwork.reserve(work[j].rows());
      for(int k = 0; k < work[j].rows(); k ++)
        pwork.emplace_back(work[j].row(k));
      auto py(predv<num_t>(pwork));
      vector<SimpleVector<num_t> > qwork;
      qwork.reserve(work[j].cols());
      for(int k = 0; k < work[j].cols(); k ++)
        qwork.emplace_back(work[j].col(k));
      auto px(predv<num_t>(qwork));
      SimpleMatrix<num_t> swork(work[j].rows() + py.first.size() * 2,
                                work[j].cols() + px.first.size() * 2);
      swork.O().setMatrix(py.first.size(), px.first.size(), work[j]);
      for(int k = 0; k < py.first.size(); k ++)
        for(int kk = 0; kk < py.first[k].size(); kk ++) {
          swork(py.first.size() - k - 1, px.first.size() + kk) =
            revertProgramInvariant<num_t>(make_pair(
              py.first[k][kk] / py.first[k][py.first[k].size() - 1],
              num_t(int(1)) ));
          swork(py.first.size() + work[j].rows() + k, px.first.size() + kk) =
            revertProgramInvariant<num_t>(make_pair(
              py.second[k][kk] / py.second[k][py.second[k].size() - 1],
              num_t(int(1)) ));
        }
      for(int k = 0; k < px.first.size(); k ++)
        for(int kk = 0; kk < px.first[k].size(); kk ++) {
          swork(py.first.size() + kk, px.first.size() - k - 1) =
            revertProgramInvariant<num_t>(make_pair(
              px.first[k][kk] / px.first[k][px.first[k].size() - 1],
              num_t(int(1)) ));
          swork(py.first.size() + kk, px.first.size() + work[j].cols() + k) =
            revertProgramInvariant<num_t>(make_pair(
              px.second[k][kk] / px.second[k][px.second[k].size() - 1],
              num_t(int(1)) ));
        }
      std::swap(work[j], swork);
    }
    if(! savep2or3<num_t>(argv[i], work) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

