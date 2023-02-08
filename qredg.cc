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
    vector<SimpleVector<num_t> > pwork;
    vector<SimpleVector<num_t> > qwork;
    pwork.resize(work[0].rows(),
      SimpleVector<num_t>(work[0].cols() * work.size()).O());
    qwork.resize(work[0].cols(),
      SimpleVector<num_t>(work[0].rows() * work.size()).O());
    for(int j = 0; j < work.size(); j ++)
      for(int k = 0; k < work[j].rows(); k ++)
        pwork[k].setVector(j * work[j].rows() * work[j].cols() +
                           k * work[j].cols(), work[j].row(k));
    for(int j = 0; j < work.size(); j ++)
      for(int k = 0; k < work[j].cols(); k ++)
        qwork[k].setVector(j * work[j].rows() * work[j].cols() +
                           k * work[j].rows(), work[j].col(k));
    auto py(predv<num_t>(pwork));
    auto px(predv<num_t>(qwork));
    vector<SimpleMatrix<num_t> > swork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + py.first.size() * 2,
                          work[0].cols() + px.first.size() * 2).O());
    for(int j = 0; j < work.size(); j ++)
      swork[j].O().setMatrix(py.first.size(), px.first.size(), work[j]);
    for(int k = 0; k < py.first.size(); k ++)
      for(int kk = 0; kk < py.first[k].size() / work.size(); kk ++)
        for(int j = 0; j < work.size(); j ++) {
          swork[j](py.first.size() - k - 1, px.first.size() + kk) =
            revertProgramInvariant<num_t>(make_pair(
              py.first[k][j * py.first[k].size() / work.size() + kk] /
                py.first[k][py.first[k].size() - 1],
            num_t(int(1)) ));
          swork[j](py.first.size() + work[j].rows() + k, px.first.size() + kk) =
            revertProgramInvariant<num_t>(make_pair(
              py.second[k][j * py.first[k].size() / work.size() + kk] /
                py.second[k][py.second[k].size() - 1],
            num_t(int(1)) ));
        }
    for(int k = 0; k < px.first.size(); k ++)
      for(int kk = 0; kk < px.first[k].size() / work.size(); kk ++)
        for(int j = 0; j < work.size(); j ++) {
          swork[j](py.first.size() + kk, px.first.size() - k - 1) =
            revertProgramInvariant<num_t>(make_pair(
              px.first[k][j * px.first[k].size() / work.size() + kk] /
                px.first[k][px.first[k].size() - 1],
              num_t(int(1)) ));
          swork[j](py.first.size() + kk, px.first.size() + work[j].cols() + k) =
            revertProgramInvariant<num_t>(make_pair(
              px.second[k][j * px.first[k].size() / work.size() + kk] /
                px.second[k][px.second[k].size() - 1],
              num_t(int(1)) ));
        }
    if(! savep2or3<num_t>(argv[i], normalize<num_t>(swork)) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

