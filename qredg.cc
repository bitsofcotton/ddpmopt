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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(2 < argc);
  for(int i = 2; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    vector<SimpleVector<num_t> > pwork;
    pwork.resize(work[0].rows(),
      SimpleVector<num_t>(work[0].cols() * work.size()).O());
    for(int k = 0; k < work[0].rows(); k ++)
      for(int j = 0; j < work.size(); j ++)
        pwork[k].setVector(j * work[j].cols(), work[j].row(k));
    auto p(predv<num_t>(pwork, std::atoi(argv[1]) ));
    vector<SimpleMatrix<num_t> > swork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + p.first.size() * 2,
        work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      swork[j].setMatrix(p.first.size(), 0, work[j] / num_t(int(2)));
    for(int k = 0; k < p.first.size(); k ++)
      for(int kk = 0; kk < p.first[k].size() / work.size(); kk ++)
        for(int j = 0; j < work.size(); j ++) {
          swork[j](p.first.size() - k - 1, kk) =
            revertProgramInvariant<num_t>(make_pair(
              p.first[k][j * p.first[k].size() / work.size() + kk] /
                p.first[k][p.first[k].size() - 1],
            num_t(int(1)) ));
          swork[j](p.first.size() + work[j].rows() + k, kk) =
            revertProgramInvariant<num_t>(make_pair(
              p.second[k][j * p.first[k].size() / work.size() + kk] /
                p.second[k][p.second[k].size() - 1],
            num_t(int(1)) ));
        }
    if(! savep2or3<num_t>(argv[i], normalize<num_t>(swork)) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

