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
#define int int64_t
//#define int int32_t
  assert(1 < argc);
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    vector<vector<SimpleVector<num_t> > > pwork;
    pwork.resize(work[0].rows());
    for(int i = 0; i < pwork.size(); i ++)
      for(int j = 0; j < work.size(); j ++)
        pwork[i].emplace_back(work[j].row(i));
    auto p(predvResizeVec<num_t>(pwork));
    vector<SimpleMatrix<num_t> > swork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + p.first.size() + p.second.size(),
        work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      swork[j].setMatrix(p.first.size(), 0, work[j]);
    for(int k = 0; k < p.first.size(); k ++) {
      for(int j = 0; j < work.size(); j ++) {
        swork[j].row(p.second.size() - k - 1) =
          p.second[k][j].subVector(0, swork[j].cols());
        swork[j].row(p.first.size()  + work[j].rows() + k) =
          p.first[k][j].subVector(0, swork[j].cols());
      }
    }
    if(! savep2or3<num_t>(argv[i], normalize<num_t>(swork)) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

