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

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
//#define int int32_t
  assert(1 < argc);
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  for(int i0 = 1; i0 < argc; i0 ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i0])) continue;
    work = normalize<num_t>(work.size() == 3 ? rgb2xyz<num_t>(work) : work);
    vector<vector<SimpleVector<num_t> > > pwork;
    pwork.resize(work[0].rows());
    for(int i = 0; i < pwork.size(); i ++) {
      pwork[i].reserve(work.size());
      for(int j = 0; j < work.size(); j ++)
        pwork[i].emplace_back(work[j].row(i));
    }
    auto p(predVec<num_t>(pwork));
    vector<SimpleMatrix<num_t> > wwork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + p.size(), work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      wwork[j].setMatrix(0, 0, work[j]);
    for(int i = 0; i < p.size(); i ++)
      for(int j = 0; j < p[i].size(); j ++)
        wwork[j].row(i - p.size() + wwork[j].rows()) = move(p[i][j]);
    if(! savep2or3<num_t>(argv[i0], normalize<num_t>(wwork.size() == 3 ?
      xyz2rgb<num_t>(wwork) : wwork) ) )
        cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

