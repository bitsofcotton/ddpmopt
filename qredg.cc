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
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  for(int i0 = 1; i0 < argc; i0 ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i0])) continue;
    if(work.size() == 3)
      work = normalize<num_t>(rgb2xyz<num_t>(work));
    vector<vector<SimpleVector<num_t> > > pwork;
    pwork.resize(work[0].rows());
    for(int i = 0; i < pwork.size(); i ++)
      for(int j = 0; j < work.size(); j ++)
        pwork[i].emplace_back(work[j].row(i));
    auto p(predVec<num_t, _PERSISTENT_>(pwork));
    vector<SimpleMatrix<num_t> > wwork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + p.first.size() + p.second.size(),
        work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      wwork[j].setMatrix(p.first.size(), 0, work[j]);
    for(int j = 0; j < work.size(); j ++)
      for(int k = 0; k < p.first.size(); k ++) {
        wwork[j].row(p.first.size() - k - 1) = move(p.second[k][j]);
        wwork[j].row(work[j].rows() + p.first.size() + k) = move(p.first[k][j]);
      }
    if(! savep2or3<num_t>(argv[i0], wwork.size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(wwork)) : wwork, work[0].rows()) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

