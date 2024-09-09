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
    auto p(predVec<num_t, _PNOISE_>(pwork));
    vector<SimpleMatrix<num_t> > wwork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + 1, work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      wwork[j].setMatrix(0, 0, work[j]);
    for(int j = 0; j < p.size(); j ++)
      wwork[j].row(wwork[j].rows() - 1) = move(p[j]);
    num_t norm2(int(0));
    for(int j = 0; j < wwork.size(); j ++)
      for(int k = 0; k < wwork[j].rows() - 1; k ++)
        norm2 += wwork[j].row(k).dot(wwork[j].row(k));
    norm2 = norm2 / num_t(wwork[0].rows() - 1);
    num_t norm2_m1(int(0));
    for(int j = 0; j < wwork.size(); j ++)
      norm2_m1 += wwork[j].row(wwork[j].rows() - 1).dot(
        wwork[j].row(wwork[j].rows() - 1) );
    for(int j = 0; j < wwork.size(); j ++)
      wwork[j].row(wwork[j].rows() - 1) *=
        sqrt(norm2 / norm2_m1) * num_t(int(3)) / num_t(int(2));
    if(! savep2or3<num_t>(argv[i0], normalize<num_t>(wwork.size() == 3 ?
      xyz2rgb<num_t>(wwork) : wwork), work[0].rows()) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

