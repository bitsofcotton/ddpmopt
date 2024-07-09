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

#define IMG_BITS 8

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
      for(int j = 0; j < work.size() * IMG_BITS; j ++) {
        pwork[i].emplace_back(work[j / IMG_BITS].row(i));
        for(int k = 0; k < pwork[i][j].size(); k ++) {
          pwork[i][j][k] *= num_t(j % IMG_BITS ? int(1) << (j % IMG_BITS) : int(1));
          pwork[i][j][k] -= floor(pwork[i][j][k]);
          pwork[i][j][k] *= num_t(int(2));
          pwork[i][j][k] -= floor(pwork[i][j][k]);
          pwork[i][j][k] += num_t(int(1));
          pwork[i][j][k] /= num_t(int(2));
        }
      }
    auto p(predVec<num_t>(pwork));
    vector<SimpleMatrix<num_t> > wwork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + 2, work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      wwork[j].setMatrix(1, 0, work[j]);
    for(int j = 0; j < p.first.size(); j ++)
      if(! (j % IMG_BITS)) {
        wwork[j / IMG_BITS].row(0) = move(p.second[j]);
        wwork[j / IMG_BITS].row(wwork[j / IMG_BITS].rows() - 1) =
          move(p.first[j]);
      } else {
        wwork[j / IMG_BITS].row(0) +=
          p.second[j] / num_t(int(1) << (j % IMG_BITS));
        wwork[j / IMG_BITS].row(wwork[j / IMG_BITS].rows() - 1) +=
          p.first[ j] / num_t(int(1) << (j % IMG_BITS));
      }
    num_t norm2(int(0));
    for(int j = 0; j < wwork.size(); j ++)
      for(int k = 1; k < wwork[j].rows() - 1; k ++)
        norm2 += wwork[j].row(k).dot(wwork[j].row(k));
    norm2 = norm2 / num_t(wwork.size() * (wwork[0].rows() - 2));
    for(int j = 0; j < wwork.size(); j ++) {
      wwork[j].row(0) *= sqrt(norm2 / wwork[j].row(0).dot(wwork[j].row(0)));
      wwork[j].row(wwork[j].rows() - 1) *=
        sqrt(norm2 / wwork[j].row(wwork[j].rows() - 1).dot(
          wwork[j].row(wwork[j].rows() - 1) ) );
    }
    if(! savep2or3<num_t>(argv[i0], normalize<num_t>(wwork.size() == 3 ?
      xyz2rgb<num_t>(wwork) : wwork), work[0].rows() * IMG_BITS) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

