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
    for(int i = 0; i < pwork.size(); i ++) {
      for(int j = 0; j < work.size() * IMG_BITS; j ++) {
        pwork[i].emplace_back(work[j / IMG_BITS].row(i));
        for(int k = 0; k < pwork[i][j].size(); k ++) {
          pwork[i][j][k] *= num_t(j % IMG_BITS ? int(1) << (j % IMG_BITS) : int(1));
          pwork[i][j][k] -= floor(pwork[i][j][k]);
        }
      }
    }
    vector<SimpleMatrix<num_t> > wwork;
    const auto color(max(int(1), min(int(65535), int(8 * work[0].rows()) )) );
          auto p(predVec<num_t>(pwork));
    if(! p.first.size() || ! p.second.size()) break;
    wwork.resize(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + 2, work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      wwork[j].setMatrix(1, 0, work[j]);
    for(int j = 0; j < work.size() * IMG_BITS; j ++) {
      for(int k = 0; k < p.first[j].size(); k ++) {
        p.first[j][k]  = sgn<num_t>(p.first[j][k]);
        p.second[j][k] = sgn<num_t>(p.second[j][k]);
      }
      if(! (j % IMG_BITS)) {
        wwork[j / IMG_BITS].row(0) =
          p.second[j].subVector(0, wwork[j / IMG_BITS].cols());
        wwork[j / IMG_BITS].row(work[j / IMG_BITS].rows() + 1) =
          p.first[j].subVector(0, wwork[j / IMG_BITS].cols());
      } else {
        wwork[j / IMG_BITS].row(0) +=
          p.second[j].subVector(0, wwork[j / IMG_BITS].cols()) /
            num_t(int(1) << (j % IMG_BITS));
        wwork[j / IMG_BITS].row(work[j / IMG_BITS].rows() + 1) +=
          p.first[j].subVector(0, wwork[j / IMG_BITS].cols()) /
            num_t(int(1) << (j % IMG_BITS));
      }
    }
    if(! savep2or3<num_t>(argv[i0], wwork.size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(wwork)) : normalize<num_t>(wwork), color) )
        cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

