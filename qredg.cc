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
    SimpleVector<int> cwork;
    int color(65535);
    for(int j0 = 1; 0 < j0; j0 ++) {
      const auto p(predVec<num_t>(pwork, j0));
      color = min(color, max(int(1), min(int(65535), 8 * work[0].rows() / j0
        - int(p.first.size()) )) );
      if(! p.first.size () || ! p.second.size()) break;
      vector<SimpleMatrix<num_t> > swork(work.size(),
        SimpleMatrix<num_t>(work[0].rows() + p.first.size() + p.second.size(),
          work[0].cols()).O());
      for(int j = 0; j < work.size(); j ++)
        swork[j].setMatrix(p.first.size(), 0, work[j]);
      if(j0 == 1) {
        wwork = vector<SimpleMatrix<num_t> >(swork);
        cwork.resize(wwork[0].rows());
        cwork.O(0);
      }
      for(int k = 0; k < p.first.size(); k ++)
        for(int j = 0; j < work.size() * IMG_BITS; j ++)
          if(! (j % IMG_BITS)) {
            swork[j / IMG_BITS].row(p.second.size() - k - 1) =
              p.second[k][j].subVector(0, swork[j / IMG_BITS].cols());
            swork[j / IMG_BITS].row(p.first.size() +
                work[j / IMG_BITS].rows() + k) =
              p.first[k][j].subVector(0, swork[j / IMG_BITS].cols());
          } else {
            swork[j / IMG_BITS].row(p.second.size() - k - 1) +=
              p.second[k][j].subVector(0, swork[j / IMG_BITS].cols()) /
              num_t(int(1) << (j % IMG_BITS));
            swork[j / IMG_BITS].row(p.first.size() +
                work[j / IMG_BITS].rows() + k) +=
              p.first[k][j].subVector(0, swork[j / IMG_BITS].cols()) /
              num_t(int(1) << (j % IMG_BITS));
          }
      const auto ustart((cwork.size() - swork[0].rows()) / 2);
      for(int j = ustart; j < cwork.size() - ustart; j ++) {
        for(int m = 0; m < wwork.size(); m ++)
          wwork[m].row(j) += swork[m].row(j - ustart);
        cwork[j] ++;
      }
    }
    for(int j = 0; j < cwork.size(); j ++)
      for(int m = 0; m < wwork.size(); m ++)
        wwork[m].row(j) /= num_t(cwork[j]);
    if(! savep2or3<num_t>(argv[i0], wwork.size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(wwork)) : wwork, color) )
        cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

