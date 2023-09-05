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
    num_t norm2(int(0));
    pwork.emplace_back(vector<SimpleVector<num_t> >());
    for(int ii = 0; 0 <= ii; ii ++) {
      if(ii) pwork[ii] = pwork[ii - 1];
      else pwork[ii].resize(work[0].rows() * 2 - 1,
             SimpleVector<num_t>(work[0].cols() * work.size()).O());
      int ok_cnt(0);
      if(ii)
        for(int k = 0; k < pwork[ii].size(); k += 2) {
          const auto ga(revertProgramInvariant<num_t>(make_pair(makeProgramInvariant<num_t>(pwork[ii - 1][k]).second, num_t(int(1)) )) );
          if(abs(ga) <= sqrt(SimpleMatrix<num_t>().epsilon()) ) ok_cnt ++
          for(int kk = 0; kk < pwork[ii][k].size(); kk ++)
            pwork[ii][k][kk] -= ga;
          if(0 < k)
            pwork[ii][k - 1] = (pwork[ii][k] + pwork[ii][k - 2]) / num_t(int(2));
        }
      else
        for(int k = 0; k < work[0].rows(); k += 2) {
          for(int j = 0; j < work.size(); j ++)
            pwork[ii][k].setVector(j * work[j].cols(), work[j].row(k / 2));
          if(0 < k)
            pwork[ii][k - 1] = (pwork[ii][k] + pwork[ii][k - 2]) / num_t(int(2));
          norm2 += pwork[ii][k].dot(pwork[ii][k]);
        }
      if(ok_cnt < pwork[ii].size() / 2)
        pwork.emplace_back(vector<SimpleVector<num_t> >());
      else
        break;
    }
    norm2 /= num_t(int(pwork[0].size()) / 2);
    cerr << "needs total " << pwork.size() << " loop" << endl;
    auto p(predv<num_t>(pwork[0]));
    vector<SimpleMatrix<num_t> > swork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + p.first.size() / 2 * 2,
        work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      swork[j].setMatrix(p.first.size() / 2, 0, work[j]);
    for(int ii = 1; ii < pwork.size(); ii ++) {
      cerr << "loop#" << ii << endl;
      auto q(predv<num_t>(pwork[ii]));
      for(int k = 1; k < p.first.size(); k += 2) {
        p.first[k]  += q.first[k];
        p.second[k] += q.second[k];
      }
    }
    for(int k = 1; k < p.first.size(); k += 2) {
      p.first[k]  *= sqrt(norm2 / p.first[k].dot( p.first[k]));
      p.second[k] *= sqrt(norm2 / p.second[k].dot(p.second[k]));
      for(int j = 0; j < work.size(); j ++) {
        swork[j].row(p.first.size() / 2 - k / 2 - 1) =
          p.second[k].subVector(j * p.first[k].size() / work.size(),
            swork[j].cols());
        swork[j].row(p.first.size() / 2 + work[j].rows() + k / 2) =
          p.first[k].subVector(j * p.first[k].size() / work.size(),
            swork[j].cols());
      }
    }
    if(! savep2or3<num_t>(argv[i], normalize<num_t>(swork)) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

