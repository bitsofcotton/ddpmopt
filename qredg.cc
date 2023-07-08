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
    vector<SimpleVector<num_t> > pwork;
    pwork.resize(work[0].rows() * 2 - 1,
      SimpleVector<num_t>(work[0].cols() * work.size()).O());
    num_t norm2(int(0));
    for(int k = 0; k < work[0].rows(); k += 2) {
      for(int j = 0; j < work.size(); j ++)
        pwork[k].setVector(j * work[j].cols(), work[j].row(k / 2));
      if(0 < k) {
        pwork[k - 1] = (pwork[k] + pwork[k - 2]) / num_t(int(2));
        norm2 += pwork[k - 1].dot(pwork[k - 1]);
      }
      norm2 += pwork[k].dot(pwork[k]);
    }
    norm2 /= num_t(work[0].rows() * 2 - 1);
    auto p(predv<num_t>(pwork));
    vector<SimpleMatrix<num_t> > swork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + p.first.size() / 2 * 2,
        work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      swork[j].setMatrix(p.first.size() / 2, 0, work[j]);
    num_t pnorm2(int(0));
    for(int j = 0; j < p.first.size(); j ++) {
      p.first[j] = autoLevel<num_t>(p.first[j], int(num_t(p.first[j].size()) / num_t(int(3)) ));
      pnorm2 += p.first[j].dot(p.first[j]);
    }
    for(int j = 0; j < p.second.size(); j ++) {
      p.second[j] = autoLevel<num_t>(p.second[j], int(num_t(p.second[j].size()) / num_t(int(3)) ));
      pnorm2 += p.second[j].dot(p.first[j]);
    }
    pnorm2 /= num_t(p.first.size() + p.second.size());
    norm2  /= pnorm2;
    norm2   = sqrt(norm2);
    for(int k = 1; k < p.first.size(); k += 2)
      for(int j = 0; j < work.size(); j ++) {
        swork[j].row(p.first.size() - k / 2 - 1) =
          p.second[k].subVector(j * p.first[k].size() / work.size(),
            swork[j].cols()) * norm2;
        swork[j].row(p.first.size() + work[j].rows() + k / 2) =
          p.first[k].subVector(j * p.first[k].size() / work.size(),
            swork[j].cols()) * norm2;
      }
    if(! savep2or3<num_t>(argv[i], normalize<num_t>(swork)) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

