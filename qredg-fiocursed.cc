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

const char* input[] = {
#include "input.ppm.string.array.txt"
};

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
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  for(int i0 = 0; i0 < sizeof(input) / sizeof(char*); i0 ++) {
    vector<SimpleMatrix<num_t> > work;
    istringstream line(input[i0]);
    if(! loadp2or3<num_t>(work, static_cast<istream&>(line))) continue;
    if(work.size() == 3)
      work = normalize<num_t>(rgb2xyz<num_t>(work));
    vector<vector<SimpleVector<num_t> > > pwork;
    pwork.resize(work[0].rows());
    for(int i = 0; i < pwork.size(); i ++)
      for(int j = 0; j < work.size(); j ++)
        pwork[i].emplace_back(work[j].row(i));
    const auto color(work[0].rows());
          auto p(predVec<num_t>(pwork));
    if(! p.first.size () || ! p.second.size()) continue;
    vector<SimpleMatrix<num_t> > wwork(work.size(),
      SimpleMatrix<num_t>(work[0].rows() + p.first.size() + p.second.size(),
        work[0].cols()).O());
    for(int j = 0; j < work.size(); j ++)
      wwork[j].setMatrix(p.first.size(), 0, work[j]);
    for(int k = 0; k < p.first.size(); k ++)
      for(int j = 0; j < work.size(); j ++) {
        wwork[j].row(p.second.size() - k - 1) = move(p.second[k][j]);
        wwork[j].row(p.first.size() + work[j].rows() + k) = move(p.first[k][j]);
      }
    if(! savep2or3<num_t>((std::string("qredg-") + std::to_string(i0) + std::string(".ppm")).c_str(), wwork.size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(wwork)) : wwork, color) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

