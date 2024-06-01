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
#include <unistd.h>

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
using std::move;
using std::pair;

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
//#define int int32_t
  assert(1 < argc);
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  vector<vector<SimpleMatrix<num_t> > > in;
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    if(work.size() == 3)
      work = rgb2xyz<num_t>(work);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int ii = 0; ii < work.size(); ii ++)
      for(int j = 0; j < work[ii].rows(); j ++)
        for(int k = 0; k < work[ii].cols(); k ++) {
          work[ii](j, k) += num_t(int(1)) / num_t(int(65536));
          work[ii](j, k) /= num_t(int(1)) + num_t(int(1)) / num_t(int(256));
        }
    in.emplace_back(work);
  }
  in = normalize<num_t>(in);
  for(int i0 = 1; 0 < i0; i0 ++) {
    const int  color0(8 * in.size() / i0);
    const auto p(predMat<num_t, P0J>(in, i0));
    if(! p.first.size() || ! p.second.size()) break;
    for(int i = 0; i < p.first.size(); i ++) {
      const int color(max(int(1), min(int(65535), color0 - i)) );
      if(! savep2or3<num_t>((std::string("predg-forward-") + std::to_string(i) + std::string("-") + std::to_string(i0) + std::string(".ppm")).c_str(), p.first[i].size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(p.first[i])) : p.first[i], color) )
        cerr << "failed to save." << endl;
      if(! savep2or3<num_t>((std::string("predg-backward-") + std::to_string(i) + std::string("-") + std::to_string(i0) + std::string(".ppm")).c_str(), p.second[i].size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(p.second[i])) : p.second[i], color) )
        cerr << "failed to save." << endl;
    }
  }
  cerr << " Done" << endl;
  return 0;
}

