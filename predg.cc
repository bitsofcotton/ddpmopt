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

#define IMG_BITS 8

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
    vector<SimpleMatrix<num_t> > work2;
    work2.resize(work.size() * IMG_BITS);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int ii = 0; ii < work2.size(); ii ++) {
      work2[ii].resize(work[ii / IMG_BITS].rows(), work[ii / IMG_BITS].cols());
      work2[ii].O();
      for(int j = 0; j < work2[ii].rows(); j ++)
        for(int k = 0; k < work2[ii].cols(); k ++) {
          work2[ii](j, k) = work[ii / IMG_BITS](j, k) * num_t(ii % IMG_BITS ? int(1) << (ii % IMG_BITS) : int(1));
          work2[ii](j, k) -= floor(work2[ii](j, k));
        }
    }
    in.emplace_back(work2);
  }
  in = normalize<num_t>(in);
  for(int i = 0; i < in.size(); i ++)
    for(int j = 0; j < in[i].size(); j ++)
      in[i][j] /= num_t(int(2));
  pair<vector<vector<SimpleMatrix<num_t> > >, vector<vector<SimpleMatrix<num_t> > > > pw;
  const int  color(IMG_BITS * in.size());
        auto p(predMat<num_t>(in));
  for(int i = 0; i < p.first.size(); i ++) {
    vector<SimpleMatrix<num_t> > bm;
    vector<SimpleMatrix<num_t> > fm;
    bm.resize(p.first[i].size() / IMG_BITS);
    fm.resize(p.second[i].size() / IMG_BITS);
    for(int j = 0; j < p.first[i].size(); j ++)
      if(! (j % IMG_BITS)) {
        bm[j / IMG_BITS] = p.first[i][j];
        fm[j / IMG_BITS] = p.second[i][j];
      } else {
        bm[j / IMG_BITS] += p.first[i][j]  / num_t(int(1) << (j % IMG_BITS));
        fm[j / IMG_BITS] += p.second[i][j] / num_t(int(1) << (j % IMG_BITS));
      }
    pw.first.emplace_back(move(bm));
    pw.second.emplace_back(move(fm));
  }
  for(int i = 0; i < pw.first.size(); i ++) {
    if(! savep2or3<num_t>((std::string("predg-forward-") + std::to_string(i) + std::string(".ppm")).c_str(), pw.first[i].size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(pw.first[i])) : pw.first[i], color) )
      cerr << "failed to save." << endl;
    if(! savep2or3<num_t>((std::string("predg-backward-") + std::to_string(i) + std::string(".ppm")).c_str(), pw.second[i].size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(pw.second[i])) : pw.second[i], color) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

