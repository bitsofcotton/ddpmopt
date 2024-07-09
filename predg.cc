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
          work2[ii](j, k) *= num_t(int(2));
          work2[ii](j, k) -= floor(work2[ii](j, k));
          work2[ii](j, k) += num_t(int(1));
          work2[ii](j, k) /= num_t(int(2));
        }
    }
    in.emplace_back(move(work2));
  }
  in = normalize<num_t>(in);
  const auto p(predMat<num_t>(in));
  pair<vector<SimpleMatrix<num_t> >, vector<SimpleMatrix<num_t> > > pw;
  pw.first.resize( p.first.size()  / IMG_BITS);
  pw.second.resize(p.second.size() / IMG_BITS);
  for(int j = 0; j < p.first.size(); j ++)
    if(! ( j % IMG_BITS)) {
      pw.first[ j / IMG_BITS] = move(p.first[ j]);
      pw.second[j / IMG_BITS] = move(p.second[j]);
    } else {
      pw.first[ j / IMG_BITS] += p.first[ j] / num_t(int(1) << (j % IMG_BITS));
      pw.second[j / IMG_BITS] += p.second[j] / num_t(int(1) << (j % IMG_BITS));
    }
  if(! savep2or3<num_t>("predg-f.ppm",
    normalize<num_t>(p.first.size() == 3 ? xyz2rgb<num_t>(p.first)
      : p.first), in.size() * IMG_BITS) )
    cerr << "failed to save." << endl;
  if(! savep2or3<num_t>("predg-b.ppm",
    normalize<num_t>(p.second.size() == 3 ? xyz2rgb<num_t>(p.second)
      : p.second), in.size() * IMG_BITS) )
    cerr << "failed to save." << endl;
  cerr << " Done" << endl;
  return 0;
}

