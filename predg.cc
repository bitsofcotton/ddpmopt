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
using std::move;

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
//#define int int32_t
  assert(1 < argc);
  vector<SimpleVector<num_t> > in;
  int color(0);
  int rows(0);
  int cols(0);
  in.reserve(argc - 1);
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    assert(! color || color == work.size());
    assert(! rows  || rows  == work[0].rows());
    assert(! cols  || cols  == work[0].cols());
    color = work.size();
    rows  = work[0].rows();
    cols  = work[0].cols();
    in.emplace_back(SimpleVector<num_t>(work.size() * work[0].rows() * work[0].cols()));
    for(int j = 0; j < work.size(); j ++)
      for(int k = 0; k < work[j].rows(); k ++)
        in[i - 1].setVector(j * work[j].rows() * work[j].cols() +
          k * work[j].cols(), work[j].row(k));
  }
  vector<SimpleVector<num_t> > out;
  auto p(predv<num_t>(in));
  out = move(p.first);
  out.insert(out.end(), p.second.begin(), p.second.end());
  vector<SimpleMatrix<num_t> > outs;
  outs.resize(color);
  for(int i = 0; i < out.size(); i ++) {
    for(int j = 0; j < outs.size(); j ++) {
      outs[j].resize(rows, cols);
      outs[j].O();
    }
    for(int j = 0; j < outs.size(); j ++)
      for(int k = 0; k < outs[j].rows() * outs[j].cols(); k ++)
        outs[j](k / outs[j].cols(), k % outs[j].cols()) =
          out[i][j * outs[j].rows() * outs[j].cols() + k];
    if(! savep2or3<num_t>((std::string("predg-") + std::to_string(i) + std::string(".ppm")).c_str(), normalize<num_t>(outs)) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

