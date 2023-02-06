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

#if defined(_OPENMP)
#include <omp.h>
#endif

//#define int int64_t
#define int int32_t
#include "lieonn.hh"
typedef myfloat num_t;
#include "util.hh"
#include "p0.hh"
#include "p1.hh"
#include "catg.hh"
#include "p2.hh"

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
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  vector<vector<SimpleVector<num_t> > > in;
  int rows(0);
  int cols(0);
  in.reserve(argc - 1);
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    in.emplace_back(vector<SimpleVector<num_t> >());
    in[i - 1].resize(work.size());
    for(int j = 0; j < work.size(); j ++) {
      in[i - 1][j].resize(work[j].rows() * work[j].cols());
      for(int k = 0; k < work[j].rows(); k ++)
        in[i - 1][j].setVector(k * work[j].cols(), work[j].row(k));
      rows = work[j].rows();
      cols = work[j].cols();
      assert(in[0][0].size() == in[i - 1][j].size());
    }
    assert(in[0].size() == work.size());
  }
  vector<vector<SimpleVector<num_t> > > out;
  out.resize(in[0].size());
  for(int i = 0; i < out.size(); i ++) {
    vector<SimpleVector<num_t> > d;
    d.reserve(in.size());
    for(int k = 0; k < in.size(); k ++)
      d.emplace_back(move(in[k][i]));
    auto p(predv<num_t>(d));
    out[i] = move(p.first);
    out[i].insert(out[i].end(), p.second.begin(), p.second.end());
  }
  vector<SimpleMatrix<num_t> > outs;
  outs.resize(out.size());
  for(int i = 0; i < out[0].size(); i ++) {
    for(int j = 0; j < out.size(); j ++) {
      outs[j].resize(rows, cols);
      for(int k = 0; k < outs[j].rows() * outs[j].cols(); k ++)
        outs[j](k / outs[j].cols(), k % outs[j].cols()) =
          revertProgramInvariant<num_t>(make_pair(
            out[j][i][k] / out[j][i][out[j][i].size() - 1], num_t(int(1)) ));
    }
    if(! savep2or3<num_t>((std::string("predg-") + std::to_string(i) + std::string(".ppm")).c_str(), normalize<num_t>(outs)) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

