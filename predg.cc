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
  assert(2 < argc);
  vector<vector<SimpleVector<num_t> > > in;
  int color(0);
  int rows(0);
  int cols(0);
  in.resize(std::atoi(argv[1]));
  for(int ii = 0; ii < in.size(); ii ++) {
    if(ii) in[ii] = in[ii - 1];
    else in[ii].reserve(argc * 2 - 3);
    for(int i = 2; i < argc; i ++) {
      if(ii) {
        num_t mm(in[ii][(i - 2) * 2][0]);
        auto  MM(mm);
        for(int j = 0; j < in[ii][(i - 2) * 2].size(); j ++) {
          mm = min(mm, in[ii][(i - 2) * 2][j]);
          MM = max(MM, in[ii][(i - 2) * 2][j]);
        }
        num_t ga(int(0));
        for(int j = 0; j < in[ii][(i - 2) * 2].size(); j ++)
          ga += log(in[ii][(i - 2) * 2][j] - mm + (MM - mm) / num_t(int(2)) );
        ga = exp(ga / num_t(int(in[ii][(i - 2) * 2].size())) );
        for(int j = 0; j < in[ii][(i - 2) * 2].size(); j ++)
          in[ii][(i - 2) * 2][j] = (in[ii][(i - 2) * 2][j] - mm + (MM - mm) / num_t(int(2)) ) - ga + (mm - (MM - mm) / num_t(int(2)) );
        if(2 < i) in[ii][(i - 2) * 2 - 1] = (in[ii][(i - 2) * 2] + in[ii][(i - 2) * 2 - 2]) / num_t(int(2));
      } else {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i])) continue;
        assert(! color || color == work.size());
        assert(! rows  || rows  == work[0].rows());
        assert(! cols  || cols  == work[0].cols());
        color = work.size();
        rows  = work[0].rows();
        cols  = work[0].cols();
        if(2 < i) in[ii].emplace_back(SimpleVector<num_t>(work.size() * work[0].rows() * work[0].cols()));
        in[ii].emplace_back(SimpleVector<num_t>(work.size() * work[0].rows() * work[0].cols()));
        for(int j = 0; j < work.size(); j ++)
          for(int k = 0; k < work[j].rows(); k ++)
            in[ii][in[ii].size() - 1].setVector(j * work[j].rows() * work[j].cols() +
                k * work[j].cols(), work[j].row(k));
        if(2 < i) in[ii][in[ii].size() - 2] = (in[ii][in[ii].size() - 3] + in[ii][in[ii].size() - 1]) / num_t(int(2));
      }
    }
  }
  vector<SimpleVector<num_t> > out;
  auto p(predv<num_t>(in[0]));
  out = move(p.first);
  if(out.size() & 1) out.erase(out.end() - 1);
  out.insert(out.end(), p.second.begin(), p.second.end());
  for(int ii = 1; ii < in.size(); ii ++) {
    auto q(predv<num_t>(in[ii]));
    for(int i = 1; i < q.first.size(); i += 2) out[i / 2] += q.first[i];
    for(int i = 1; i < q.second.size(); i += 2)
      out[(q.first.size() / 2) * 2 + (i / 2)] += q.second[i];
  }
  vector<SimpleMatrix<num_t> > outs;
  outs.resize(color);
  for(int i = 1; i < out.size(); i += 2) {
    for(int j = 0; j < outs.size(); j ++) {
      outs[j].resize(rows, cols);
      outs[j].O();
    }
    for(int j = 0; j < outs.size(); j ++)
      for(int k = 0; k < outs[j].rows(); k ++)
        outs[j].row(k) = out[i].subVector(j * outs[j].rows() * outs[j].cols() +
          k * outs[j].cols(), outs[j].cols());
    if(! savep2or3<num_t>((std::string("predg-") + std::to_string(i / 2) + std::string(".ppm")).c_str(), normalize<num_t>(outs) ))
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

