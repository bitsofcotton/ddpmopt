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
  vector<vector<SimpleVector<num_t> > > in;
  int color(0);
  int rows(0);
  int cols(0);
  int rrows(0);
  int rcols(0);
  in.emplace_back(vector<SimpleVector<num_t> >());
  SimpleMatrix<num_t> resizeL, resizeR;
  for(int ii = 0; 0 <= ii && ii <= argc * 2 - 3; ii ++) {
    if(ii) in[ii] = in[ii - 1];
    else in[ii].reserve(argc * 2 - 3);
    int ok_cnt(0);
    for(int i = 1; i < argc; i ++)
      if(ii) {
        const auto ga(revertProgramInvariant<num_t>(make_pair(makeProgramInvariant<num_t>(in[ii - 1][(i - 1) * 2]).second, num_t(int(1)) )) );
        if(abs(ga) <= num_t(int(1)) / num_t(int(2))) ok_cnt ++;
        else cerr << abs(ga) << endl;
        for(int j = 0; j < in[ii][(i - 1) * 2].size(); j ++)
          in[ii][(i - 1) * 2][j] -= ga;
        if(1 < i) in[ii][(i - 1) * 2 - 1] = (in[ii][(i - 1) * 2] + in[ii][(i - 1) * 2 - 2]) / num_t(int(2));
      } else {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i])) continue;
        assert(! color || color == work.size());
        assert(! rows  || rows  == work[0].rows());
        assert(! cols  || cols  == work[0].cols());
        color = work.size();
        rows  = work[0].rows();
        cols  = work[0].cols();
        if(! rrows) {
          const auto pixels(num_t(argc * 2 - 3) / num_t(color) );
          rrows = int(pixels * num_t(rows) * sqrt(num_t(rows * rows + cols * cols)) / num_t(rows * cols));
          rcols = int(pixels * num_t(cols) * sqrt(num_t(rows * rows + cols * cols)) / num_t(rows * cols));
          assert(rrows && rcols);
          resizeL = (dft<num_t>(- rrows).subMatrix(0, 0, rrows, min(rrows, rows)) * dft<num_t>(rows).subMatrix(0, 0, min(rrows, rows), rows)).template real<num_t>();
          resizeR = (dft<num_t>(- rcols).subMatrix(0, 0, rcols, min(rcols, cols)) * dft<num_t>(cols).subMatrix(0, 0, min(rcols, cols), cols)).template real<num_t>().transpose();
        }
        for(int j = 0; j < work.size(); j ++)
          work[j] = resizeL * work[j] * resizeR;
        if(1 < i) in[ii].emplace_back(SimpleVector<num_t>(work.size() * work[0].rows() * work[0].cols()));
        in[ii].emplace_back(SimpleVector<num_t>(work.size() * work[0].rows() * work[0].cols()));
        for(int j = 0; j < work.size(); j ++)
          for(int k = 0; k < work[j].rows(); k ++)
            in[ii][in[ii].size() - 1].setVector(j * work[j].rows() * work[j].cols() +
                k * work[j].cols(), work[j].row(k));
        if(1 < i) in[ii][in[ii].size() - 2] = (in[ii][in[ii].size() - 3] + in[ii][in[ii].size() - 1]) / num_t(int(2));
      }
    if(ii < argc * 2 - 3 && ok_cnt < argc - 1)
      in.emplace_back(vector<SimpleVector<num_t> >());
    else
      break;
  }
  cerr << "needs total " << in.size() << " loop" << endl;
  vector<SimpleVector<num_t> > out;
  auto p(predv<num_t>(in[0]));
  out = move(p.first);
  if(out.size() & 1) out.erase(out.end() - 1);
  out.insert(out.end(), p.second.begin(), p.second.end());
  for(int ii = 1; ii < in.size(); ii ++) {
    cerr << "loop#" << ii << endl;
    auto q(predv<num_t>(in[ii]));
    for(int i = 1; i < q.first.size(); i += 2) out[i / 2] += q.first[i];
    for(int i = 1; i < q.second.size(); i += 2)
      out[(q.first.size() / 2) * 2 + (i / 2)] += q.second[i];
  }
  vector<SimpleMatrix<num_t> > outs;
  outs.resize(color);
  for(int i = 1; i < out.size(); i += 2) {
    for(int j = 0; j < outs.size(); j ++) {
      outs[j].resize(rrows, rcols);
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

