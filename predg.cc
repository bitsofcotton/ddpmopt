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
using std::pair;

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
//#define int int32_t
  assert(1 < argc);
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(num_t(2)) << endl;
  vector<vector<SimpleMatrix<num_t> > > in;
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    in.emplace_back(work);
  }
  // N.B. we need mipmap for gaining global contexts.
  vector<vector<vector<SimpleMatrix<num_t> > > > mipmap;
  mipmap.resize(int(log(num_t(int(min(in[0][0].rows(), in[0][0].cols() ))) ) / log(num_t(int(2)) ) ) );
  for(int i = 0; i < mipmap.size(); i ++) {
    if(i) {
      mipmap[i].resize(mipmap[i - 1].size());
      for(int j = 0; j < mipmap[i - 1].size(); j ++)
        for(int k = 0; k < mipmap[i - 1][j].size(); k ++) {
          SimpleMatrix<num_t> work((mipmap[i - 1][j][k].rows() + 1) / 2,
                                   (mipmap[i - 1][j][k].cols() + 1) / 2);
          for(int ii = 0; ii < work.rows(); ii ++)
            for(int jj = 0; jj < work.cols(); jj ++) {
              work(ii, jj) = num_t(int(0));
              int cnt(0);
              for(int iii = ii * 2;
                      iii < min((ii + 1) * 2, int(mipmap[i - 1][j][k].rows()) );
                      iii ++)
                for(int jjj = jj * 2;
                        jjj < min((jj + 1) * 2, int(mipmap[i - 1][j][k].cols()) );
                        jjj ++, cnt ++)
                  work(ii, jj) += mipmap[i - 1][j][k](iii, jjj);
              if(cnt) work(ii, jj) /= num_t(cnt);
            }
          mipmap[i][j].emplace_back(work);
        }
    } else
      mipmap[0] = move(in);
  }
  pair<vector<vector<SimpleMatrix<num_t> > >, vector<vector<SimpleMatrix<num_t> > > > p;
  for(int i = mipmap.size() - 1; 0 <= i; i --) {
    cerr << "mipmap: " << mipmap.size() - i << " / " << mipmap.size() << endl;
    auto work(mipmap[i]);
    if(i < mipmap.size() - 1) {
      for(int j = 0; j < work.size(); j ++)
        for(int k = 0; k < work[j].size(); k ++)
          for(int ii = 0; ii < work[j][k].rows(); ii ++)
            for(int jj = 0; jj < work[j][k].cols(); jj ++)
              work[j][k](ii, jj) -= mipmap[i + 1][j][k](ii / 2, jj / 2);
    }
    auto pk(predMat<num_t>(work));
    if(p.first.size()) {
      for(int i = 0; i < min(p.first.size(), pk.first.size()); i ++)
        for(int j = 0; j < min(p.first[i].size(), pk.first[i].size()); j ++) {
          auto first( p.first[ i][j]);
          auto second(p.second[i][j]);
          p.first[ i][j].resize(pk.first[ i][j].rows(), pk.first[ i][j].cols());
          p.second[i][j].resize(pk.second[i][j].rows(), pk.second[i][j].cols());
          for(int ii = 0; ii < p.first[i][j].rows(); ii ++)
            for(int jj = 0; jj < p.first[i][j].cols(); jj ++) {
              p.first[ i][j](ii, jj) = first( ii / 2, jj / 2);
              p.second[i][j](ii, jj) = second(ii / 2, jj / 2);
            }
          p.first[ i][j] += pk.first[ i][j];
          p.second[i][j] += pk.second[i][j];
        }
    } else p = move(pk);
    for(int i = 0; i < p.first.size(); i ++) {
      if(! savep2or3<num_t>((std::string("predg-forward-") + std::to_string(i) + std::string(".ppm")).c_str(), normalize<num_t>(autoGamma<num_t>(normalize<num_t>(p.first[i]) ))) )
        cerr << "failed to save." << endl;
      if(! savep2or3<num_t>((std::string("predg-backward-") + std::to_string(i) + std::string(".ppm")).c_str(), normalize<num_t>(autoGamma<num_t>(normalize<num_t>(p.second[i]) ))) )
        cerr << "failed to save." << endl;
    }
  }
  cerr << " Done" << endl;
  return 0;
}

