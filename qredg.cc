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
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(num_t(2)) << endl;
  for(int i0 = 1; i0 < argc; i0 ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i0])) continue;
    vector<vector<SimpleVector<num_t> > > pwork;
    pwork.resize(work[0].rows());
    for(int i = 0; i < pwork.size(); i ++)
      for(int j = 0; j < work.size(); j ++)
        pwork[i].emplace_back(work[j].row(i));
    vector<vector<vector<SimpleVector<num_t> > > > mipwork;
    mipwork.resize(int(log(num_t(int(pwork[0].size() ))) / log(num_t(int(2)) ) ) );
    for(int i = 0; i < mipwork.size(); i ++) {
      if(i) {
        mipwork[i].reserve(mipwork[i - 1].size());
        for(int j = 0; j < mipwork[i - 1].size(); j ++)
          for(int k = 0; k < mipwork[i - 1][j].size(); k ++) {
            SimpleVector<num_t> work((mipwork[i - 1][j][k].size() + 1) / 2);
            for(int ii = 0; ii < work.size(); ii ++) {
              work[ii] = num_t(int(0));
              int cnt(0);
              for(int iii = ii * 2;
                      iii < min((ii + 1) * 2, int(mipwork[i - 1][j][k].size()) );
                      iii ++, cnt ++) work[ii] += mipwork[i - 1][j][k][iii];
              if(cnt) work[ii] /= num_t(cnt);
            }
            mipwork[i][j].emplace_back(work);
          }
      } else mipwork[0] = move(pwork);
    }
    pair<vector<vector<SimpleVector<num_t> > >, vector<vector<SimpleVector<num_t> > > > p;
    for(int i = mipwork.size() - 1; 0 <= i; i --) {
      std::cerr << "mipwork: " << mipwork.size() - i << " / " << mipwork.size() << std::endl;
      auto mwork(mipwork[i]);
      if(i < mipwork.size() - 1) {
        for(int j = 0; j < mwork.size(); j ++)
          for(int k = 0; k < mwork[j].size(); k ++)
            for(int ii = 0; ii < mwork[j][k].size(); ii ++)
              mwork[j][k][ii] -= mipwork[i + 1][j][k][ii / 2];
      }
      auto pk(predVec<num_t>(mwork));
      if(p.first.size()) {
        for(int i = 0; i < min(p.first.size(), pk.first.size()); i ++)
          for(int j = 0; j < min(p.first[i].size(), pk.first[i].size()); j ++) {
            auto first( p.first[ i][j]);
            auto second(p.second[i][j]);
            p.first[ i][j].resize(pk.first[ i][j].size());
            p.second[i][j].resize(pk.second[i][j].size());
            for(int ii = 0; ii < p.first[i][j].size(); ii ++) {
              p.first[ i][j][ii] = first[ii / 2];
              p.second[i][j][ii] = second[ii / 2];
            }
            p.first[ i][j] += pk.first[ i][j];
            p.second[i][j] += pk.second[i][j];
          }
      } else p = move(pk);
      vector<SimpleMatrix<num_t> > swork(work.size(),
        SimpleMatrix<num_t>(work[0].rows() + p.first.size() + p.second.size(),
          work[0].cols()).O());
      for(int j = 0; j < work.size(); j ++)
        swork[j].setMatrix(p.first.size(), 0, work[j]);
      for(int k = 0; k < p.first.size(); k ++) {
        for(int j = 0; j < work.size(); j ++) {
          swork[j].row(p.second.size() - k - 1) =
            p.second[k][j].subVector(0, swork[j].cols());
          swork[j].row(p.first.size()  + work[j].rows() + k) =
            p.first[k][j].subVector(0, swork[j].cols());
        }
      }
      if(! savep2or3<num_t>(argv[i0], normalize<num_t>(swork)) )
        cerr << "failed to save." << endl;
    }
  }
  cerr << " Done" << endl;
  return 0;
}

