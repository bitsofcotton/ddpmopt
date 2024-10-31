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

#define ACCURACY 8

using std::cout;
using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::to_string;

template <typename T> vector<vector<SimpleMatrix<T> > > convert(const vector<vector<SimpleMatrix<T> > >& in) {
  vector<vector<SimpleMatrix<num_t> > > inL, inD, inR;
  inL.resize(in.size());
  inD.resize(in.size());
  inR.resize(in.size());
  cerr << "SVD";
  for(int i = 0; i < in.size(); i ++) {
    cerr << "." << flush;
    const auto& work(in[i]);
    inL[i].resize(work.size());
    inD[i].resize(work.size());
    inR[i].resize(work.size());
    for(int j = 0; j < work.size(); j ++) {
      inL[i][j] = work[j].SVD();
      inR[i][j] = work[j].transpose().SVD().transpose();
      inD[i][j] = inL[i][j] * work[j] * inR[i][j];
    }
  }
  inL = normalize<num_t>(inL);
  inD = normalize<num_t>(inD);
  inR = normalize<num_t>(inR);
  assert(inL.size() == inD.size() && inD.size() == inR.size());
  vector<vector<SimpleMatrix<num_t> > > ins;
  ins.resize(inL.size());
  for(int i = 0; i < ins.size(); i ++) {
    ins[i].resize(inL[i].size() * ACCURACY * 3);
    for(int j = 0; j < inL[i].size() * ACCURACY; j ++) {
      ins[i][j].resize(inL[i][0].rows(), inL[i][0].cols());
      ins[i][j].O();
      for(int ii = 0; ii < inL[i][0].rows(); ii ++)
        for(int jj = 0; jj < inL[i][0].cols(); jj ++) {
          inL[i][j / ACCURACY](ii, jj) *= num_t(int(2));
          if(abs(inL[i][j / ACCURACY](ii, jj)) >= num_t(int(1)))
            ins[i][j](ii, jj) = sgn<num_t>(inL[i][j / ACCURACY](ii, jj));
          inL[i][j / ACCURACY](ii, jj) -=
            floor(inL[i][j / ACCURACY](ii, jj));
        }
    }
    for(int j = 0; j < inD[i].size() * ACCURACY; j ++) {
      ins[i][j + inL[i].size() * ACCURACY].resize(inD[i][0].rows(), inD[i][0].cols());
      ins[i][j + inL[i].size() * ACCURACY].O();
      for(int ii = 0; ii < inD[i][0].rows(); ii ++)
        for(int jj = 0; jj < inD[i][0].cols(); jj ++) {
          inD[i][j / ACCURACY](ii, jj) *= num_t(int(2));
          if(abs(inD[i][j / ACCURACY](ii, jj)) >= num_t(int(1)))
            ins[i][j + inL[i].size() * ACCURACY](ii, jj) =
              sgn<num_t>(inD[i][j / ACCURACY](ii, jj));
          inD[i][j / ACCURACY](ii, jj) -=
            floor(inD[i][j / ACCURACY](ii, jj));
        }
    }
    for(int j = 0; j < inR[i].size() * ACCURACY; j ++) {
      ins[i][j + inL[i].size() * ACCURACY * 2].resize(inR[i][0].rows(), inR[i][0].cols());
      ins[i][j + inL[i].size() * ACCURACY * 2].O();
      for(int ii = 0; ii < inR[i][0].rows(); ii ++)
        for(int jj = 0; jj < inR[i][0].cols(); jj ++) {
          inR[i][j / ACCURACY](ii, jj) *= num_t(int(2));
          if(abs(inR[i][j / ACCURACY](ii, jj)) >= num_t(int(1)))
            ins[i][j + inL[i].size() * ACCURACY * 2](ii, jj) =
              sgn<num_t>(inR[i][j / ACCURACY](ii, jj));
          inR[i][j / ACCURACY](ii, jj) -=
            floor(inR[i][j / ACCURACY](ii, jj));
        }
    }
  }
  return ins;
}

template <typename T> vector<SimpleMatrix<T> > revert(const vector<SimpleMatrix<T> >& in) {
  vector<SimpleMatrix<T> > res;
  res.resize(in.size() / ACCURACY / 3);
  for(int i = 0; i < res.size(); i ++) {
    auto L(in[i * ACCURACY * 3]);
    auto D(in[i * ACCURACY * 3 + ACCURACY]);
    auto R(in[i * ACCURACY * 3 + ACCURACY * 2]);
    for(int ii = 0; ii < L.rows(); ii ++)
      for(int jj = 0; jj < L.cols(); jj ++)
        L(ii, jj) -= num_t(int(1)) / num_t(int(2));
    for(int ii = 0; ii < D.rows(); ii ++)
      for(int jj = 0; jj < D.cols(); jj ++)
        D(ii, jj) -= num_t(int(1)) / num_t(int(2));
    for(int ii = 0; ii < R.rows(); ii ++)
      for(int jj = 0; jj < R.cols(); jj ++)
        R(ii, jj) -= num_t(int(1)) / num_t(int(2));
    for(int j = 1; j < ACCURACY; j ++) {
      for(int ii = 0; ii < L.rows(); ii ++)
        for(int jj = 0; jj < L.cols(); jj ++)
          L(ii, jj) += (in[i * ACCURACY * 3 + j](ii, jj) - num_t(int(1)) / num_t(int(2)) ) * pow(num_t(int(2)), - num_t(int(j)) );
      for(int ii = 0; ii < D.rows(); ii ++)
        for(int jj = 0; jj < D.cols(); jj ++)
          D(ii, jj) += (in[i * ACCURACY * 3 + ACCURACY + j](ii, jj) - num_t(int(1)) / num_t(int(2)) ) * pow(num_t(int(2)), - num_t(int(j)) );
      for(int ii = 0; ii < R.rows(); ii ++)
        for(int jj = 0; jj < R.cols(); jj ++)
          R(ii, jj) += (in[i * ACCURACY * 3 + ACCURACY * 2 + j](ii, jj) - num_t(int(1)) / num_t(int(2)) ) * pow(num_t(int(2)), - num_t(int(j)) );
    }
    res[i] = L * D * R;
  }
  return res;
}

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
    in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
  }
  in = normalize<num_t>(convert<num_t>(in));
  {
    auto ins(in);
    auto p(predMat<num_t>(ins, - 1, 1));
    assert(p.size() == 1);
    p[0] = revert<num_t>(p[0]);
    if(! savep2or3<num_t>("predg.ppm",
        normalize<num_t>(p[0].size() == 3 ? xyz2rgb<num_t>(p[0]) : p[0])) )
          cerr << "failed to save." << endl;
  }
  {
    auto p(predMat<num_t>(in));
    for(int i = 0; i < p.size(); i ++) {
      p[i] = revert<num_t>(p[i]);
      if(! savep2or3<num_t>((string("predg-") + to_string(i) +
        string(".ppm")).c_str(),
          normalize<num_t>(p[i].size() == 3 ? xyz2rgb<num_t>(p[i]) : p[i])) )
            cerr << "failed to save." << endl;
    }
  }
  cerr << " Done" << endl;
  return 0;
}

