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

//#define int int64_t
#define int int32_t
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

#include <stdlib.h>

template <typename T> static inline vector<vector<SimpleMatrix<T> > > shrinken(const vector<vector<SimpleMatrix<T> > >& in, const int& r = 3) {
  auto shrink(in);
  for(int i = 0; i < shrink.size(); i ++)
    for(int j = 0; j < shrink[i].size(); j ++) {
      shrink[i][j] = SimpleMatrix<T>((in[0][0].rows() + r - 1) / r,
                                     (in[0][0].cols() + r - 1) / r).O();
      for(int ii = 0; ii < shrink[i][j].rows(); ii ++)
        for(int jj = 0; jj < shrink[i][j].cols(); jj ++) {
          shrink[i][j](ii, jj) = T(int(0));
          int cnt(0);
          for(int iii = max(0, ii * r - r / 2);
                  iii < max(0, min(in[i][j].rows(),
                    ii * r + r / 2)); iii ++)
            for(int jjj = max(0, jj * r - r / 2);
                    jjj < max(0, min(in[i][j].cols(),
                      jj * r + r / 2)); jjj ++, cnt ++)
              shrink[i][j](ii, jj) += in[i][j](iii, jjj);
          if(cnt) shrink[i][j](ii, jj) /= T(cnt);
        }
    }
  return shrink;
}

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  const auto sz(3);
  const auto m(argv[1][0]);
  if(m == '-') {
    vector<vector<SimpleVector<num_t> > > L;
    L.resize(sz * sz);
    for(int i = 0; i < L.size(); i ++) {
      auto sz0(0);
      std::cin >> sz0;
      L[i].reserve(sz0 * 2);
      for(int j = 0; j < sz0; j ++) {
        SimpleVector<num_t> b;
        std::cin >> b;
        auto bb(b);
        bb[sz * sz] = num_t(int(0));
        L[i].emplace_back(bb /= sqrt(bb.dot(bb)));
        b /= - num_t(b[sz * sz]);
        b[sz * sz] = num_t(int(0));
        L[i].emplace_back(b);
        assert(L[0][0].size() == L[i][j].size());
      }
      assert(L[i].size());
    }
    for(int i = 2; i < argc; i ++) {
      cerr << i - 2 << " / " << argc - 2 << endl;
      vector<SimpleMatrix<num_t> > out;
      if(! loadp2or3<num_t>(out, argv[i])) return - 1;
      auto shrink(out);
      {
        vector<vector<SimpleMatrix<num_t> > > buf;
        buf.emplace_back(out);
        shrink = shrinken(buf, sz)[0];
      }
      auto outs(out);
      for(int n = 0; n < outs.size(); n ++)
        outs[n] = SimpleMatrix<num_t>(out[n].rows() / sz * sz,
                                      out[n].cols() / sz * sz).O();
      SimpleVector<num_t> buf(shrink.size() * shrink[0].rows() * shrink[0].cols() * sz * sz);
      SimpleVector<num_t> v(sz * sz + 1);
      SimpleVector<num_t> vv(v.size() + 1);
      buf.O();
      for(int j = 0; j < shrink.size(); j ++)
        for(int idx = 0; idx < sz * sz; idx ++)
          for(int m = 0; m < shrink[j].rows() * shrink[j].cols(); m ++) {
            int   Midx(- 1);
            num_t M(int(0));
            for(int ii = 0; ii < sz; ii ++)
              for(int jj = 0; jj < sz; jj ++)
                if((m / shrink[j].cols()) + ii < shrink[j].rows() &&
                   (m % shrink[j].cols()) + jj < shrink[j].cols())
                  v[ii * sz + jj] = shrink[j](
                    (m / shrink[j].cols()) + ii,
                    (m % shrink[j].cols()) + jj);
                else goto next0;
            v[sz * sz] = num_t(int(0));
            vv = makeProgramInvariant<num_t>(v).first;
            for(int k = 0; k < L[idx].size(); k += 2) {
              const auto lM(abs(L[idx][k].dot(vv) ));
              if(Midx < 0 || lM < M) {
                M    = lM;
                Midx = k;
              }
            }
            assert(0 <= Midx && Midx < L[idx].size());
            buf[m + idx * shrink[j].rows() * shrink[j].cols() +
                    j * sz * sz * shrink[j].rows() * shrink[j].cols()] =
              L[idx][++ Midx].dot(vv);
           next0:
            ;
          }
      buf = revertProgramInvariant<num_t>(buf);
      SimpleVector<num_t> normal(sz * sz);
      for(int idx = 0; idx < normal.size(); idx ++) {
        int cnt(0);
        for(int j = 0; j < shrink.size(); j ++)
          for(int m = 0; m < shrink[j].rows() * shrink[j].cols(); m ++)
            if((m / shrink[j].cols()) * sz + idx / sz < outs[j].rows() &&
               (m % shrink[j].cols()) * sz + idx % sz < outs[j].cols()) {
              const auto& b(buf[m + idx * shrink[j].rows() * shrink[j].cols() +
                j * sz * sz * shrink[j].rows() * shrink[j].cols()]);
              normal[idx] += b * b;
              cnt ++;
            }
        normal[idx] = sqrt(normal[idx] / num_t(cnt));
      }
      for(int j = 0; j < shrink.size(); j ++)
        for(int idx = 0; idx < sz * sz; idx ++)
          for(int m = 0; m < shrink[j].rows() * shrink[j].cols(); m ++)
            if((m / shrink[j].cols()) * sz + idx / sz < outs[j].rows() &&
               (m % shrink[j].cols()) * sz + idx % sz < outs[j].cols())
              outs[j]((m / shrink[j].cols()) * sz + idx / sz,
                      (m % shrink[j].cols()) * sz + idx % sz) =
                buf[m + idx * shrink[j].rows() * shrink[j].cols() +
                        j * sz * sz * shrink[j].rows() * shrink[j].cols()] /
                  normal[idx];
      if(! savep2or3<num_t>(argv[i], normalize<num_t>(outs)) )
        cerr << "failed to save." << endl;
    }
  } else if(m == '+') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.resize(argc - 2);
    for(int i = 2; i < argc; i ++) {
      if(! loadp2or3<num_t>(in[i - 2], argv[i])) continue;
      assert(in[0][0].rows() == in[i - 2][0].rows() &&
             in[0][0].cols() == in[i - 2][0].cols());
    }
    auto shrink(shrinken<num_t>(in, sz));
    for(int i0 = 0; i0 < sz * sz; i0 ++) {
      cerr << i0 << " / " << sz * sz << endl;
      vector<SimpleVector<num_t> > v;
      v.reserve(in.size() * in[0].size() * in[0][0].rows() * in[0][0].cols());
      SimpleVector<num_t> tv(sz * sz + 1);
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < in[i].size(); j ++)
          for(int m = 0; m < shrink[i][j].rows() * shrink[i][j].cols(); m ++) {
            for(int ii = 0; ii < sz; ii ++)
              for(int jj = 0; jj < sz; jj ++)
                if((m / shrink[i][j].cols()) + ii < shrink[i][j].rows() &&
                   (m % shrink[i][j].cols()) + jj < shrink[i][j].cols())
                  tv[ii * sz + jj] = shrink[i][j](
                    (m / shrink[i][j].cols()) + ii,
                    (m % shrink[i][j].cols()) + jj);
                else goto next1;
            tv[sz * sz] = in[i][j](
               min(in[i][j].rows() - 1,
                 (m / shrink[i][j].cols()) * sz + i0 / sz),
               min(in[i][j].cols() - 1,
                 (m % shrink[i][j].cols()) * sz + i0 % sz) );
            v.emplace_back(tv);
           next1:
            ;
          }
      auto c(crush(v));
      int  cnt(0);
      for(int i = 0; i < c.size(); i ++)
        if(! (c[i].first.size() < v[0].size() + 1 + 1) ) cnt ++;
      cout << cnt << endl;
      for(int i = 0; i < c.size(); i ++) {
        SimpleMatrix<num_t> lc(c[i].first.size(), v[0].size() + 1);
        for(int j = 0; j < lc.rows(); j ++)
          lc.row(j) = makeProgramInvariant<num_t>(c[i].first[j]).first;
        if(lc.rows() <= lc.cols()) {
          /* */
          ;
        } else
          cout << linearInvariant<num_t>(lc);
      }
    }
  }
  return 0;
}

