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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  const auto sz(2);
  const auto m(argv[1][0]);
  if(m == '-') {
    vector<SimpleVector<num_t> > L;
    std::string s;
    while(std::getline(std::cin, s, '\n')) {
      SimpleVector<num_t> l;
      std::stringstream ins(s);
      ins >> l;
      if(l.size() != sz * sz + 1) break;
      L.emplace_back(l / sqrt(l.dot(l)));
      l /= - num_t(l[3]);
      l[3] = num_t(int(0));
      L.emplace_back(l);
    }
    assert(L.size() && ! (L.size() & 1));
    for(int i0 = 2; i0 < argc; i0 ++) {
      cerr << i0 - 2 << " / " << argc - 2 << endl;
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i0])) return - 1;
      if(in.size() != 3) {
        std::cerr << argv[i0] << " doesn't include 3 colors" << std::endl;
        continue;
      }
      vector<SimpleMatrix<num_t> > out;
      out.emplace_back(in[0]);
      out[0].O();
      for(int i = 0; i < in[0].rows(); i ++)
        for(int j = 0; j < in[0].cols(); j ++) {
          SimpleVector<num_t> work(4);
          for(int m = 0; m < in.size(); m ++)
            work[m] = in[m](i, j);
          work[3] = num_t(int(1)) / num_t(int(2));
          auto work2(makeProgramInvariant<num_t>(work, - num_t(int(1)), true).first);
          assert(work2.size() == L[0].size());
          int idx(0);
          for(int m = 2; m < L.size(); m += 2) {
            assert(L[m].size()   == work2.size());
            assert(L[idx].size() == work2.size());
            if(abs(L[idx].dot(work2)) <= abs(L[m].dot(work2)))
              idx = m;
          }
          auto last(sqrt(work.dot(work)));
          for(int ii = 0;
                  ii < 2 * int(- log(SimpleMatrix<num_t>().epsilon()) / log(num_t(int(2))) )
                  && sqrt(work.dot(work) * SimpleMatrix<num_t>().epsilon()) <
                       abs(work[3] - last); ii ++) {
            last = work[3];
            const auto work2(makeProgramInvariant<num_t>(work, - num_t(int(1)), true));
            work[3] = revertProgramInvariant<num_t>(make_pair(L[idx + 1].dot(work2.first) * sgn<num_t>(L[idx].dot(work2.first)), work2.second), true);
          }
          out[0](i, j) = work[3];
        }
      if(! savep2or3<num_t>((std::string(argv[i0]) + std::string(".pgm")).c_str(), out) )
        cerr << "failed to save." << endl;
    }
  } else if(m == '+') {
    vector<vector<SimpleMatrix<num_t> > > in;
    vector<SimpleMatrix<num_t> > out;
    assert(! ((argc - 2) & 1));
    in.resize((argc - 2) / 2);
    out.resize((argc - 2) / 2);
    int cnt(0);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      if(! (i & 1)) {
        assert(work.size() == 1);
        out[i / 2 - 1] = move(work[0]);
        cnt += out[i / 2 - 1].rows() * out[i / 2 - 1].cols();
      } else {
        in[i / 2 - 1] = move(work);
        assert(in[i / 2 - 1].size() == 3);
        assert(out[i / 2 - 1].rows() == in[i / 2 - 1][0].rows() &&
               out[i / 2 - 1].cols() == in[i / 2 - 1][0].cols());
      }
    }
    assert(in.size() == out.size());
    vector<SimpleVector<num_t> > v;
    v.reserve(cnt);
    for(int i = 0; i < in.size(); i ++)
      for(int j = 0; j < out[i].rows(); j ++)
        for(int k = 0; k < out[i].cols(); k ++) {
          SimpleVector<num_t> work(4);
          for(int m = 0; m < 3; m ++)
            work[m] = in[i][m](j, k);
          work[3] = out[i](j, k);
          v.emplace_back(move(work));
        }
    const auto c(crush<num_t>(v));
    for(int i = 0; i < c.size(); i ++) {
      if(! c[i].first.size()) continue;
      auto vv(makeProgramInvariant<num_t>(c[i].first[0], - num_t(int(1)), true).first);
      for(int j = 1; j < c[i].first.size(); j ++)
        vv += makeProgramInvariant<num_t>(c[i].first[j], - num_t(int(1)), true).first;
      vv /= num_t(c[i].first.size());
      if(vv.dot(vv) != num_t(int(0))) cout << vv;
    }
    cout << endl;
  }
  return 0;
}

