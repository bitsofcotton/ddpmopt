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

static inline num_t rng() {
  myuint res(0);
  // XXX: we don't trust system or compiler PRNG.
  // static std::random_device rd;
  // XXX: we want natural, deterministic, better PRNG, however,
  //      we don't search deepinside of this PRNG.
  //      (might be predecessor exists.)
  static uint64_t t(1);
  assert(t && "rng() should not be periodical.");
#if defined(_FLOAT_BITS_)
  for(int i = 0; i < _FLOAT_BITS_ / sizeof(uint32_t) / 8; i ++) {
#else
  for(int i = 0; i < 2; i ++) {
#endif
    res <<= sizeof(uint32_t) * 8;
#if defined(_FLOAT_BITS_)
typedef SimpleFloat<uint64_t, unsigned __int128, 64, int64_t> thisfl;
#else
typedef long double thisfl;
#endif
    auto buf(sin(thisfl(t ++)) * pow(thisfl(int(2)), thisfl(int(32))));
    buf  -= floor(buf);
    res  |= uint32_t(int(buf * pow(thisfl(int(2)), thisfl(int(32)) )));
#undef thisfl
    // res  |= uint32_t(rd());
  }
  return max(num_t(int(0)), min(num_t(int(1)), num_t(res) / num_t(~ myuint(0)) ));
}

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  assert(1 < argc);
  const auto len(std::atoi(argv[1]));
  if(len < 0) {
    int size(0);
    std::cin >> size;
    std::string s;
    vector<SimpleVector<num_t> > work;
    work.resize(size);
    for(int j = 0; j < size; j ++)
      std::cin >> work[j];
    while(std::getline(std::cin, s, '\n')) {
      SimpleVector<num_t> lwork(work[0].size() - 1);
      for(int j = 0; j < abs(len); j ++) {
        for(int i = 0; i < lwork.size(); i ++) lwork[i] = rng();
        for(int i = 0; i < min(lwork.size() - 1, int(s.size())); i ++)
          lwork[i - min(lwork.size() - 1, int(s.size())) + lwork.size()]
            = num_t(s[i - min(lwork.size() - 1, int(s.size())) + int(s.size())]) / num_t(int(256));
        lwork[lwork.size() - 1] = num_t(int(0));
        int   Midx(- 1);
        num_t M(int(0));
        for(int k = 0; k < work.size(); k ++) {
          const auto lM(abs(work[k].dot(makeProgramInvariant<num_t>(lwork).first) / sqrt(work[k].dot(work[k])) ));
          if(Midx < 0 || M < lM) {
            M = lM;
            Midx = k;
          }
        }
        s += char(int(abs(work[Midx].dot(makeProgramInvariant<num_t>(lwork).first)) * num_t(int(256))));
      }
      std::cout << s << std::endl;
    }
  } else if(0 < len) {
    string s;
    string t;
    while(std::getline(std::cin, s, '\n')) t += s;
    vector<SimpleVector<num_t> > work;
    work.reserve(t.size() - len + 1);
    for(int i = 0; i <= t.size() - len; i ++) {
      SimpleVector<num_t> lwork(len);
      for(int j = 0; j < len; j ++)
        lwork[j] = num_t(t[i + j]) / num_t(int(256));
      work.emplace_back(lwork);
    }
    auto vwork(crush<num_t>(work));
    std::cout << vwork.size() << std::endl;
    for(int i = 0; i < vwork.size(); i ++) {
      SimpleMatrix<num_t> lwork(vwork[i].first.size(), len + 1);
      for(int j = 0; j < lwork.rows(); j ++) {
        auto llwork(makeProgramInvariant<num_t>(vwork[i].first[j]));
        lwork.row(j)  = move(llwork.first);
        lwork.row(j) *=
          pow(llwork.second, ceil(- log(lwork.epsilon()) ));
      }
      if(lwork.rows() <= lwork.cols() + 1) {
        for(int j = 1; j < lwork.rows(); j ++)
          lwork.row(0) += lwork.row(j);
        lwork.row(0) /= - num_t(lwork(0, lwork.cols() - 2));
        lwork(0, lwork.cols() - 2) = num_t(int(0));
        cout << lwork.row(0);
      } else {
        auto vvwork(linearInvariant(lwork));
        vvwork /= - num_t(vvwork[vvwork.size() - 2]);
        vvwork[vvwork.size() - 2] = num_t(int(0));
        cout << vvwork;
      }
    }
  }
  return 0;
}

