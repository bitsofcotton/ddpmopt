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
  auto len(1 < argc ? std::atoi(argv[1]) : 0);
  if(len < 0) {
    int size(0);
    std::cin >> size;
    std::string s;
    vector<SimpleVector<num_t> > work;
    work.reserve(size);
    for(int j = 0; j < size; j ++) {
      SimpleVector<num_t> vwork;
      std::cin >> vwork;
      vwork[vwork.size() - 2] = num_t(int(0));
      work.emplace_back(vwork /= sqrt(vwork.dot(vwork)));
    }
    while(std::getline(std::cin, s, '\n')) {
      SimpleVector<num_t> lwork(work[0].size() - 1);
      for(int j = 0; j < abs(len); j ++) {
        for(int i = 0; i < lwork.size(); i ++) lwork[i] = rng();
        for(int i = 0; i < min(lwork.size() - 1, int(s.size())); i ++)
          lwork[i - min(lwork.size(), int(s.size()) + 1) + lwork.size()]
            = num_t(s[i - min(lwork.size() - 1, int(s.size())) + int(s.size())]) / num_t(int(256));
        lwork[lwork.size() - 1] = num_t(int(0));
        int   Midx(- 1);
        num_t M(int(0));
        const auto pinv(makeProgramInvariant<num_t>(lwork));
        for(int k = 0; k < work.size(); k ++) {
          const auto lM(abs(work[k].dot(pinv.first) ));
          if(Midx < 0 || M < lM) {
            M    = lM;
            Midx = k;
          }
        }
        assert(0 <= Midx && Midx < work.size() );
        s += char(int(revertProgramInvariant<num_t>(make_pair(work[Midx].dot(pinv.first), num_t(int(1)) /* pinv.second */)) * num_t(int(256)) ));
      }
      std::cout << s << std::endl;
    }
  } else if(0 <= len) {
    string s;
    string t;
    while(std::getline(std::cin, s, '\n')) t += s;
    if(! len) len = min(int(6), int(ceil(pow(num_t(int(t.size())), num_t(int(1)) / num_t(int(6)) ))) + 1);
    vector<SimpleVector<num_t> > work;
    work.reserve(t.size() - len + 1);
    for(int i = 0; i <= t.size() - len; i ++) {
      SimpleVector<num_t> lwork(len);
      for(int j = 0; j < len; j ++)
        lwork[j] = num_t(t[i + j]) / num_t(int(256));
      work.emplace_back(lwork);
    }
    auto vwork(crush<num_t>(work, work[0].size(), work.size()));
    std::cout << vwork.size() << std::endl;
    for(int i = 0; i < vwork.size(); i ++) {
      if(! (vwork[i].first.size() <= len + 1)) cerr << "!" << i << endl;
      SimpleVector<num_t> lwork(len + 1);
      lwork.O();
      for(int j = 0; j < vwork[i].first.size(); j ++)
        lwork += makeProgramInvariant<num_t>(vwork[i].first[j]).first;
      cout << lwork;
    }
  }
  return 0;
}

