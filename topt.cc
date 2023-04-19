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
    std::string s;
    SimpleVector<num_t> work;
    std::cin >> work;
    while(std::getline(std::cin, s, '\n')) {
      SimpleVector<num_t> lwork(work.size() - 1);
      for(int j = 0; j < abs(len); j ++) {
        for(int i = 0; i < lwork.size(); i ++) lwork[i] = rng();
        for(int i = 0; i < min(lwork.size() - 1, int(s.size())); i ++)
          lwork[i - min(lwork.size() - 1, int(s.size())) + lwork.size()]
            = num_t(s[i - min(lwork.size() - 1, int(s.size())) + int(s.size())]) / num_t(int(256));
        lwork[lwork.size() - 1] = num_t(int(0));
        const char res(int(abs(work.dot(makeProgramInvariant<num_t>(lwork).first)) * num_t(int(256))));
        s += res;
      }
      std::cout << s << std::endl;
    }
  } else if(0 < len) {
    string s;
    string t;
    while(std::getline(std::cin, s, '\n')) t += s;
    SimpleMatrix<num_t> work(t.size() - len + 1, len + 1);
    for(int i = 0; i <= t.size() - len; i ++) {
      SimpleVector<num_t> lwork(len);
      for(int j = 0; j < len; j ++)
        lwork[j] = num_t(t[i + j]) / num_t(int(256));
      auto llwork(makeProgramInvariant<num_t>(lwork));
      work.row(i)  = move(llwork.first);
      work.row(i) *=
        pow(llwork.second, ceil(- log(work.epsilon()) ));
    }
    auto vwork(linearInvariant(work));
    vwork /= - num_t(vwork[vwork.size() - 2]);
    vwork[vwork.size() - 2] = - num_t(int(1));
    cout << vwork;
  }
  return 0;
}

