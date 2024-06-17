#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>
#include <assert.h>

#include "lieonn.hh"
typedef myfloat num_t;

int main(int argc, const char* argv[]) {
  assert(2 < argc);
  const auto len(std::atoi(argv[1]));
  const auto step(std::atoi(argv[2]));
  const auto r(argc <= 3 ? int(2) : std::atoi(argv[3]));
  assert(3 <= len && len <= 7);
  const auto  d(diff<num_t>(len).row(len - 1));
  const auto& n(pnextcacher<num_t>(len, step, r));
  num_t sumn(int(0));
  auto  sumd(sumn);
  auto  n2n(sumn);
  auto  n2d(sumn);
  for(int i = 0; i < len; i ++) {
    std::cout << n[i] << ", " << d[i] << std::endl;
    sumn += n[i];
    sumd += d[i];
    n2n  += n[i] * n[i];
    n2d  += d[i] * d[i];
  }
  std::cout << sumn << ", " << sumd << ", " << sqrt(n2n) << ", " << sqrt(n2d) << std::endl;
  return 0;
}

