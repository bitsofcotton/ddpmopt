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

using std::cout;
using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::to_string;

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
//#define int int32_t
  assert(1 < argc);
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  vector<vector<SimpleMatrix<num_t> > > in;
  in.reserve(argc - 1);
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
  }
  {
    const auto p(predMat<num_t>(in = normalize<num_t>(in)));
    assert(p.size() == 1);
    if(! savep2or3<num_t>("predg.ppm",
        normalize<num_t>(p[0].size() == 3 ? xyz2rgb<num_t>(p[0]) : p[0])) )
          cerr << "failed to save." << endl;
  }
  in.reserve(argc - 1);
  for(int i = 1; i < argc; i ++) {
    vector<SimpleMatrix<num_t> > work;
    if(! loadp2or3<num_t>(work, argv[i])) continue;
    in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
  }
  {
    const auto p(predMat<num_t>(in = normalize<num_t>(in)));
    for(int i = 0; i < p.size(); i ++)
      if(! savep2or3<num_t>((string("predg-") + to_string(i) +
        string(".ppm")).c_str(),
          normalize<num_t>(p[i].size() == 3 ? xyz2rgb<num_t>(p[i]) : p[i])) )
            cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

