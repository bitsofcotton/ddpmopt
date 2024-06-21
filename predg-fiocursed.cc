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
using std::vector;
using std::sort;
using std::binary_search;
using std::make_pair;
using std::istringstream;
using std::move;
using std::pair;

const char* input[] = {
#include "input.ppm.string.array.txt"
};

#undef int
int main(int argc, const char* argv[]) {
#define int int64_t
//#define int int32_t
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  vector<vector<SimpleMatrix<num_t> > > in;
  for(int i = 0; i < sizeof(input) / sizeof(char*); i ++) {
    vector<SimpleMatrix<num_t> > work;
    istringstream line(input[i]);
    if(! loadp2or3<num_t>(work, static_cast<istream&>(line))) continue;
    if(work.size() == 3)
      work = rgb2xyz<num_t>(work);
    in.emplace_back(move(work));
  }
  const int  color(in.size());
        auto p(predMat<num_t>(normalize<num_t>(in)));
  for(int i = 0; i < p.first.size(); i ++) {
    if(! savep2or3<num_t>((std::string("predg-forward-") + std::to_string(i) + std::string(".ppm")).c_str(), p.first[i].size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(p.first[i])) : p.first[i], color) )
      cerr << "failed to save." << endl;
    if(! savep2or3<num_t>((std::string("predg-backward-") + std::to_string(i) + std::string(".ppm")).c_str(), p.second[i].size() == 3 ? normalize<num_t>(xyz2rgb<num_t>(p.second[i])) : p.second[i], color) )
      cerr << "failed to save." << endl;
  }
  cerr << " Done" << endl;
  return 0;
}

