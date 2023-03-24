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
  const auto m(argv[1][0]);
  if(m == '-') {
    vector<SimpleMatrix<num_t> > L;
    int sz0(0);
    int h(0);
    int w(0);
    std::cin >> sz0;
    std::cin >> h;
    std::cin >> w;
    assert(0 < sz0 && 0 < h && 0 < w);
    L.reserve(3);
    for(int j = 0; j < 3; j ++) {
      SimpleMatrix<num_t> wL(h * w, sz0 * sz0 + 2);
      for(int i = 0; i < wL.rows(); i ++)
        std::cin >> wL.row(i);
      // L.emplace_back(move(wL));
      L.emplace_back(wL);
      assert(L[0].rows() == L[j].rows() && L[0].cols() == L[j].cols());
    }
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > out;
      if(! loadp2or3<num_t>(out, argv[i])) return - 1;
      assert(out[0].rows() * out[0].cols() == sz0 * sz0);
      auto outs(out);
      for(int n = 0; n < outs.size(); n ++)
        outs[n] = SimpleMatrix<num_t>(h, w);
      auto rin(out[0]);
      for(int n = 0; n < rin.rows(); n ++)
        for(int nn = 0; nn < rin.cols(); nn ++)
          rin(n, nn) = argv[1][1] == '0' ? num_t(int(1)) : rng();
      if(argv[1][1] != '0')
        rin = (dft<num_t>(- rin.rows()) * rin.template cast<complex<num_t> >() * dft<num_t>(- rin.cols())).template real<num_t>();
      for(int j = 0; j < out.size(); j ++) {
        cerr << j << " / " << out.size() << " over " << i - 2 << " / " << argc - 2 << endl;
        SimpleVector<num_t> vwork0(out[j].rows() * out[j].cols() + 1);
        for(int n = 0; n < out[j].rows(); n ++)
          for(int nn = 0; nn < out[j].cols(); nn ++)
            vwork0[n * out[j].cols() + nn] = out[j](n, nn) * rin(n, nn);
        vwork0[vwork0.size() - 1] = num_t(int(0));
        auto outwork(L[j] * makeProgramInvariant<num_t>(vwork0).first);
        for(int n = 0; n < outwork.size(); n ++)
          outs[j](n / outs[j].cols(), n % outs[j].cols()) =
            revertProgramInvariant<num_t>(make_pair(outwork[n], num_t(int(1)) ));
      }
      if(! savep2or3<num_t>(argv[i], outs) )
        cerr << "failed to save." << endl;
    }
  } else if(m == '+' || m == '0') {
    vector<vector<SimpleMatrix<num_t> > > in;
    vector<vector<SimpleMatrix<num_t> > > noise;
    in.resize(argc - 2);
    noise.resize(in.size());
          int sz(0);
    const int num(argv[1][1] == '+' ? sqrt(num_t(in.size())) : (argv[1][0] == '0' ? num_t(int(1)) : log(num_t(in.size())) / log(num_t(int(2)))));
    for(int i = 2; i < argc; i ++) {
      if(! loadp2or3<num_t>(in[i - 2], argv[i])) continue;
      assert(in[0][0].rows() == in[i - 2][0].rows() &&
             in[0][0].cols() == in[i - 2][0].cols());
      if(i == 2) sz = int(sqrt(num_t(min(int(in.size()), int(in[i - 2][0].rows())))));
      noise[i - 2].resize(num, SimpleMatrix<num_t>(sz, sz));
      for(int j = 0; j < num; j ++) {
        for(int n = 0; n < sz; n ++)
          for(int nn = 0; nn < sz; nn ++)
            noise[i - 2][j](n, nn) = argv[1][0] == '0' && argv[1][1] == '0' ?
              num_t(int(1)) : rng();
        if(argv[1][1] != '0')
          noise[i - 2][j] = (dft<num_t>(- noise[i - 2][j].rows()) * noise[i - 2][j].template cast<complex<num_t> >() * dft<num_t>(- noise[i - 2][j].cols())).template real<num_t>();
      }
    }
    cout << sz << endl;
    cout << in[0][0].rows() << endl;
    cout << in[0][0].cols() << endl;
    auto shrink(in);
    for(int i = 0; i < shrink.size(); i ++)
      for(int j = 0; j < shrink[i].size(); j ++) {
        shrink[i][j] = SimpleMatrix<num_t>(sz, sz).O();
        for(int ii = 0; ii < sz; ii ++)
          for(int jj = 0; jj < sz; jj ++) {
            int cnt(0);
            for(int iik = 0;
                iik <= min(in[i][j].rows() / sz - 1,
                  in[i][j].rows() - ii * (in[i][j].rows() / sz)); iik ++)
              for(int jjk = 0; jjk <= min(in[i][j].cols() / sz - 1,
                    in[i][j].cols() - jj * (in[i][j].cols() / sz));
                  jjk ++, cnt ++)
                shrink[i][j](ii, jj) +=
                  in[i][j](ii * (in[i][j].rows() / sz) + iik,
                           jj * (in[i][j].cols() / sz) + jjk);
            shrink[i][j](ii, jj) /= num_t(cnt);
          }
      }
    for(int j = 0; j < in[0].size(); j ++)
      for(int m = 0; m < in[0][0].rows() * in[0][0].cols(); m ++){
        cerr << j * in[0][0].rows() * in[0][0].cols() + m << " / " << in[0][0].rows() * in[0][0].cols() * in[0].size() << endl;
        SimpleMatrix<num_t> work(num * in.size(), shrink[0][0].rows() * shrink[0][0].cols() + 2);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
        for(int i = 0; i < in.size(); i ++)
          for(int jj = 0; jj < num; jj ++) {
            SimpleVector<num_t> vwork(shrink[i][j].rows() * shrink[i][j].cols() + 1);
            for(int n = 0; n < shrink[i][j].rows(); n ++)
              for(int nn = 0; nn < shrink[i][j].cols(); nn ++)
                vwork[n * shrink[i][j].cols() + nn] = shrink[i][j](n, nn) * noise[i][jj](n, nn);
            vwork[vwork.size() - 1] = in[i][j](m / in[0][0].cols(), m % in[0][0].cols());
            auto mpi(makeProgramInvariant<num_t>(vwork));
            work.row(i * num + jj)  = move(mpi.first);
            work.row(i * num + jj) *=
              pow(mpi.second, ceil(- log(in[0][0].epsilon()) ));
          }
        auto vwork(linearInvariant(work));
        vwork /= - num_t(vwork[vwork.size() - 2]);
        vwork[vwork.size() - 2] = num_t(int(0));
        cout << vwork;
      }
  }
  return 0;
}

