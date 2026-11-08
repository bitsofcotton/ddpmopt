#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <cctype>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(_MIMALLOC_)
#include <mimalloc.h>
#endif

#if defined(_P_VULKAN_)
#include <vulkan/vulkan.h>
#endif

#if !defined(_OLDCPP_) && defined(_PERSISTENT_)
# if !defined(_FLOAT_BITS_)
#  define int ssize_t
# elif _FLOAT_BITS_ == 64
#  define int int32_t
# elif _FLOAT_BITS_ == 128
#  define int int64_t
# endif
#else
# define int int64_t
#endif
#include "lieonn.hh"
typedef myfloat num_t;
lieonn_t lieonn;

using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::vector;

#include <stdlib.h>

template <typename T> SimpleMatrix<T> unOffsetHalf(const SimpleMatrix<T>& m) {
  SimpleMatrix<T> res(m);
  res.entity = unOffsetHalf<T>(res.entity);
  return res;
}

template <typename T> vector<SimpleMatrix<T> > unOffsetHalf(const vector<SimpleMatrix<T> >& m) {
  vector<SimpleMatrix<T> > res(m);
  for(int i = 0; i < res.size(); i ++) res[i] = unOffsetHalf<T>(res[i]);
  return res;
}

#undef int
int main(int argc, const char* argv[]) {
#if !defined(_OLDCPP_) && defined(_PERSISTENT_)
# if !defined(_FLOAT_BITS_)
#  define int ssize_t
# elif _FLOAT_BITS_ == 64
#  define int int32_t
# elif _FLOAT_BITS_ == 128
#  define int int64_t
# endif
#else
# define int int64_t
#endif
  const char& m(argv[1][0]);
  lieonnStaticInit();
  if(argc <= 1) goto usage;
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  if(m == 'p' || m == 'T') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 1);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in.emplace_back(move(work));
      assert(in[in.size() - 1][0].rows() == in[0][0].rows() &&
        in[in.size() - 1][0].cols() == in[0][0].cols());
    }
    if(in[0].size() == 3)
      for(int i = 0; i < in.size(); i ++) in[i] = rgb2xyz<num_t>(in[i]);
    in = normalize<num_t>(in);
    if(in[0].size() == 3)
      for(int i = 0; i < in.size(); i ++) in[i] = xyz2rgb<num_t>(in[i]);
    if(m == 'T') {
      if(! savep2or3<num_t>("testg.ppm", predMat<num_t, 20>(move(in), 3,
        string(" test")) ))
          cerr << "failed to save test." << endl;
    } else {
      if(! savep2or3<num_t>("predg.ppm", normalize<num_t>(stretchPred<num_t>(
        normalize<num_t>(predMat<num_t, 20>(move(in), 3, string(" pred"))) )) ))
          cerr << "failed to save whole pred." << endl;
    }
  } else if(m == 'q') {
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      work = normalize<num_t>(work);
      int seed(loop22<num_t>() + infBase() + 1);
      int seed1;
      for(seed1 = 1;
        seed1 <= int(log(num_t(seed + seed1)) / log(num_t(2))); seed1 ++) ;
      seed += seed1;
      const int px(work[0].rows() / seed / 2);
      assert(0 < px);
      vector<SimpleMatrix<num_t> > wwork(work.size());
      for(int j = 0; j < wwork.size(); j ++)
        wwork[j].resize(work[0].rows() + px, work[0].cols()).O(
          ).setMatrix(0, 0, work[j]);
      SimpleVector<vector<SimpleVector<num_t> > > v(work[0].rows());
      for(int i = 0; i < v.size(); i ++) {
        v[i].resize(work.size());
        for(int j = 0; j < v[i].size(); j ++) v[i][j] = move(work[j].row(i));
      }
      for(int j = 0; j < px; j ++) {
        SimpleVector<SimpleVector<num_t> > q(predVec<num_t, 20, true>(skipX<
          vector<SimpleVector<num_t> > >(v, j + 1), 2, to_string(j) +
            string("/") + to_string(px) ) );
        for(int i = 0; i < wwork.size(); i ++)
          wwork[i].row(work[i].rows() + j) = move(q[i]);
      }
      if(! savep2or3<num_t>(argv[i0], move(wwork)) )
        cerr << "failed to save." << endl;
    }
  } else if(m == 'h') {
    bool met(false);
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      if(work.size() != 1) work[0] = rgb2d<num_t>(work);
      work.resize(1);
      work = normalize<num_t>(work);
      if(! met) {
        std::cout << "const int cg_width = " << work[0].cols() << ";" << endl;
        std::cout << "const int cg_height = " << work[0].rows() << ";" << endl;
        std::cout << "const int cg_n = " << argc - 2 << ";" << endl;
        std::cout << "const char cg[" << argc - 2 << "][" << work[0].rows() * work[0].cols() << "] = {" << endl;
        std::cout << "{";
      } else std::cout << endl << "{";
      for(int i = 0; i < work[0].rows(); i ++)
        for(int j = 0; j < work[0].cols(); j ++)
          std::cout << int(work[0](i, j) * num_t(int(255))) << ",";
      std::cout << "},";
      met = true;
    }
    std::cout << "};" << endl;
  } else if(m == 'H') {
    // cf. sox in.mp3 -c 1 -b 16 -e signed -r 65536 a.raw
    //     ./ddpmopt H ... < a.raw > b.raw
    //     sox -M -b 16 -e signed -r 65536 a.raw -b 16 -e signed -r 65536
    //       b.raw out.wav
    // cf. "hemi-sync" on some surface search.
    // XXX: political or patent matter on publish?
    const int blocks(65536);
    int shift(atoi(argv[2]));
    SimpleVector<int16_t> v(blocks);
    SimpleVector<complex(num_t)> f;
    while(! std::cin.eof() && ! std::cin.bad()) {
      std::cin.read(reinterpret_cast<char*>(&v[0]), sizeof(int16_t) * blocks);
      f = fft<num_t>(v.template cast<num_t>().template cast<complex(num_t)>());
      for(int i = 0; i < f.size() - shift; i ++)
        f[f.size() - 1 - i] = f[f.size() - shift - i - 1];
      v = ifft<num_t>(f).template real<num_t>().template cast<int>().template cast<int16_t>();
      std::cout.write(reinterpret_cast<char*>(&v[0]), sizeof(int16_t) * blocks);
    }
  } else if(m == 'G') {
    // cf. sox in.mp3 -c 1 -b 16 -e signed -r 65536 a.raw
    //     ./ddpmopt G < a.raw > b.raw
    //     sox -b 16 -e signed -r 65536 b.raw out.wav
    // cf. 440 Hz : 432 Hz
    // XXX: political or patent matter on publish?
    const int blocks(65536);
    SimpleVector<int16_t> v(blocks);
    SimpleVector<complex(num_t)> f;
    while(! std::cin.eof() && ! std::cin.bad()) {
      std::cin.read(reinterpret_cast<char*>(&v[0]), sizeof(int16_t) * blocks);
      f = fft<num_t>(v.template cast<num_t>().template cast<complex(num_t)>());
      const SimpleVector<complex(num_t)> f0(f);
      int i;
      for(i = 0; i < f.size(); i ++) {
        const int idx(exp(log(num_t(i + 1)) * log(num_t(55)) / log(num_t(54)) ));
        if(f.size() <= idx) break;
        f[i] = f0[idx];
      }
      for( ; i < f.size(); i ++)
        f[i] = complexctor(num_t)(num_t(int(0)));
      v = ifft<num_t>(f).template real<num_t>().template cast<int>().template cast<int16_t>();
      std::cout.write(reinterpret_cast<char*>(&v[0]), sizeof(int16_t) * blocks);
    }
  } else if(m == 'C') {
    const int len(atoi(argv[2]));
    if(!len) {
      std::cout << "const float sqe = " << sqrt(SimpleMatrix<num_t>().epsilon() ) << ";" << endl;
      std::cout << "const float denom = " << (num_t(int(1)) + sqrt(sqrt(SimpleMatrix<num_t>().epsilon() )) ) << ";" << endl;
    } else {
      std::cout << "float[" << len << "][" << len << "](" << endl;
      for(int i = 0; i < len; i ++) {
        const SimpleVector<num_t>& pn(pnextcacher<num_t>(i + 1, 1));
        int j;
        std::cout << "float[" << len << "](";
        for(j = 0; j < pn.size() - 1; j ++) std::cout << pn[j] << ", ";
        std::cout << pn[j ++];
        if(pn.size() != len) {
          std::cout << ", ";
          for( ; j < len - 1; j ++) std::cout << num_t(int(0)) << ", ";
          std::cout << num_t(int(0));
        }
        std::cout << ")," << endl;
      }
      std::cout << ");" << endl << flush;
    }
  } else if(m == 'd') {
    vector<SimpleMatrix<num_t> > in0;
    vector<SimpleMatrix<num_t> > in1;
    loadp2or3<num_t>(in0, argv[2]);
    loadp2or3<num_t>(in1, argv[3]);
    assert(in0.size() == in1.size());
    for(int i = 0; i < in0.size(); i ++)
      for(int j = 0; j < in0[i].rows(); j ++) {
        for(int k = 0; k < in1[i].cols(); k ++)
          in1[i](j, k) = sgn<num_t>(unOffsetHalf<num_t>(in1[i](j, k)));
        in0[i].row(j) = offsetHalf<num_t>(unOffsetHalf<num_t>(in0[i].row(j)
          ).tdot(in1[i].row(j)) );
      }
    savep2or3("tdot.ppm", in0);
  } else goto usage;
  cerr << "Done" << endl;
  lieonnStaticDestroy();
  return 0;
 usage:
  lieonnStaticDestroy();
  cerr << "Usage:" << endl;
  cerr << "# predict following image" << endl;
  cerr << argv[0] << " p <in0.ppm> ..." << endl;
  cerr << "# check predict last image" << endl;
  cerr << argv[0] << " T <in0.ppm> ..." << endl;
  cerr << "# predict down scanlines" << endl;
  cerr << argv[0] << " q <in0out.ppm> ..." << endl;
  return - 1;
}

