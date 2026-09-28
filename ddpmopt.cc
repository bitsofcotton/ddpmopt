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
#define MIMALLOC_OVERRIDE_H
#define MIMALLOC_NEW_DELETE_H
#include <mimalloc.h>
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
  if(m == 'p') {
    vector<vector<SimpleMatrix<num_t> > > in0;
    in0.reserve(argc - 1);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in0.emplace_back(move(work));
    }
    int len(absceil(log(num_t(int(in0.size() ))) / log(num_t(int(2))) ));
    int py(1);
    int px(1);
    int ppy(1);
    int ppx(1);
    for(int i = 1; i < max(in0[0][0].rows(), in0[0][0].cols()) || !py || !px;
        i ++) {
      ppy = py; ppx = px;
      if(in0[0][0].rows() < in0[0][0].cols()) {
        py = i * in0[0][0].rows() / in0[0][0].cols();
        px = i;
      } else {
        py = i;
        px = i * in0[0][0].cols() / in0[0][0].rows();
      }
      if(! (px * py <= len)) break;
    }
    py = ppy; px = ppx;
    cerr << "internal(" << py << ", " << px << ")" << endl;
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(in0.size());
    for(int i = 0; i < in0.size(); i ++) {
      vector<SimpleMatrix<num_t> > work;
      work.reserve(in0[i].size());
      for(int j = 0; j < in0[i].size(); j ++) work.emplace_back( (
        dftcache<num_t>(- py) * dftcache<num_t>(in0[i][j].rows()).subMatrix(0,
          0, py, in0[i][j].rows()) * in0[i][j].template cast<complex(num_t)>(
          ) * (dftcache<num_t>(- px) * dftcache<num_t>(in0[i][j].cols()
            ).subMatrix(0, 0, px, in0[i][j].cols() )).transpose()
              ).template real<num_t>());
      in.emplace_back(move(work));
    }
    in0 = normalize<num_t>(in0);
    in  = normalize<num_t>(in);
    vector<vector<SimpleMatrix<num_t> > > win(const_cast<vector<vector<
      SimpleMatrix<num_t> > >&>(in));
    if(! savep2or3<num_t>("testg.ppm", predMatTangleLast<num_t, 20, true>(
      move(win), 3, string(" ") + string(argv[0]) + string(" ") +
        string(argv[1]) ) ) )
      cerr << "failed to save test." << endl;
    win.resize(0);
    const int ry(in0[0][0].rows() / py);
    const int rx(in0[0][0].cols() / px);
    const int row(in0[0][0].rows());
    const int col(in0[0][0].cols());
    const int rowp(in[0][0].rows());
    const int colp(in[0][0].cols());
    if(! savep2or3<num_t>("predg.ppm", normalize<num_t>(stretchPred<num_t>(
      enlargeApply<num_t>(ry, rx, rowp, colp, col, normalize<num_t>(
        predMat<num_t, 20>(move(in), 3, string(" ") + string(argv[0]) +
          string(" ") + string(argv[1]) ) ), enlargeApply0<num_t>(
            enlargePrep<num_t, 40>(ry, rx, row, col, move(in0) )) ),
              in[0].size() == 1 ? 15 : 5)) ))
      cerr << "failed to save whole pred." << endl;
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
      const int u_step(work[0].rows() / seed / 2);
      assert(0 < u_step);
      const int px(absceil(log(num_t(u_step)) / log(int(2)) ));
      assert(0 < px);
      SimpleVector<vector<SimpleVector<num_t> > > pwork(work[0].rows());
      for(int i = 0; i < pwork.size(); i ++) {
        pwork[i].reserve(work.size());
        for(int j = 0; j < work.size(); j ++)
          pwork[i].emplace_back(((dft<num_t>(- px) * dft<num_t>(work[j].cols()
            ).subMatrix(0, 0, px, work[j].cols()) ) * work[j].row(i
              ).template cast<complex(num_t)>() ).template real<num_t>());
      }
      vector<SimpleMatrix<num_t> > pred;
      pred.resize(work.size(),
        SimpleMatrix<num_t>(u_step, pwork[0][0].size()).O());
      for(int j = 0; j < pred[0].rows(); j ++) {
        SimpleVector<SimpleVector<num_t> > q(predVec<num_t, 20, true>(
          skipX<vector<SimpleVector<num_t> > >(pwork, j + 1), 2, to_string(j) +
            string("/") + to_string(pred[0].rows()) ) );
        for(int i = 0; i < pred.size(); i ++) pred[i].row(j) = move(q[i]);
      }
      const int ry(1);
      const int rx(work[0].cols() / px);
      const int row(work[0].rows());
      const int col(work[0].cols());
      const int rowp(pred[0].rows());
      const int colp(pred[0].cols());
      vector<SimpleMatrix<num_t> > wwork(work.size());
      for(int j = 0; j < wwork.size(); j ++)
        wwork[j].resize(work[0].rows() + pred[0].rows(), work[0].cols()).O(
          ).setMatrix(0, 0, work[j]);
      vector<vector<SimpleMatrix<num_t> > > ww;
      ww.resize(1, move(work));
      pred = enlargeApply<num_t>(ry, rx, rowp, colp, col, normalize<num_t>(
        pred), enlargeApply0<num_t>(enlargePrep<num_t, 40>(ry, rx, rowp, col,
          move(ww) )) );
      for(int j = 0; j < wwork.size(); j ++)
        wwork[j].setMatrix(row + j, 0, pred[j]);
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
  } else goto usage;
  cerr << "Done" << endl;
  lieonnStaticDestroy();
  return 0;
 usage:
  lieonnStaticDestroy();
  cerr << "Usage:" << endl;
  cerr << "# predict following image" << endl;
  cerr << argv[0] << " p <in0.ppm> ..." << endl;
  cerr << "# predict down scanlines" << endl;
  cerr << argv[0] << " q <in0out.ppm> ..." << endl;
  return - 1;
}

