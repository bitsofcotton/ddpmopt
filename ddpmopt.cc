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

#if defined(_P_VULKAN_)
#include <vulkan/vulkan.h>
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
  if(m == '-') {
    vector<vector<SimpleVector<num_t> > > L;
    string s;
    for(int len = 0; 0 <= len && std::getline(std::cin, s, '\n'); ) {
      if(s[0] == '*') {
        L.emplace_back(vector<SimpleVector<num_t> >());
        len ++;
        continue;
      } else if(! s.size()) continue;
      SimpleVector<num_t> l;
      stringstream ins(s);
      ins >> l;
      L[len - 1].emplace_back(move(l));
    }
    L = enlargeApply0<num_t>(move(L));
    for(int i0 = 2; i0 < argc; i0 ++) {
      cerr << i0 - 2 << " / " << argc - 2 << endl;
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i0])) return - 1;
      if(argv[1][1] == '\0' && in.size() != 3) {
        cerr << argv[i0] << " doesn't include 3 colors" << endl;
        continue;
      }
      int sq(sqrt(num_t(L.size() / in.size())));
      if((sq + 1) * (sq + 1) == L.size()) sq ++;
      vector<SimpleMatrix<num_t> > out(enlargeApply(atoi(&argv[1][1]),
        atoi(&argv[1][1]), sq / atoi(&argv[1][1]), sq / atoi(&argv[1][1]),
          in, L));
      if(! savep2or3<num_t>((string(argv[i0]) + string(".pgm")).c_str(), out) )
        cerr << "failed to save." << endl;
    }
  } else if(m == '+') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc);
    for(int i = 3; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in.emplace_back(move(work));
    }
    vector<vector<SimpleVector<num_t> > > v(enlargePrep(atoi(&argv[1][1]),
      atoi(&argv[1][1]), atoi(argv[2]), atoi(argv[2]), move(in)) );
    for(int i0 = 0; i0 < v.size(); i0 ++) {
      std::cout << "*" << endl;
      for(int i = 0; i < v[i0].size(); i ++) std::cout << v[i0][i];
      std::cout << endl;
    }
  } else if(m == 'p' || m == 'T') {
    vector<vector<SimpleMatrix<num_t> > > in0;
    in0.reserve(argc - 1);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in0.emplace_back(move(work));
    }
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(in0.size());
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
    if(m == 'T') {
      if(! savep2or3<num_t>("testg.ppm", predMatTangleLast<num_t, 20>(
        move(in), 3, string(" ") + string(argv[0]) + string(" ") +
          string(argv[1]) ) ) )
        cerr << "failed to save." << endl;
    } else {
      const int ry(in0[0][0].rows() / py);
      const int rx(in0[0][0].cols() / px);
      const int row(in0[0][0].rows());
      const int col(in0[0][0].cols());
      const int rowp(in[0][0].rows());
      const int colp(in[0][0].cols());
      if(! savep2or3<num_t>("predg.ppm", normalize<num_t>(stretchPred<num_t>(
        enlargeApply<num_t>(ry, rx, rowp, colp, normalize<num_t>(
          predMat<num_t, 20>(move(in), 3, string(" ") + string(argv[0]) +
            string(" ") + string(argv[1]) ) ), enlargeApply0<num_t>(
              enlargePrep<num_t>(ry, rx, row, col, move(in0) )) ),
                in[0].size() == 1 ? 15 : 5)) ))
        cerr << "failed to save." << endl;
    }
  } else if(m == 'q') {
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      work = normalize<num_t>(work);
      SimpleVector<vector<SimpleVector<num_t> > > pwork(work[0].rows());
      for(int i = 0; i < pwork.size(); i ++) {
        pwork[i].reserve(work.size());
        for(int j = 0; j < work.size(); j ++)
          pwork[i].emplace_back(work[j].row(i));
      }
      const int step(work[0].rows() / (loop22<num_t>() + infBase() + 1) );
      vector<SimpleMatrix<num_t> > wwork;
      wwork.resize(work.size(),
        SimpleMatrix<num_t>(work[0].rows() + step, work[0].cols()).O());
      for(int j = 0; j < wwork.size(); j ++)
        wwork[j].setMatrix(0, 0, work[j]);
      for(int j = 0; j < wwork[0].rows() - work[0].rows(); j ++) {
        SimpleVector<SimpleVector<num_t> > q(predVec<num_t, 20, true>(
          skipX<vector<SimpleVector<num_t> > >(pwork, j + 1), 2, to_string(j) +
            string("/") + to_string(wwork[0].rows() - work[0].rows()) )  );
        for(int i = 0; i < wwork.size(); i ++)
          wwork[i].row(work[0].rows() + j) = move(q[i]);
      }
      if(! savep2or3<num_t>(argv[i0], move(wwork)) )
        cerr << "failed to save." << endl;
    }
  } else if(m == 'x' || m == 'y' || m == 'i' || m == 't') {
    vector<num_t> score;
    score.resize(argc + 1, num_t(int(0)));
    switch(argv[1][0]) {
    case 'x':
    case 'i':
      for(int i0 = 2; i0 < argc; i0 ++) {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i0])) continue;
        for(int i = 0; i < work.size(); i ++)
          for(int ii = 0; ii < work[i].rows(); ii ++) {
            idFeeder<num_t> w(3);
            for(int jj = 0; jj < work[i].cols(); jj ++) {
              if(w.full) {
                const num_t pp(p0maxNext<num_t>(w.res));
                score[i0] += (pp - work[i](ii, jj)) * (pp - work[i](ii, jj));
              } 
              w.next(work[i](ii, jj));
            }
          }
        if(argv[1][0] == 'x')
          score[i0] /= num_t(work[0].rows() * work[0].cols() * work.size());
      }
      if(argv[1][0] == 'x') break;
    case 'y':
      for(int i0 = 2; i0 < argc; i0 ++) {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i0])) continue;
        for(int i = 0; i < work.size(); i ++)
          for(int jj = 0; jj < work[i].cols(); jj ++) {
            idFeeder<num_t> w(3);
            for(int ii = 0; ii < work[i].rows(); ii ++) {
              if(w.full) {
                const num_t pp(p0maxNext<num_t>(w.res));
                score[i0] += (pp - work[i](ii, jj)) * (pp - work[i](ii, jj));
              } 
              w.next(work[i](ii, jj));
            }
          }
        score[i0] /= num_t(work[0].rows() * work[0].cols() * work.size());
      }
      break;
    case 't':
      {
        vector<vector<SimpleMatrix<num_t> > > b;
        b.resize(argc + 1);
        for(int i = 2; i < argc; i ++) {
          vector<SimpleMatrix<num_t> > work;
          if(! loadp2or3<num_t>(work, argv[i])) continue;
          b[i] = move(work);
          assert(b[i].size() == b[2].size() &&
            b[i][0].rows() == b[2][0].rows() &&
            b[i][0].cols() == b[2][0].cols());
        }
        int cnt(0);
        for(int i = 0; i < b[2][0].rows(); i ++)
          for(int j = 0; j < b[2][0].cols(); j ++)
            for(int k = 0; k < b[2].size(); k ++) {
              idFeeder<num_t> w(3);
              for(int ii = 2; ii < b.size(); ii ++, cnt ++) {
                if(! b[ii].size()) continue;
                if(w.full) {
                  const num_t pp(p0maxNext<num_t>(w.res));
                  score[0] += (pp - b[ii][k](i, j)) * (pp - b[ii][k](i, j));
                }
                w.next(b[ii][k](i, j));
              }
            }
        score[0] /= num_t(cnt);
      }
      break;
    }
    if(argv[1][0] == 'x' || argv[1][0] == 'y' || argv[1][0] == 'i')
      for(int i = 2; i < argc; i ++)
        std::cout << sqrt(score[i]) << ", " << argv[i] << endl;
    else
      std::cout << sqrt(score[0]) << ", whole image index" << endl;
    return 0;
  } else if(m == 'c') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 2 + 1);
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      in.emplace_back(move(work));
      assert(in[0].size() == in[in.size() - 1].size() &&
             in[0][0].rows() == in[in.size() - 1][0].rows() &&
             in[0][0].cols() == in[in.size() - 1][0].cols() );
    }
    in = normalize<num_t>(zcollect<num_t>(in));
    for(int i = 0; i < in.size(); i ++)
      if(! savep2or3<num_t>((string(argv[i + 2]) + string("-c3.ppm")).c_str(), in[i]) )
        cerr << "failed to save." << endl;
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
  } else if(m == '?' || m == '!') {
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      num_t wavg(int(0));
      num_t wavgn0(int(0));
      num_t wavgn1(int(0));
      for(int i = 0; i < work.size(); i ++) {
        SimpleMatrix<complex(num_t) > lwork(dft<num_t>(work[i].rows()) *
          work[i].template cast<complex(num_t) >() *
            dft<num_t>(work[i].cols()).transpose() );
        work[i].resize(work[i].rows(), work[i].cols() * 2);
        for(int j = 0; j < work[i].rows(); j ++)
          for(int k = 0; k < lwork.cols(); k ++) {
            work[i](j, k) = abs(lwork(j, k));
            work[i](j, k + lwork.cols()) = arg(lwork(j, k));
          }
        SimpleMatrix<num_t> llwork(work[i].subMatrix(0, 0, work[i].rows(), lwork.cols() ));
        llwork.entity = normalizeS<num_t>(llwork.entity).first;
        work[i].setMatrix(0, 0, llwork);
        llwork = work[i].subMatrix(0, lwork.cols(), work[i].rows(), lwork.cols());
        llwork.entity = normalizeS<num_t>(llwork.entity).first;
        work[i].setMatrix(0, lwork.cols(), llwork);
        for(int j = 0; j < work[i].rows(); j ++)
          for(int k = 0; k < work[i].cols() / 2; k ++) {
            const num_t weight(sqrt(num_t(abs(j - work[i].rows() / 2) *
              abs(k - work[i].cols() / 4) ) /
                num_t(work[i].rows() / 2 * work[i].cols() / 4) ));
            wavg   += work[i](j, k) * weight;
            wavgn0 += work[i](j, k) * work[i](j, k);
            wavgn1 += weight * weight;
          }
        work[i].entity = offsetHalf<num_t>(work[i].entity);
      }
      std::cout << wavg / sqrt(wavgn0 * wavgn1) << endl;
      if(! savep2or3<num_t>((string(argv[i0]) + string("-ex.ppm")).c_str(),
        work) )  cerr << "failed to save." << endl;
    }
  } else if(m == 'L') {
    vector<SimpleMatrix<num_t> > in, inc;
    if(! loadp2or3<num_t>(in, argv[2])) return - 1;
    if(! loadp2or3<num_t>(inc, argv[3])) return - 1;
    assert(in.size() == 1 && inc.size() == 3);
    vector<SimpleMatrix<num_t> > out(inc);
    for(int i = 0; i < out[0].rows(); i ++)
      for(int j = 0; j < out[0].cols(); j ++) {
        num_t cavg(int(0));
        for(int k = 0; k < inc.size(); k ++) cavg += inc[k](i, j);
        cavg /= num_t(int(3));
        for(int k = 0; k < inc.size(); k ++) inc[k](i, j) *= in[0](i, j) / cavg;
      }
    if(! savep2or3<num_t>((string(argv[2]) + string("-color.ppm")
      ).c_str(), normalize<num_t>(out)) ) cerr << "failed to save." << endl;
  } else goto usage;
  cerr << "Done" << endl;
  lieonnStaticDestroy();
  return 0;
 usage:
  lieonnStaticDestroy();
  cerr << "Usage:" << endl;
  cerr << "# copy enlarge structure" << endl;
  cerr << argv[0] << " +<ratio> <pixels> <in0.ppm> ... > cache.txt" << endl;
  cerr << "# apply enlarge structure" << endl;
  cerr << argv[0] << " -<ratio> <in0.ppm> ... < cache.txt" << endl;
  cerr << "# predict following image" << endl;
  cerr << argv[0] << " p <in0.ppm> ..." << endl;
  cerr << "# predict down scanlines" << endl;
  cerr << argv[0] << " q <in0out.ppm> ..." << endl;
  cerr << "# show continuity" << endl;
  cerr << argv[0] << " [xyit] <in0.ppm> ..." << endl;
  cerr << "# some of the volume curvature like transform" << endl;
  cerr << argv[0] << " c <in0.ppm> ..." << endl;
  return - 1;
}

