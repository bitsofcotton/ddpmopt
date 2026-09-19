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

using std::cout;
using std::cerr;
using std::endl;
using std::atoi;
using std::string;
using std::vector;
using std::sort;
using std::binary_search;
using std::make_pair;
using std::istringstream;

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
    std::string s;
    for(int len = 0; 0 <= len && std::getline(std::cin, s, '\n'); ) {
      if(s[0] == '*') {
        L.emplace_back(vector<SimpleVector<num_t> >());
        len ++;
        continue;
      } else if(! s.size()) continue;
      SimpleVector<num_t> l;
      std::stringstream ins(s);
      ins >> l;
      L[len - 1].emplace_back(l / sqrt(l.dot(l)));
      l /= - l[0];
      l[0] = num_t(int(0));
      L[len - 1].emplace_back(l);
    }
    for(int i0 = 2; i0 < argc; i0 ++) {
      cerr << i0 - 2 << " / " << argc - 2 << endl;
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i0])) return - 1;
      if(argv[1][1] == '\0' && in.size() != 3) {
        std::cerr << argv[i0] << " doesn't include 3 colors" << std::endl;
        continue;
      }
      vector<SimpleMatrix<num_t> > out;
      if(argv[1][1] == '\0' || !std::atoi(&argv[1][1])) {
        out.emplace_back(in[0]);
        out[0].O();
      } else out.resize(in.size());
      for(int i = 0; i < in[0].rows(); i ++)
        for(int j = 0; j < in[0].cols(); j ++) if(argv[1][1] == '\0') {
          SimpleVector<num_t> work(4);
          for(int m = 0; m < in.size(); m ++)
            work[m + 1] = in[m](i, j);
          work[0] = num_t(int(1)) / num_t(int(2));
          SimpleVector<num_t> work2(makeProgramInvariant<num_t>(work).first);
          int idx(0);
          for(int m = 2; m < L[0].size(); m += 2)
            if(abs(L[0][idx].dot(work2)) <= abs(L[0][m].dot(work2))) idx = m;
          pair<SimpleVector<num_t>, num_t> vdp(makeProgramInvariant<num_t>(
            work, num_t(int(1)) ));
          vdp.first[0] = - L[0][idx + 1].dot(vdp.first) * sgn<num_t>(L[0][idx
            ].dot(vdp.first));
          out[0](i, j) = revertProgramInvariant(make_pair(vdp.first[0],
            vdp.second));
        } else {
          const int ratio(std::atoi(&argv[1][1]));
          int sq(sqrt(num_t(L.size() )));
          if((sq + 1) * (sq + 1) == L.size()) sq ++;
          for(int i = 0; i < in.size(); i ++) {
            out[i].resize(in[i].rows() * ratio, in[i].cols() * ratio);
            out[i].O();
            SimpleMatrix<int> cnt(out[i].rows(), out[i].cols());
            cnt.O();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
            for(int i1 = 0; i1 < in[i].rows(); i1 += sq / ratio)
              for(int j1 = 0; j1 < in[i].cols(); j1 += sq / ratio) {
                i1 = min(i1, int(in[i].rows() - sq / ratio));
                j1 = min(j1, int(in[i].cols() - sq / ratio));
                SimpleMatrix<num_t> w(in[i].subMatrix(i1, j1, sq / ratio, sq / ratio));
                SimpleMatrix<num_t> Q(w.QR());
                SimpleMatrix<num_t> R(Q * w);
                num_t MM(int(0));
                for(int k1 = 0; k1 < R.rows(); k1 ++)
                  for(int m1 = 0; m1 < R.cols(); m1 ++)
                    MM = max(MM, abs(R(k1, m1)));
                R /= MM;
                SimpleVector<num_t> wv(Q.rows() * Q.cols() +
                  (R.rows() + 1) * (R.cols() + 2) / 2 + 1);
                wv.O();
                for(int k1 = 0; k1 < Q.rows(); k1 ++)
                  wv.setVector(1 + k1 * Q.cols(), Q.row(k1));
                for(int k1 = 0, m1 = 0; k1 < R.rows(); k1 ++) {
                  wv.setVector(1 + Q.rows() * Q.cols() + m1, R.row(k1
                    ).subVector(k1, R.cols() - k1) );
                  m1 += R.cols() - k1;
                }
                wv = offsetHalf<num_t>(wv);
                pair<SimpleVector<num_t>, num_t> vdp(makeProgramInvariant<
                  num_t>(wv) );
                for(int k1 = 0; k1 < L.size(); k1 ++) {
                  if(cnt(i1 * ratio + k1 / sq, j1 * ratio + k1 % sq)) continue;
                  int idx(0);
                  for(int m = 2; m < L[k1].size(); m += 2)
                    if(abs(L[k1][idx].dot(vdp.first)) <=
                      abs(L[k1][m].dot(vdp.first))) idx = m;
                  SimpleVector<num_t> wwv(vdp.first);
                  wwv[0] = - L[k1][idx + 1].dot(vdp.first) * sgn<num_t>(
                    L[k1][idx].dot(vdp.first));
                  out[i](i1 * ratio + k1 / sq, j1 * ratio + k1 % sq) +=
                    revertProgramInvariant(make_pair(wwv, vdp.second))[0];
                  cnt(i1 * ratio + k1 / sq, j1 * ratio + k1 % sq) ++;
                }
              }
          }
          out = normalize<num_t>(out);
        }
      if(! savep2or3<num_t>((std::string(argv[i0]) + std::string(".pgm")).c_str(), out) )
        cerr << "failed to save." << endl;
    }
  } else if(m == '+') {
    vector<vector<SimpleVector<num_t> > > v;
    if(argv[1][1] == '\0' || !std::atoi(&argv[1][1])) {
      vector<vector<SimpleMatrix<num_t> > > in;
      vector<SimpleMatrix<num_t> > out;
      assert(! ((argc - 2) & 1));
      in.resize((argc - 2) / 2);
      out.resize((argc - 2) / 2);
      int cnt(0);
      for(int i = 2; i < argc; i ++) {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i])) continue;
        if(! (i & 1)) {
          assert(work.size() == 1);
          out[i / 2 - 1] = move(work[0]);
          cnt += out[i / 2 - 1].rows() * out[i / 2 - 1].cols();
        } else {
          in[i / 2 - 1] = move(work);
          assert(in[i / 2 - 1].size() == 3);
          assert(out[i / 2 - 1].rows() == in[i / 2 - 1][0].rows() &&
                 out[i / 2 - 1].cols() == in[i / 2 - 1][0].cols());
        }
      }
      assert(in.size() == out.size());
      v.resize(1);
      v[0].reserve(cnt);
      for(int i = 0; i < in.size(); i ++)
        for(int j = 0; j < out[i].rows(); j ++)
          for(int k = 0; k < out[i].cols(); k ++) {
            SimpleVector<num_t> work(4);
            for(int m = 0; m < 3; m ++) work[m + 1] = in[i][m](j, k);
            work[0] = out[i](j, k);
            v[0].emplace_back(move(work));
          }
    } else {
      const int ratio(std::atoi(&argv[1][1]));
      v.resize(std::atoi(argv[2]) * std::atoi(argv[2]));
      for(int i = 3; i < argc; i ++) {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i])) continue;
        // N.B. we bet shrinked image is near distance when pixel offset
        //      changed. however this isn't suppose clustering counts.
        for(int j = 0; j < work.size(); j ++)
         while(std::atoi(argv[2]) <= work[j].rows() &&
           std::atoi(argv[2]) <= work[j].cols()) {
          for(int i1 = 0; i1 <= work[j].rows() - std::atoi(argv[2]); i1 ++)
            for(int j1 = 0; j1 <= work[j].cols() - std::atoi(argv[2]); j1 ++) {
              SimpleMatrix<num_t> ww(work[j].subMatrix(i1, j1,
                std::atoi(argv[2]), std::atoi(argv[2]) ));
              SimpleMatrix<num_t> sw(ww.rows() / ratio, ww.cols() / ratio);
              sw.O();
              for(int k1 = 0; k1 < sw.rows(); k1 ++)
                for(int m1 = 0; m1 < sw.cols(); m1 ++)
                  for(int k2 = 0; k2 < ratio; k2 ++) for(int m2 = 0; m2 < ratio; m2 ++)
                    sw(k1, m1) += ww(k1 * ratio + k2, m1 * ratio + m2);
              sw /= num_t(ratio * ratio);
              SimpleMatrix<num_t> Q(sw.QR());
              SimpleMatrix<num_t> R(Q * sw);
              SimpleVector<num_t> wv(Q.rows() * Q.cols() +
                (R.rows() + 1) * (R.cols() + 2) / 2 + 1);
              num_t MM(int(0));
              for(int k1 = 0; k1 < R.rows(); k1 ++)
                for(int m1 = 0; m1 < R.cols(); m1 ++)
                  MM = max(MM, abs(R(k1, m1)));
              R /= MM;
              wv.O();
              for(int k1 = 0; k1 < Q.rows(); k1 ++)
                wv.setVector(1 + k1 * Q.cols(), Q.row(k1));
              for(int k1 = 0, m1 = 0; k1 < R.rows(); k1 ++) {
                wv.setVector(1 + Q.rows() * Q.cols() + m1, R.row(k1).subVector(
                  k1, R.cols() - k1) );
                m1 += R.cols() - k1;
              }
              for(int k1 = 0; k1 < ww.rows(); k1 ++)
                for(int m1 = 0; m1 < ww.cols(); m1 ++) {
                  wv[0] = ww(k1, m1);
                  v[k1 * ww.cols() + m1].emplace_back(offsetHalf<num_t>(wv));
                }
            }
          SimpleMatrix<num_t> w2(work[j].rows() / 2, work[j].cols() / 2);
          w2.O();
          for(int i0 = 0; i0 < w2.rows(); i0 ++)
            for(int j0 = 0; j0 < w2.cols(); j0 ++)
              for(int i1 = 0; i1 < 2; i1 ++) for(int j1 = 0; j1 < 2; j1 ++)
                w2(i0, j0) += work[j](min(i0 * 2 + i1, int(work[j].rows() - 1)),
                  min(j0 + j1, int(work[j].cols() - 1) ));
          work[j] = move(w2 /= num_t(int(4)));
         }
      }
    }
    for(int i0 = 0; i0 < v.size(); i0 ++) {
      vector<pair<vector<SimpleVector<num_t> >, vector<int> > > c(
        crush<num_t, true>(v[i0]));
      cout << "*" << endl;
      for(int i = 0; i < c.size(); i ++) {
        if(! c[i].first.size()) continue;
        SimpleVector<num_t> vv(makeProgramInvariant<num_t>(c[i].first[0]).first);
        for(int j = 1; j < c[i].first.size(); j ++)
          vv += makeProgramInvariant<num_t>(c[i].first[j]).first;
        vv /= num_t(c[i].first.size());
        if(vv.dot(vv) != num_t(int(0))) cout << vv;
      }
      cout << endl;
    }
  } else if(m == 'p' || m == 'T') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 1);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in.emplace_back(move(work));
    }
    vector<SimpleMatrix<num_t> > p(m == 'T' ? predMatTangleLast<num_t, 20>(
      move(in), 3, string(" ") + string(argv[0]) + string(" ") +
        string(argv[1]) ) : predMat<num_t, 20>(move(in), 3, string(" ") +
          string(argv[0]) + string(" ") + string(argv[1]) ) );
    if(! savep2or3<num_t>("predg.ppm", m == 'T' ? p : normalize<num_t>(p) ) )
      cerr << "failed to save." << endl;
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
        cout << sqrt(score[i]) << ", " << argv[i] << endl;
    else
      cout << sqrt(score[0]) << ", whole image index" << endl;
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
        std::cout << "const int cg_width = " << work[0].cols() << ";" << std::endl;
        std::cout << "const int cg_height = " << work[0].rows() << ";" << std::endl;
        std::cout << "const int cg_n = " << argc - 2 << ";" << std::endl;
        std::cout << "const char cg[" << argc - 2 << "][" << work[0].rows() * work[0].cols() << "] = {" << std::endl;
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
    int shift(std::atoi(argv[2]));
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
    const int len(std::atoi(argv[2]));
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
      std::vector<SimpleMatrix<num_t> > work;
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
      std::cout << wavg / sqrt(wavgn0 * wavgn1) << std::endl;
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
    if(! savep2or3<num_t>((std::string(argv[2]) + std::string("-color.ppm")
      ).c_str(), normalize<num_t>(out)) ) cerr << "failed to save." << endl;
  } else goto usage;
  cerr << "Done" << endl;
  lieonnStaticDestroy();
  return 0;
 usage:
  lieonnStaticDestroy();
  cerr << "Usage:" << endl;
  cerr << "# copy color structure" << endl;
  cerr << argv[0] << " + <in0out.pgm> <in0in.ppm> ... > cache.txt" << endl;
  cerr << "# apply color structure" << endl;
  cerr << argv[0] << " - <in0.ppm> ... < cache.txt" << endl;
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

