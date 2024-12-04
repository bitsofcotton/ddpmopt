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

#undef int
int main(int argc, const char* argv[]) {
//#define int int64_t
#define int int32_t
  const auto  sz(2);
  const auto& m(argv[1][0]);
  if(argc <= 1 || argv[1][1] != '\0') goto usage;
  cerr << "Coherent: sqrt(2): " << sqrt<num_t>(Complex<num_t>(num_t(2))) << endl;
  if(m == '-') {
    vector<SimpleVector<num_t> > L;
    std::string s;
    while(std::getline(std::cin, s, '\n')) {
      SimpleVector<num_t> l;
      std::stringstream ins(s);
      ins >> l;
      if(l.size() != sz * sz + 1) break;
      L.emplace_back(l / sqrt(l.dot(l)));
      l /= - num_t(l[3]);
      l[3] = num_t(int(0));
      L.emplace_back(l);
    }
    assert(L.size() && ! (L.size() & 1));
    for(int i0 = 2; i0 < argc; i0 ++) {
      cerr << i0 - 2 << " / " << argc - 2 << endl;
      vector<SimpleMatrix<num_t> > in;
      if(! loadp2or3<num_t>(in, argv[i0])) return - 1;
      if(in.size() != 3) {
        std::cerr << argv[i0] << " doesn't include 3 colors" << std::endl;
        continue;
      }
      vector<SimpleMatrix<num_t> > out;
      out.emplace_back(in[0]);
      out[0].O();
      for(int i = 0; i < in[0].rows(); i ++)
        for(int j = 0; j < in[0].cols(); j ++) {
          SimpleVector<num_t> work(4);
          for(int m = 0; m < in.size(); m ++)
            work[m] = in[m](i, j);
          work[3] = num_t(int(1)) / num_t(int(2));
          auto work2(makeProgramInvariant<num_t>(work, - num_t(int(1)), true).first);
          assert(work2.size() == L[0].size());
          int idx(0);
          for(int m = 2; m < L.size(); m += 2) {
            assert(L[m].size()   == work2.size());
            assert(L[idx].size() == work2.size());
            if(abs(L[idx].dot(work2)) <= abs(L[m].dot(work2)))
              idx = m;
          }
          auto last(sqrt(work.dot(work)));
          for(int ii = 0;
                  ii < 2 * int(- log(SimpleMatrix<num_t>().epsilon()) / log(num_t(int(2))) )
                  && sqrt(work.dot(work) * SimpleMatrix<num_t>().epsilon()) <
                       abs(work[3] - last); ii ++) {
            last = work[3];
            const auto work2(makeProgramInvariant<num_t>(work, - num_t(int(1)), true));
            work[3] = revertProgramInvariant<num_t>(make_pair(L[idx + 1].dot(work2.first) * sgn<num_t>(L[idx].dot(work2.first)), work2.second), true);
          }
          out[0](i, j) = work[3];
        }
      if(! savep2or3<num_t>((std::string(argv[i0]) + std::string(".pgm")).c_str(), out) )
        cerr << "failed to save." << endl;
    }
  } else if(m == '+') {
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
    vector<SimpleVector<num_t> > v;
    v.reserve(cnt);
    for(int i = 0; i < in.size(); i ++)
      for(int j = 0; j < out[i].rows(); j ++)
        for(int k = 0; k < out[i].cols(); k ++) {
          SimpleVector<num_t> work(4);
          for(int m = 0; m < 3; m ++)
            work[m] = in[i][m](j, k);
          work[3] = out[i](j, k);
          v.emplace_back(move(work));
        }
    const auto c(crush<num_t>(v));
    for(int i = 0; i < c.size(); i ++) {
      if(! c[i].first.size()) continue;
      auto vv(makeProgramInvariant<num_t>(c[i].first[0], - num_t(int(1)), true).first);
      for(int j = 1; j < c[i].first.size(); j ++)
        vv += makeProgramInvariant<num_t>(c[i].first[j], - num_t(int(1)), true).first;
      vv /= num_t(c[i].first.size());
      if(vv.dot(vv) != num_t(int(0))) cout << vv;
    }
    cout << endl;
  } else if(m == '0' || m == 'a') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 1);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
    }
    if(m == '0') {
      auto p(predMat<num_t>(in = normalize<num_t>(in), - 1, 1));
      assert(p.size() == 1);
      if(! savep2or3<num_t>("predg.ppm",
          normalize<num_t>(p[0].size() == 3 ? xyz2rgb<num_t>(p[0]) : p[0])) )
            cerr << "failed to save." << endl;
      in.reserve(argc - 1);
      for(int i = 2; i < argc; i ++) {
        vector<SimpleMatrix<num_t> > work;
        if(! loadp2or3<num_t>(work, argv[i])) continue;
        in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
      }
    }
    {
      auto p(m == 'a' ?
        predMat<num_t>(in = normalize<num_t>(in), 0, 0) :
        predMat<num_t>(in = normalize<num_t>(in)));
      for(int i = 0; i < p.size(); i ++) {
        if(! savep2or3<num_t>((string("predg-") + to_string(i) +
          string(".ppm")).c_str(),
            normalize<num_t>(p[i].size() == 3 ? xyz2rgb<num_t>(p[i]) : p[i])) )
              cerr << "failed to save." << endl;
      }
    }
  } else if(m == 'w') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 1);
    int wp(0);
    int c(0);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      wp += work[0].rows() * work[0].cols();
      c   = max(c, int(work.size()));
      in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
    }
    vector<vector<SimpleVector<num_t> > > ep;
    ep.resize(wp);
    int pp(0);
    for(int i = 0; i < in.size(); i ++)
      for(int j = 0; j < in[i][0].rows(); j ++)
        for(int k = 0; k < in[i][0].cols(); k ++, pp ++) {
          ep[pp].resize(1, SimpleVector<num_t>(c * 2).O());
          for(int m = 0; m < c; m ++)
            ep[pp][0][m] = in[i][m % in[i].size()](j, k);
        }
    assert(pp == wp);
    for(pp = 0; pp < ep.size(); pp ++)
      for(int m = 0; m < c; m ++)
        ep[ep.size() - pp - 1][0][c + m] = ep[pp][0][m];
    assert(pp == wp);
    auto p(predVec<num_t>(ep = normalize<num_t>(ep)));
    assert(in[0][0].rows() * in[0][0].cols() <= p.size());
    for(int i = 0, pp = 0; i < p.size() / (in[0][0].rows() * in[0][0].cols()); i ++) {
      vector<SimpleMatrix<num_t> > outf(c);
      for(int j = 0; j < outf.size(); j ++) {
        outf[j].resize(in[0][0].rows(), in[0][0].cols());
        outf[j].O();
      }
      for(int ii = 0; ii < outf[0].rows(); ii ++)
        for(int jj = 0; jj < outf[0].cols() && pp < p.size(); jj ++, pp ++)
          for(int j = 0; j < outf.size(); j ++)
            outf[j](ii, jj) = p[pp][0][j];
      if(! savep2or3<num_t>((string("predgw-") + to_string(i) + string(".ppm")).c_str(),
          normalize<num_t>(outf.size() == 3 ? xyz2rgb<num_t>(outf) : outf)) )
        cerr << "failed to save." << endl;
      vector<SimpleMatrix<num_t> > outwf(c);
      for(int j = 0; j < outwf.size(); j ++) {
        outwf[j].resize(1, 4);
        outwf[j].O();
      }
      for(int ii = 0; ii < 4 && pp < p.size(); ii ++, pp ++)
        for(int j = 0; j < outwf.size(); j ++)
          outwf[j](0, ii) = p[pp][0][j];
      if(! savep2or3<num_t>((string("predgw-4-") + to_string(i) + string(".ppm")).c_str(),
          normalize<num_t>(outwf.size() == 3 ? xyz2rgb<num_t>(outwf) : outwf)) )
        cerr << "failed to save." << endl;
    }
  } else if(m == 'q') {
    for(int i0 = 1; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      work = normalize<num_t>(work.size() == 3 ? rgb2xyz<num_t>(work) : work);
      vector<vector<SimpleVector<num_t> > > pwork;
      pwork.resize(work[0].rows());
      for(int i = 0; i < pwork.size(); i ++) {
        pwork[i].reserve(work.size());
        for(int j = 0; j < work.size(); j ++)
          pwork[i].emplace_back(work[j].row(i));
      }
      auto p(predVec<num_t>(pwork));
      vector<SimpleMatrix<num_t> > wwork(work.size(),
        SimpleMatrix<num_t>(work[0].rows() + p.size(), work[0].cols()).O());
      for(int j = 0; j < work.size(); j ++)
        wwork[j].setMatrix(0, 0, work[j]);
      for(int i = 0; i < p.size(); i ++)
        for(int j = 0; j < p[i].size(); j ++)
          wwork[j].row(i - p.size() + wwork[j].rows()) = move(p[i][j]);
      if(! savep2or3<num_t>(argv[i0], normalize<num_t>(wwork.size() == 3 ?
        xyz2rgb<num_t>(wwork) : wwork) ) )
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
        P0maxRank<num_t> p;
        if(! loadp2or3<num_t>(work, argv[i0])) continue;
        for(int i = 0; i < work.size(); i ++)
          for(int ii = 0; ii < work[i].rows(); ii ++) {
            idFeeder<num_t> w(3);
            for(int jj = 0; jj < work[i].cols(); jj ++) {
              if(w.full) {
                const auto pp(p.next(w.res));
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
        P0maxRank<num_t> p;
        if(! loadp2or3<num_t>(work, argv[i0])) continue;
        for(int i = 0; i < work.size(); i ++)
          for(int jj = 0; jj < work[i].cols(); jj ++) {
            idFeeder<num_t> w(3);
            for(int ii = 0; ii < work[i].rows(); ii ++) {
              if(w.full) {
                const auto pp(p.next(w.res));
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
        P0maxRank<num_t> p;
        int cnt(0);
        for(int i = 0; i < b[2][0].rows(); i ++)
          for(int j = 0; j < b[2][0].cols(); j ++)
            for(int k = 0; k < b[2].size(); k ++) {
              idFeeder<num_t> w(3);
              for(int ii = 2; ii < b.size(); ii ++, cnt ++) {
                if(! b[ii].size()) continue;
                if(w.full) {
                  const auto pp(p.next(w.res));
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
    cerr << "y" << flush;
    auto jy(in);
    auto left(diff<num_t>(jy[0][0].rows()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < in.size(); i ++) {
      cerr << "." << flush;
      for(int j = 0; j < in[i].size(); j ++)
        jy[i][j] = left * in[i][j];
    }
    cerr << "x" << flush;
    auto jx(in);
    auto right(diff<num_t>(jy[0][0].cols()).transpose());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < in.size(); i ++) {
      cerr << "." << flush;
      for(int j = 0; j < in[i].size(); j ++)
        jx[i][j] = in[i][j] * right;
    }
    cerr << "z" << flush;
    auto jz(in);
    auto middle(diff<num_t>(jz.size()));
    for(int i = 0; i < in[0].size(); i ++) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
      for(int j = 0; j < in[0][0].rows(); j ++) {
        cerr << "." << flush;
        for(int k = 0; k < in[0][0].cols(); k ++) {
          SimpleVector<num_t> work(in.size());
          for(int m = 0; m < work.size(); m ++)
            work[m] = in[m][i](j, k);
          work = middle * work;
          for(int m = 0; m < work.size(); m ++)
            jz[m][i](j, k) = work[m];
        }
      }
    }
    cerr << endl << "eigen: " << flush;
    // [[1, 0, 0, jx], [0, 1, 0, jy], [0, 0, 1, jz], [jx, jy, jz, z]]
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < in.size(); i ++) {
      cerr << "." << flush;
      for(int j = 0; j < in[i].size(); j ++)
        for(int k = 0; k < in[i][j].rows(); k ++)
          for(int m = 0; m < in[i][j].cols(); m ++) {
            SimpleMatrix<num_t> work(4, 4);
            work.I();
            work(3, 0) = work(0, 3) = jx[i][j](k, m);
            work(3, 1) = work(1, 3) = jy[i][j](k, m);
            work(3, 2) = work(2, 3) = jz[i][j](k, m);
            work(3, 3) = work(3, 3) = in[i][j](k, m);
            in[i][j](k, m) = work.determinant();
          }
    }
    for(int i = 0; i < in.size(); i ++)
      if(! savep2or3<num_t>((string(argv[i + 2]) + string("-c3.ppm")).c_str(), normalize<num_t>(in[i]) ) )
        cerr << "failed to save." << endl;
  } else goto usage;
  cerr << "Done" << endl;
  return 0;
 usage:
  cerr << "Usage:" << endl;
  cerr << "# copy color structure" << endl;
  cerr << argv[0] << " + <in0out.pgm> <in0in.ppm> ... > cache.txt" << endl;
  cerr << "# apply color structure" << endl;
  cerr << argv[0] << " - <in0.ppm> ... < cache.txt" << endl;
  cerr << "# predict next image mode === '0' for normal, mode == 'a' to get all." << endl;
  cerr << argv[0] << " [0a] <in0.ppm> ..." << endl;
  cerr << "# predict with whole pixel context" << endl;
  cerr << argv[0] << " w <in0.ppm> <in0-4.ppm> ..." << endl;
  cerr << "# predict down scanlines." << endl;
  cerr << argv[0] << " q <in0out.ppm> ..." << endl;
  cerr << "# show continuity" << endl;
  cerr << argv[0] << " [xyit] <in0.ppm> ..." << endl;
  cerr << "# some of the volume curvature like transform" << endl;
  cerr << argv[0] << " c <in0.ppm> ..." << endl;
  return - 1;
}

