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
#include <cctype>
#include <assert.h>
#include <stdlib.h>

#define int int32_t
//#define int int64_t
#include "lieonn.hh"
typedef myfloat num_t;

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

#undef int
int main(int argc, const char* argv[]) {
#define int int32_t
//#define int int64_t
  const int   sz(2);
  const char& m(argv[1][0]);
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
          SimpleVector<num_t> work2(makeProgramInvariant<num_t>(work).first);
          assert(work2.size() == L[0].size());
          int idx(0);
          for(int m = 2; m < L.size(); m += 2) {
            assert(L[m].size()   == work2.size());
            assert(L[idx].size() == work2.size());
            if(abs(L[idx].dot(work2)) <= abs(L[m].dot(work2)))
              idx = m;
          }
          num_t last(sqrt(work.dot(work)));
          for(int ii = 0;
                  ii < 2 * int(- log(SimpleMatrix<num_t>().epsilon()) / log(num_t(int(2))) )
                  && sqrt(work.dot(work) * SimpleMatrix<num_t>().epsilon()) <
                       abs(work[3] - last); ii ++) {
            last = work[3];
            const pair<SimpleVector<num_t>, num_t> work2(makeProgramInvariant<num_t>(work));
            work[3] = revertProgramInvariant<num_t>(make_pair(L[idx + 1].dot(work2.first) * sgn<num_t>(L[idx].dot(work2.first)), work2.second) );
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
    const vector<pair<vector<SimpleVector<num_t> >, vector<int> > > c(crush<num_t>(v));
    for(int i = 0; i < c.size(); i ++) {
      if(! c[i].first.size()) continue;
      SimpleVector<num_t> vv(makeProgramInvariant<num_t>(c[i].first[0]).first);
      for(int j = 1; j < c[i].first.size(); j ++)
        vv += makeProgramInvariant<num_t>(c[i].first[j]).first;
      vv /= num_t(c[i].first.size());
      if(vv.dot(vv) != num_t(int(0))) cout << vv;
    }
    cout << endl;
  } else if(m == 'p') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 1);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
    }
    // N.B. with good spreaded input, we can suppose original as a 'T' command
    //      case.
    vector<vector<SimpleMatrix<num_t> > > p(predMat<num_t>(in = normalize<num_t>(in)));
    for(int i = 0; i < p.size(); i ++)
      if(! savep2or3<num_t>(
        (string("predg") + to_string(i) + string(".ppm")).c_str(),
        p[i].size() == 3 ? xyz2rgb<num_t>(p[i]) : move(p[i]) ))
          cerr << "failed to save." << endl;
  } else if(m == 'P') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 1);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
    }
    // N.B. 10 + 1 * 2 < work.size() / step for predv[pq].
    for(int j = 1; j < in.size() / 12; j ++) {
      SimpleVector<vector<SimpleMatrix<num_t> > > work;
      work.entity = skipX<vector<SimpleMatrix<num_t> > >(in, j);
      vector<vector<SimpleMatrix<num_t> > > p(
        predMat<num_t>(work.entity = normalize<num_t>(work.subVector(work.size() - 12, 12).entity)));
      for(int i = 0; i < p.size(); i ++)
        if(! savep2or3<num_t>(
          (string("predg") + to_string(j) + string("-") + to_string(i)
            + string(".ppm")).c_str(), 
          p[i].size() == 3 ? xyz2rgb<num_t>(p[i]) : move(p[i]) ))
          cerr << "failed to save." << endl;
    }
  } else if(m == 'w') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 2);
    for(int i = 2; i < argc; i ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i])) continue;
      in.emplace_back(work.size() == 3 ? rgb2xyz<num_t>(work) : move(work));
    }
    in = normalize<num_t>(in);
    vector<SimpleVector<num_t> > work;
    work.resize(in.size());
    for(int i = 0; i < in.size(); i ++) {
      work[i].resize(in[i].size() * in[i][0].rows() * in[i][0].cols());
      for(int j = 0; j < in[i].size(); j ++)
        for(int k = 0; k < in[i][j].rows(); k ++)
          work[i].setVector(j * in[i][0].rows() * in[i][0].cols() +
            k * in[i][0].cols(), in[i][j].row(k));
    }
    SimpleVector<num_t> vp(predv4<num_t, 20>(work));
    vector<SimpleMatrix<num_t> > p;
    p.resize(in[1].size());
    for(int i = 0; i < p.size(); i ++) {
      p[i].resize(in[1][0].rows(), in[1][0].cols());
      for(int j = 0; j < p[i].rows(); j ++)
        p[i].row(j) = vp.subVector(i * p[0].rows() * p[0].cols() +
          j * p[0].cols(), p[0].cols());
    }
    if(! savep2or3<num_t>("predgw.ppm",
      normalize<num_t>(p.size() == 3 ? xyz2rgb<num_t>(p) : move(p))) )
        cerr << "failed to save." << endl;
  } else if(m == 'q' || m == 'Q') {
    for(int i0 = 2; i0 < argc; i0 ++) {
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
      // N.B. same as 'p' cmd, we can suppose original as 'T' command input
      //      with long range but not in general.
      // N.B. 10 + 1 * 2 < work[0].rows() / step for PP0.
      const int ext(work[0].rows() / 12);
      vector<vector<SimpleMatrix<num_t> > > wwork;
      for(int i = 0; i < ext; i ++) {
        SimpleVector<vector<SimpleVector<num_t> > > w;
        w.entity = skipX<vector<SimpleVector<num_t> > >(pwork, i + 1);
        vector<vector<SimpleVector<num_t> > > n(
          predVec<num_t>(m == 'q' ? w.entity : w.subVector(w.size() - 12, 12).entity));
        if(! i) {
          wwork.resize(n.size());
          for(int k = 0; k < n.size(); k ++) {
            wwork[k].resize(work.size(),
              SimpleMatrix<num_t>(work[0].rows() + ext, work[0].cols()).O());
            for(int j = 0; j < work.size(); j ++)
              wwork[k][j].setMatrix(0, 0, work[j]);
          }
        }
        for(int k = 0; k < n.size(); k ++)
          for(int j = 0; j < wwork[k].size(); j ++)
            wwork[k][j].row(work[0].rows() + i) = move(n[k][j]);
      }
      for(int i = 0; i < wwork.size(); i ++)
        if(! savep2or3<num_t>(
          (string(argv[i0]) + to_string(i) + string(".ppm")).c_str(),
          wwork[i].size() == 3 ? xyz2rgb<num_t>(wwork[i]) : move(wwork[i]) ) )
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
    cerr << "y" << flush;
    vector<vector<SimpleMatrix<num_t> > > jy(in);
    SimpleMatrix<num_t> left(diff<num_t>(jy[0][0].rows()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < in.size(); i ++) {
      cerr << "." << flush;
      for(int j = 0; j < in[i].size(); j ++)
        jy[i][j] = left * in[i][j];
    }
    cerr << "x" << flush;
    vector<vector<SimpleMatrix<num_t> > > jx(in);
    SimpleMatrix<num_t> right(diff<num_t>(jy[0][0].cols()).transpose());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
    for(int i = 0; i < in.size(); i ++) {
      cerr << "." << flush;
      for(int j = 0; j < in[i].size(); j ++)
        jx[i][j] = in[i][j] * right;
    }
    cerr << "z" << flush;
    vector<vector<SimpleMatrix<num_t> > > jz(in);
    SimpleMatrix<num_t> middle(diff<num_t>(jz.size()));
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
    in = normalize<num_t>(in);
    for(int i = 0; i < in.size(); i ++)
      if(! savep2or3<num_t>((string(argv[i + 2]) + string("-c3.ppm")).c_str(), in[i]) )
        cerr << "failed to save." << endl;
  } else if(m == 'T') {
    vector<vector<SimpleMatrix<num_t> > > in;
    in.reserve(argc - 2 + 1);
    vector<vector<SimpleMatrix<num_t> > > p;
    vector<std::pair<int, int> > pc;
    for(int i0 = 2; i0 < argc; i0 ++) {
      vector<SimpleMatrix<num_t> > work;
      if(! loadp2or3<num_t>(work, argv[i0])) continue;
      if(pc.size()) {
        cout << " --- " << in.size() - 11 << " --- " << endl;
        for(int i = 0; i < pc.size(); i ++)
          for(int j = 0; j < p[0].size(); j ++) {
            const int rr(p[0][j].rows() / pc[i].first);
            const int cc(p[0][j].cols() / pc[i].second);
            for(int k = 0; k < pc[i].first; k ++)
              for(int n = 0; n < pc[i].second; n ++) {
                vector<num_t> workr;
                num_t orig(int(0));
                int   cnt(0);
                workr.resize(p.size(), int(0));
                for(int kk = k * rr; kk < min(p[0][j].rows(), (k + 1) * rr); kk ++)
                  for(int nn = n * cc; nn < min(p[0][j].cols(), (n + 1) * cc);
                    nn ++, cnt ++) {
                    for(int m = 0; m < p.size(); m ++)
                      workr[m] += p[m][j](kk, nn);
                    orig += work[j](kk, nn);
                  }
                if(cnt) {
                  for(int m = 0; m < workr.size(); m ++) workr[m] /= num_t(cnt);
                  orig  /= num_t(cnt);
                }
                // N.B. for stricter test but we don't need such a restriction:
                // workr = (sgn<num_t>(workr - num_t(int(1)) / num_t(int(2)) )
                //   + num_t(int(1)) ) / num_t(int(2));
                // also apply this on orig.
                for(int m = 0; m < p.size() - 1; m ++)
                  cout << abs(workr[m] - orig) * num_t(int(2)) << ", ";
                cout << abs(workr[p.size() - 1] - orig) * num_t(int(2)) << endl;
              }
            cout << endl;
          }
        cout << endl;
        // N.B. output can be checked as:
        //      tail -n ... < output | python3 p2/cr.py t 1 > outR
        //      with R.app, myv <- read.csv("outR")
        //                  hist(yv[,1],breaks=seq(0,...,length.out=...))
      }
      in.emplace_back(move(work));
      if(i0 == argc - 1) break;
      if(11 < in.size()) {
        vector<vector<SimpleMatrix<num_t> > > in2(in);
        p = predMat<num_t>(in2 = normalize<num_t>(in2));
        if(! pc.size()) {
          pc.emplace_back(make_pair(p[0][0].rows(), p[0][0].cols() ));
          for(int i = 1;
            0 <= i && 1 < pc[i - 1].first && 1 < pc[i - 1].second; i ++) {
            pair<int, int> work(pc[i - 1]);
            work.first  /= 2;
            work.second /= 2;
            pc.emplace_back(move(work));
          }
        }
      }
    }
  } else goto usage;
  cerr << "Done" << endl;
  return 0;
 usage:
  cerr << "Usage:" << endl;
  cerr << "# copy color structure" << endl;
  cerr << argv[0] << " + <in0out.pgm> <in0in.ppm> ... > cache.txt" << endl;
  cerr << "# apply color structure" << endl;
  cerr << argv[0] << " - <in0.ppm> ... < cache.txt" << endl;
  cerr << "# predict following image (each bit input)" << endl;
  cerr << argv[0] << " [pP] <in0.ppm> ..." << endl;
  cerr << "# predict with whole pixel context (each bit input)" << endl;
  cerr << argv[0] << " w <in0-4.ppm> <in0.ppm> ... <addition-4.ppm>" << endl;
  cerr << "# predict down scanlines. (each bit input)" << endl;
  cerr << argv[0] << " [qQ] <in0out.ppm> ..." << endl;
  cerr << "# show continuity" << endl;
  cerr << argv[0] << " [xyit] <in0.ppm> ..." << endl;
  cerr << "# some of the volume curvature like transform" << endl;
  cerr << argv[0] << " c <in0.ppm> ..." << endl;
  cerr << "# test input series of graphics predictable or not (each bit input)" << endl;
  cerr << argv[0] << " T <in0.ppm> ..." << endl;
  return - 1;
}

