import sys
import numpy
from PIL import Image

def invariant(m):
  Q, R = numpy.linalg.qr(m, mode="reduced")
  chk = []
  lastQ = []
  for t in range(0, Q.shape[1] - 1):
    vvv = Q[0, 0:Q.shape[1]]
    for qq in range(1, Q.shape[0]):
      vvv += Q[qq, 0:Q.shape[1]]
    vvv  = numpy.reshape(vvv, (Q.shape[1], 1))
    mone = numpy.dot(Q, vvv)
    midx = 0
    nchk = True
    for qq in range(0, Q.shape[0]):
      flag = False
      for j in chk:
        if(chk == qq): flag = True
      if(flag): continue
      if(nchk or abs(mone[qq]) < abs(mone[midx])):
        midx = qq
        nchk = False
    chk.append(midx)
    lastQ = []
    for tt in range(0, Q.shape[0]):
      for ttt in range(0, Q.shape[1]):
        lastQ.append(Q[tt, ttt])
    lastQ = numpy.reshape(lastQ, Q.shape)
    mq = Q[midx, :]
    nbrk = False
    for qq in range(0, Q.shape[0]):
      dt  = 0.
      dd  = 0.
      for qqq in range(0, Q.shape[1]):
        dt += Q[qq, qqq] * mq[qqq]
        dd += mq[qqq] * mq[qqq]
      for qqq in range(0, Q.shape[1]):
        Q[qq, qqq] = Q[qq, qqq] - dt / dd * mq[qqq]
        if(not numpy.isfinite(Q[qq, qqq])):
          Q = lastQ
          nbrk = True
        if(nbrk): break
      if(nbrk): break
    if(nbrk): break
  vvv = Q[0, 0:Q.shape[1]]
  for qq in range(1, Q.shape[0]):
    vvv += Q[qq, 0:Q.shape[1]]
  #print(vvv)
  vvv  = numpy.reshape(vvv, (Q.shape[1], 1))
  try:
    return numpy.linalg.solve(R, vvv)
  except:
    pass
  vvv = [1, 1, 1, 1]
  vvv = numpy.reshape(vvv, (4, 1))
  return vvv
def p01n(v, skip):
  vv = v[- len(v) + (len(v) % skip):]
  ivts = []
  for i0 in range(0, 3):
    m = []
    for i in range(i0, i0 + int((len(vv) - 4 - skip + 2 - 1) / skip) - 2):
      for j in range(i * skip, i * skip + 4):
        m.append(numpy.float32(vv[j]))
      m[- 1] = numpy.float32(vv[i * skip + 4 + skip - 2])
    m = numpy.array(m, dtype=numpy.float32)
    m = numpy.resize(m, (int((len(vv) - 4 - skip + 2 - 1) / skip) - 2, 4))
    ivts.append(invariant(m))
  ivt = ivts[0] * .302831 + ivts[1] * .302831 + ivts[2] * .394338
  return - (vv[- 3] * ivt[0,0] + vv[- 2] * ivt[1,0] + vv[- 1] * ivt[2,0]) / ivt[3,0]
imgvec = []
imgsz  = ()
t = 0
for line in sys.argv[1:]:
  im = Image.open(line)
  if(len(imgsz) == 0):
    imgsz = im.size
  imgvec.append([])
  for y in range(0, imgsz[0]):
    for x in range(0, imgsz[1]):
      igp = im.getpixel((y, x))
      for u in range(0, len(igp)):
        imgvec[t].append(igp[u])
  t += 1
imgnorm = []
for t in range(0, len(imgvec)):
  res = 0.
  for s in range(0, len(imgvec[t])):
    if(imgvec[t][s] <= 1e-4): imgvec[t][s] = 1.
    res += numpy.log(imgvec[t][s])
  res /= len(imgvec[t])
  res  = numpy.exp(res)
  imgnorm.append(res)
rimgnorm = []
for t in range(0, len(imgnorm)):
  rimgnorm.append(imgnorm[- t - 1])
pn = []
qn = []
for k in range(0, 1):
  pn.append([])
  qn.append([])
  for j in range(0, len(imgvec[0])):
    pn[- 1].append(0.)
    qn[- 1].append(0.)
for j in range(0, len(imgvec[0])):
  print(j, " / ", len(imgvec[0]))
  nbuf = []
  for i in range(0, len(imgvec)):
    nbuf.append(imgvec[i][j] / imgnorm[i])
  bbuf = []
  for i in range(0, len(nbuf)):
    bbuf.append(nbuf[- i - 1])
  pn[0][j] = p01n(nbuf, 1)
  qn[0][j] = p01n(bbuf, 1)
  if(not numpy.isfinite(pn[0][j])): pn[0][j] = 0.
  if(not numpy.isfinite(qn[0][j])): qn[0][j] = 0.
ps = []
qs = []
for i in range(0, len(pn)):
  ps.append(p01n(imgnorm, i + 1))
  qs.append(p01n(rimgnorm, i + 1))
def floor1(x):
  return abs(x) - numpy.floor(abs(x))
for pp in range(0, len(pn)):
  imp = Image.new("RGB", imgsz, (255, 255, 255))
  imq = Image.new("RGB", imgsz, (255, 255, 255))
  for t in range(0, int(len(pn[pp]) / 3)):
    imp.putpixel((int(t / imgsz[1]), int(t % imgsz[1])), \
      (int(floor1(pn[pp][3 * t] * ps[pp]) * 255), \
       int(floor1(pn[pp][3 * t + 1] * ps[pp]) * 255), \
       int(floor1(pn[pp][3 * t + 2] * ps[pp]) * 255)) )
    imq.putpixel((int(t / imgsz[1]), int(t % imgsz[1])), \
      (int(floor1(qn[pp][3 * t] * qs[pp]) * 255), \
       int(floor1(qn[pp][3 * t + 1] * qs[pp]) * 255), \
       int(floor1(qn[pp][3 * t + 2] * qs[pp]) * 255)) )
  imp.save("./predf-" + str(pp) + ".png")
  imq.save("./predb-" + str(pp) + ".png")

