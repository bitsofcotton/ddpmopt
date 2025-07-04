#! /usr/bin/env python3
import sys
import io
import struct
import hashlib
import subprocess

if(len(sys.argv) < 2):
  print("not much arguments")
elif(sys.argv[1][0] == 'R'):
  while(True):
    a = []
    idx = 0
    for line in io.open(sys.stdin.fileno(), 'r', buffering = 1, encoding = "utf-8", closefd = False):
      a.extend(struct.pack('f', float(line.split(",")[0])))
      if(abs(int(sys.argv[2])) <= idx): break
      idx += 1
    m = hashlib.sha256()
    m.update(bytearray(a))
    a = m.digest()
    for byte in a:
      print(byte - 127.5)
    sys.stdout.flush()
elif(sys.argv[1][0] == 'm'):
  import numpy
  import mido
  from mido import Message, MidiFile, MidiTrack, MetaMessage
  tbl = []
  mid = MidiFile()
  tbl0 = [0, 2, 4, 5, 7, 9, 11]
  for u in range(0, 2 * len(tbl0)):
    # Thanks to google with search word midi C value, 60 is one of a C.
    tbl.append(60 + tbl0[u % len(tbl0)] + int(u / len(tbl0)) * 12)
  mC = []
  ctr = 0
  bw  = []
  track = []
  for line in sys.stdin:
    if(len(line.split("[")) <= 1): continue
    ff = line.split("[")[1].split("]")[0].split(",")
    if(len(track) <= 0):
      for s in range(0, len(ff)):
        # Thanks to : https://qiita.com/tjsurume/items/75a96381fd57d5350971 via search engine
        trk = MidiTrack()
        trk.append(MetaMessage('set_tempo', tempo=mido.bpm2tempo(70)))
        # cf.https://ja.wikipedia.org/wiki/General_MIDI
        #trk.append(Message('program_change', program = 117))
        track.append(trk)
        mid.tracks.append(trk)
    print(ff, bw)
    if(len(bw) <= 0):
      bw = ff
      continue
    flg = True
    for w in ff:
      if(abs(float(w)) > 1. / 20):
        flg = False
        break
    for w in bw:
      if(abs(float(w)) > 1. / 20):
        flg = False
        break
    if(flg): continue
    idx = 0
    bb  = []
    for w in range(0, len(ff)):
      if(len(bb) <= idx): bb.append(int(len(tbl) / 2))
      bb[idx] += abs(numpy.arctan(abs(numpy.tan(float(ff[w]) * 6))) * 6) - 3
      f = tbl[int(bb[idx]) % int(len(tbl) / len(ff)) + idx * int(len(tbl) / len(ff))]
      idx += 1
      track[w].append(Message('note_on', note=f, velocity=127, time=0))
      track[w].append(Message('note_off', note=f, time=(120 * int(abs(numpy.arctan(abs(numpy.tan(float(bw[w])))) * 3) + 1)) ))
    bw = []
  mid.save('rand_correct.mid')
  # XXX: following are duplicate for non unistd systems.
elif(sys.argv[1][0] == 'H'):
  p = []
  for line in io.open(sys.stdin.fileno(), 'r', buffering = 1, encoding = 'utf-8', closefd = False):
    d = line[:- 1].split(",")
    if(len(p) < len(d)):
      for t in range(len(p), len(d)):
        p.append(subprocess.Popen(["sh", "-c", sys.argv[2]], stdin = subprocess.PIPE, stdout = subprocess.PIPE))
    for t in range(0, len(d)):
      p[t].stdin.write((d[t] + "\n").encode("utf-8"))
      p[t].stdin.flush()
    J = []
    for t in range(0, len(d)):
      J.extend(p[t].stdout.readline().decode("utf-8")[:- 1].split(","))
    print(",".join(J))
    sys.stdout.flush()
elif(sys.argv[1][0] == 'D'):
  p = []
  for a in sys.argv[2:]:
    p.append(subprocess.Popen(["sh", "-c", a], stdin = subprocess.PIPE, stdout = subprocess.PIPE))
  for line in io.open(sys.stdin.fileno(), 'r', buffering = 1, encoding = 'utf-8', closefd = False):
    for pp in p:
      pp.stdin.write(line.encode("utf-8"))
      pp.stdin.flush()
    D = []
    for pp in p:
      D.extend(pp.stdout.readline().decode("utf-8")[:- 1].split(","))
    D.extend(line[:- 1].split(","))
    print(",".join(D))
    sys.stdout.flush()

