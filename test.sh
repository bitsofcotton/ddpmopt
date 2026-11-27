#! /bin/sh

while true; do
  yes 1*2^0 | p2-32 f 4 | head -n 1 | p2-32 P
  mv rand_pgm-0.pgm f.pgm
  yes 1*2^0 | p2-32 R2 1 | p2-32 f 4 | head -n 42 | p2-32 P
  ddpmoptmp p1 rand_pgm-[0-9].pgm rand_pgm-[0-9][0-9].pgm
  ddpmopt k f.pgm predg.ppm
done

