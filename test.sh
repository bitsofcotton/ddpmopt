yes 1*2^0 | p2-32 R0 1 | p2-32 f 4 | head -n 65| p2-32 P
python3 ../goki predg.ppm repeat 5 ddpmoptmp p 1 rand_pgm-[0-9].pgm rand_pgm-[1-5][0-9].pgm rand_pgm-6[0-3].pgm
convert rand_pgm-[0-9].pgm rand_pgm-[1-5][0-9].pgm rand_pgm-6[0-3].pgm -average -compress none avg.pgm
ddpmopt m+ pred*m
ddpmopt d med.ppm rand_pgm-64.pgm
ddpmopt k rand_pgm-64.pgm med.ppm avg.pgm rand_pgm-63.pgm
ddpmopt k+ rand_pgm-64.pgm med.ppm avg.pgm rand_pgm-63.pgm

