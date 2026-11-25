cp rand_pgm-64.pgm last_rand_pgm-64.pgm
yes 1*2^0 | p2-32 R0 1 | p2-32 f 4 | head -n 65| p2-32 P
python3 ../goki predg.ppm repeat 21 ddpmoptmp p1 rand_pgm-[0-9].pgm rand_pgm-[1-5][0-9].pgm rand_pgm-6[0-3].pgm
convert rand_pgm-[0-9].pgm rand_pgm-[1-5][0-9].pgm rand_pgm-6[0-3].pgm -average -compress none avg.pgm
convert med.ppm -compress none med.pgm
ddpmopt d med.pgm rand_pgm-64.pgm
ddpmopt k rand_pgm-64.pgm med.pgm avg.pgm rand_pgm-63.pgm last_rand_pgm-64.pgm
ddpmopt k last_rand_pgm-64.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/0.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/1.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/2.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/3.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/4.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/5.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/6.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/7.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/8.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/9.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/a.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/b.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/c.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/d.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/e.pgm med.pgm avg.pgm rand_pgm-63.pgm
ddpmopt k ../pgm16/f.pgm med.pgm avg.pgm rand_pgm-63.pgm
p2-32 Y- rand_pgm-64.pgm | p2-32 V | p2-32 T+ | p2-32 F | tail -n 1

