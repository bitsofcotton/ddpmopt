CXX=	clang++
#CXX=	eg++

# compiler flags.
#CXXFLAGS+=	-O0 -mtune=generic -gfull
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g3
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-pg
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lgomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o ddpmopt ddpmopt32 predg predg32 qredg qredg32 ddpmoptmp ddpmopt32mp predgmp predg32mp qredgmp qredg32mp tcont tcont32 preddgmp qreddgmp preddg32mp qreddg32mp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 predg predg32 qredg qredg32 ddpmoptmp ddpmopt32mp predgmp predg32mp qredgmp qredg32mp tcont tcont32 preddgmp preddg32mp qreddgmp qreddg32mp

ddpmopt:
	${CXX} ${CXXFLAGS} -static -o ddpmopt ddpmopt.cc
ddpmopt32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmopt32 ddpmopt.cc
ddpmopt64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmopt64 ddpmopt.cc
ddpmoptmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o ddpmoptmp ddpmopt.cc
ddpmopt32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o ddpmopt32mp ddpmopt.cc
ddpmopt64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o ddpmopt64mp ddpmopt.cc
predg:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -o predg predg.cc
preddg:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_DFT_ -o preddg predg.cc
predg32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_FLOAT_BITS_=32 -o predg32 predg.cc
preddg32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_DFT_ -D_FLOAT_BITS_=32 -o preddg32 predg.cc
predgn:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -o predgn predg.cc
predgn32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_FLOAT_BITS_=32 -o predgn32 predg.cc
predg3:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=3 -o predg3 predg.cc
predg3-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=3 -D_FLOAT_BITS_=32 -o predg3-32 predg.cc
predgn3:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=3 -o predgn3 predg.cc
predgn3-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=3 -D_FLOAT_BITS_=32 -o predgn3-32 predg.cc
predg6:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=6 -o predg6 predg.cc
predg6-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=6 -D_FLOAT_BITS_=32 -o predg6-32 predg.cc
predgn6:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=6 -o predgn6 predg.cc
predgn6-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=6 -D_FLOAT_BITS_=32 -o predgn6-32 predg.cc
predg9:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=9 -o predg9 predg.cc
predg9-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=9 -D_FLOAT_BITS_=32 -o predg9-32 predg.cc
predgn9:
	${CXX} ${CXXFLAGS} -static -D_PNOISE=true -D_PREDV_=9 -o predgn9 predg.cc
predgn9-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE=true -D_PREDV_=9 -D_FLOAT_BITS_=32 -o predgn9-32 predg.cc
predgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -o predgmp predg.cc
preddgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_DFT_ -o preddgmp predg.cc
predg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_FLOAT_BITS_=32 -o predg32mp predg.cc
preddg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_DFT_ -D_FLOAT_BITS_=32 -o preddg32mp predg.cc
predgnmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -o predgnmp predg.cc
predgn32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_FLOAT_BITS_=32 -o predgn32mp predg.cc
predg3mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=3 -o predg3mp predg.cc
predg3-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=3 -D_FLOAT_BITS_=32 -o predg3-32mp predg.cc
predgn3mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=3 -o predgn3mp predg.cc
predgn3-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=3 -D_FLOAT_BITS_=32 -o predgn3-32mp predg.cc
predg6mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=6 -o predg6mp predg.cc
predg6-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=6 -D_FLOAT_BITS_=32 -o predg6-32mp predg.cc
predgn6mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=6 -o predgn6mp predg.cc
predgn6-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=6 -D_FLOAT_BITS_=32 -o predgn6-32mp predg.cc
predg9mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=9 -o predg9mp predg.cc
predg9-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=9 -D_FLOAT_BITS_=32 -o predg9-32mp predg.cc
predgn9mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=9 -o predgn9mp predg.cc
predgn9-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=9 -D_FLOAT_BITS_=32 -o predgn9-32mp predg.cc
qredg:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -o qredg qredg.cc
qreddg:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_DFT_ -o qreddg qredg.cc
qredg32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_FLOAT_BITS_=32 -o qredg32 qredg.cc
qreddg32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_DFT_ -D_FLOAT_BITS_=32 -o qreddg32 qredg.cc
qredgn:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -o qredgn qredg.cc
qredgn32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_FLOAT_BITS_=32 -o qredgn32 qredg.cc
qredg3:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=3 -o qredg3 qredg.cc
qredg3-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=3 -D_FLOAT_BITS_=32 -o qredg3-32 qredg.cc
qredgn3:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=3 -o qredgn3 qredg.cc
qredgn3-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=3 -D_FLOAT_BITS_=32 -o qredgn3-32 qredg.cc
qredg6:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=6 -o qredg6 qredg.cc
qredg6-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=6 -D_FLOAT_BITS_=32 -o qredg6-32 qredg.cc
qredgn6:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=6 -o qredgn6 qredg.cc
qredgn6-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=6 -D_FLOAT_BITS_=32 -o qredgn6-32 qredg.cc
qredg9:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=9 -o qredg9 qredg.cc
qredg9-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=false -D_PREDV_=9 -D_FLOAT_BITS_=32 -o qredg9-32 qredg.cc
qredgn9:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=9 -o qredgn9 qredg.cc
qredgn9-32:
	${CXX} ${CXXFLAGS} -static -D_PNOISE_=true -D_PREDV_=9 -D_FLOAT_BITS_=32 -o qredgn9-32 qredg.cc
qredgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -o qredgmp qredg.cc
qreddgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_DFT_ -o qreddgmp qredg.cc
qredg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_FLOAT_BITS_=32 -o qredg32mp qredg.cc
qreddg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_DFT_ -D_FLOAT_BITS_=32 -o qreddg32mp qredg.cc
qredgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -o qredgnmp qredg.cc
qredg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_FLOAT_BITS_=32 -o qredgn32mp qredg.cc
qredg3mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=3 -o qredg3mp qredg.cc
qredg3-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=3 -D_FLOAT_BITS_=32 -o qredg3-32mp qredg.cc
qredgn3mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=3 -o qredgn3mp qredg.cc
qredgn3-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=3 -D_FLOAT_BITS_=32 -o qredgn3-32mp qredg.cc
qredg6mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=6 -o qredg6mp qredg.cc
qredg6-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=6 -D_FLOAT_BITS_=32 -o qredg6-32mp qredg.cc
qredgn6mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=6 -o qredgn6mp qredg.cc
qredgn6-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=6 -D_FLOAT_BITS_=32 -o qredgn6-32mp qredg.cc
qredg9mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=9 -o qredg9mp qredg.cc
qredg9-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=false -D_PREDV_=9 -D_FLOAT_BITS_=32 -o qredg9-32mp qredg.cc
qredgn9mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=9 -o qredgn9mp qredg.cc
qredgn9-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PNOISE_=true -D_PREDV_=9 -D_FLOAT_BITS_=32 -o qredgn9-32mp qredg.cc
tcont:
	${CXX} ${CXXFLAGS} -static -o tcont tcont.cc
tcont32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o tcont32 tcont.cc

