CXX=	clang++
#CXX=	eg++

# compiler flags.
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

CLEANFILES= *.o ddpmopt ddpmopt32 predg predg32 predgn predgn32 qredg qredg32 qredgn qredgn32 ddpmoptmp ddpmopt32mp predgmp predg32mp predgnmp predgn32mp qredgmp qredg32mp qredgnmp qredgn32mp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 predg predg32 predgn predgn32 qredg qredg32 qredgn qredgn32 ddpmoptmp ddpmopt32mp predgmp predg32mp predgnmp predgn32mp qredgmp qredg32mp qredgnmp qredgn32mp

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
	${CXX} ${CXXFLAGS} -static -o predg predg.cc
predg32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o predg32 predg.cc
predg64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o predg64 predg.cc
predgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o predgmp predg.cc
predg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o predg32mp predg.cc
predg64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o predg64mp predg.cc
predgn:
	${CXX} ${CXXFLAGS} -static -DNOCOMP -o predgn predg.cc
predgn32:
	${CXX} ${CXXFLAGS} -static -DNOCOMP -D_FLOAT_BITS_=32 -o predgn32 predg.cc
predgn64:
	${CXX} ${CXXFLAGS} -static -DNOCOMP -D_FLOAT_BITS_=64 -o predgn64 predg.cc
predgnmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DNOCOMP -o predgnmp predg.cc
predgn32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DNOCOMP -D_FLOAT_BITS_=32 -o predgn32mp predg.cc
predgn64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DNOCOMP -D_FLOAT_BITS_=64 -o predgn64mp predg.cc
qredg:
	${CXX} ${CXXFLAGS} -static -o qredg qredg.cc
qredg32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o qredg32 qredg.cc
qredg64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o qredg64 qredg.cc
qredgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o qredgmp qredg.cc
qredg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o qredg32mp qredg.cc
qredg64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o qredg64mp qredg.cc
qredgn:
	${CXX} ${CXXFLAGS} -static -DNOCOMP -o qredgn qredg.cc
qredgn32:
	${CXX} ${CXXFLAGS} -static -DNOCOMP -D_FLOAT_BITS_=32 -o qredgn32 qredg.cc
qredgn64:
	${CXX} ${CXXFLAGS} -static -DNOCOMP -D_FLOAT_BITS_=64 -o qredgn64 qredg.cc
qredgnmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DNOCOMP -o qredgnmp qredg.cc
qredgn32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DNOCOMP -D_FLOAT_BITS_=32 -o qredgn32mp qredg.cc
qredgn64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DNOCOMP -D_FLOAT_BITS_=64 -o qredgn64mp qredg.cc

