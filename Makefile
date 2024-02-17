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

CLEANFILES= *.o ddpmopt ddpmopt32 predg predg32 qredg qredg32 ddpmoptmp ddpmopt32mp predgmp predg32mp qredgmp qredg32mp predcg predcg32 qredcg qredcg32 predcgmp predcg32mp qredcgmp qredcg32mp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 predg predg32 qredg qredg32 ddpmoptmp ddpmopt32mp predgmp predg32mp qredgmp qredg32mp predcg predcg32 qredcg qredcg32 predcgmp predcg32mp qredcgmp qredcg32mp

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
predcg:
	${CXX} ${CXXFLAGS} -D_CONTINUOUS_ -static -o predcg predg.cc
predcg32:
	${CXX} ${CXXFLAGS} -D_CONTINUOUS_ -static -D_FLOAT_BITS_=32 -o predcg32 predg.cc
predcg64:
	${CXX} ${CXXFLAGS} -D_CONTINUOUS_ -static -D_FLOAT_BITS_=64 -o predcg64 predg.cc
predcgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONTINUOUS_ -o predcgmp predg.cc
predcg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONTINUOUS_ -D_FLOAT_BITS_=32 -o predcg32mp predg.cc
predcg64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONTINUOUS_ -D_FLOAT_BITS_=64 -o predcg64mp predg.cc
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
qredcg:
	${CXX} ${CXXFLAGS} -D_CONTINUOUS_ -static -o qredcg qredg.cc
qredcg32:
	${CXX} ${CXXFLAGS} -D_CONTINUOUS_ -static -D_FLOAT_BITS_=32 -o qredcg32 qredg.cc
qredcg64:
	${CXX} ${CXXFLAGS} -D_CONTINUOUS_ -static -D_FLOAT_BITS_=64 -o qredcg64 qredg.cc
qredcgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONTINUOUS_ -o qredcgmp qredg.cc
qredcg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONTINUOUS_ -D_FLOAT_BITS_=32 -o qredcg32mp qredg.cc
qredcg64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONTINUOUS_ -D_FLOAT_BITS_=64 -o qredcg64mp qredg.cc

