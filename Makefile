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

CLEANFILES= *.o ddpmopt ddpmopt32 predg predg32 predga predga32 qredg qredg32 ddpmoptmp ddpmopt32mp predgmp predg32mp predgamp predga32mp qredgmp qredg32mp tcont tcont32

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 predg predg32 predga predga32 qredg qredg32 ddpmoptmp ddpmopt32mp predgmp predg32mp predgamp predga32mp qredgmp qredg32mp tcont tcont32

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
predgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o predgmp predg.cc
predg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o predg32mp predg.cc
predg64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o predg64mp predg.cc
predga:
	${CXX} ${CXXFLAGS} -static -o predga predga.cc
predga32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o predga32 predga.cc
predgamp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o predgamp predga.cc
predga32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o predga32mp predga.cc
predga64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o predga64mp predga.cc
qredg:
	${CXX} ${CXXFLAGS} -static -o qredg qredg.cc
qredg32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o qredg32 qredg.cc
qredgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o qredgmp qredg.cc
qredg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o qredg32mp qredg.cc
tcont:
	${CXX} ${CXXFLAGS} -static -o tcont tcont.cc
tcont32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o tcont32 tcont.cc

