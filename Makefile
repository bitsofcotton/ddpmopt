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

CLEANFILES= *.o ddpmopt ddpmopt32 predg predg32 predg3 predg3-32 qredg qredg32 qredg3 qredg3-32 ddpmoptmp ddpmopt32mp predgmp predg32mp predg3mp predg3-32mp qredgmp qredg32mp qredg3mp qredg3-32mp tcont tcont32

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 predg predg32 predg3 predg3-32 qredg qredg32 qredg3 qredg3-32 ddpmoptmp ddpmopt32mp predgmp predg32mp predg3mp predg3-32mp qredgmp qredg32mp qredg3mp qredg3-32mp tcont tcont32

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
predg3:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=3 -o predg3 predg.cc
predg3-32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=3 -D_FLOAT_BITS_=32 -o predg3-32 predg.cc
predg6:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=6 -o predg6 predg.cc
predg6-32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=6 -D_FLOAT_BITS_=32 -o predg6-32 predg.cc
predg9:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=9 -o predg9 predg.cc
predg9-32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=9 -D_FLOAT_BITS_=32 -o predg9-32 predg.cc
predgp:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=-1 -o predgp predg.cc
predgp32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=-1 -D_FLOAT_BITS_=32 -o predgp32 predg.cc
predgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o predgmp predg.cc
predg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o predg32mp predg.cc
predg3mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=3 -o predg3mp predg.cc
predg3-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=3 -D_FLOAT_BITS_=32 -o predg3-32mp predg.cc
predg6mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=6 -o predg6mp predg.cc
predg6-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=6 -D_FLOAT_BITS_=32 -o predg6-32mp predg.cc
predg9mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=9 -o predg9mp predg.cc
predg9-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=9 -D_FLOAT_BITS_=32 -o predg9-32mp predg.cc
predgpmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=-1 -o predgpmp predg.cc
predgp32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=-1 -D_FLOAT_BITS_=32 -o predgp32mp predg.cc
qredg:
	${CXX} ${CXXFLAGS} -static -o qredg qredg.cc
qredg32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o qredg32 qredg.cc
qredg3:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=3 -o qredg3 qredg.cc
qredg3-32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=3 -D_FLOAT_BITS_=32 -o qredg3-32 qredg.cc
qredg6:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=6 -o qredg6 qredg.cc
qredg6-32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=6 -D_FLOAT_BITS_=32 -o qredg6-32 qredg.cc
qredg9:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=9 -o qredg9 qredg.cc
qredg9-32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=9 -D_FLOAT_BITS_=32 -o qredg9-32 qredg.cc
qredgp:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=-1 -o qredgp qredg.cc
qredgp32:
	${CXX} ${CXXFLAGS} -static -D_PREDV_=-1 -D_FLOAT_BITS_=32 -o qredgp32 qredg.cc
qredgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o qredgmp qredg.cc
qredg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o qredg32mp qredg.cc
qredg3mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=3 -o qredg3mp qredg.cc
qredg3-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=3 -D_FLOAT_BITS_=32 -o qredg3-32mp qredg.cc
qredg6mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=6 -o qredg6mp qredg.cc
qredg6-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=6 -D_FLOAT_BITS_=32 -o qredg6-32mp qredg.cc
qredg9mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=9 -o qredg9mp qredg.cc
qredg9-32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=9 -D_FLOAT_BITS_=32 -o qredg9-32mp qredg.cc
qredgpmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=-1 -o qredgpmp qredg.cc
qredgp32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PREDV_=-1 -D_FLOAT_BITS_=32 -o qredgp32mp qredg.cc
tcont:
	${CXX} ${CXXFLAGS} -static -o tcont tcont.cc
tcont32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o tcont32 tcont.cc

