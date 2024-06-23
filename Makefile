CXX=	clang++
#CXX=	eg++

# compiler flags.
#CXXFLAGS+=	-O0 -mtune=generic -gfull
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g3
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-pg
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lgomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o ddpmopt ddpmopt32 predg predgp predg32 predgp32 qredg qredgp qredg32 qredgp32 ddpmoptmp ddpmopt32mp predgmp predgpmp predg32mp predgp32mp qredgmp qredgpmp qredg32mp qredgp32mp tcont tcont32

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 predg predg32 qredg qredg32 ddpmoptmp ddpmopt32mp predgmp predg32mp qredgmp qredg32mp tcont tcont32

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
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=false -o predg predg.cc
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=true  -o predgp predg.cc
predg32:
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=false -D_FLOAT_BITS_=32 -o predg32 predg.cc
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=true  -D_FLOAT_BITS_=32 -o predgp32 predg.cc
predgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=false -o predgmp predg.cc
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=true  -o predgpmp predg.cc
predg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=false -D_FLOAT_BITS_=32 -o predg32mp predg.cc
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=true  -D_FLOAT_BITS_=32 -o predgp32mp predg.cc
qredg:
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=false -o qredg qredg.cc
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=true  -o qredgp qredg.cc
qredg32:
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=false -D_FLOAT_BITS_=32 -o qredg32 qredg.cc
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_=true  -D_FLOAT_BITS_=32 -o qredgp32 qredg.cc
qredgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=false -o qredgmp qredg.cc
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=true  -o qredgpmp qredg.cc
qredg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=false -D_FLOAT_BITS_=32 -o qredg32mp qredg.cc
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_=true  -D_FLOAT_BITS_=32 -o qredgp32mp qredg.cc
tcont:
	${CXX} ${CXXFLAGS} -static -o tcont tcont.cc
tcont32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o tcont32 tcont.cc

