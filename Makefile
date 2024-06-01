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

CLEANFILES= *.o ddpmopt ddpmopt32 predg predg32 predg0 predg032 qredg qredg32 qredg0 qredg032 ddpmoptmp ddpmopt32mp predgmp predg32mp predg0mp predg032mp qredgmp qredg32mp qredg0mp qredg032mp tcont tcont32

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 predg predg32 predg0 predg032 qredg qredg32 qredg0 qredg032 ddpmoptmp ddpmopt32mp predgmp predg32mp predg0mp predg032mp qredgmp qredg32mp qredg0mp qredg032mp tcont tcont32

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
	${CXX} ${CXXFLAGS} -static -DP0J=false -o predg predg.cc
predg32:
	${CXX} ${CXXFLAGS} -static -DP0J=false -D_FLOAT_BITS_=32 -o predg32 predg.cc
predg64:
	${CXX} ${CXXFLAGS} -static -DP0J=false -D_FLOAT_BITS_=64 -o predg64 predg.cc
predg0:
	${CXX} ${CXXFLAGS} -static -DP0J=true -o predg0 predg.cc
predg032:
	${CXX} ${CXXFLAGS} -static -DP0J=true -D_FLOAT_BITS_=32 -o predg032 predg.cc
predg064:
	${CXX} ${CXXFLAGS} -static -DP0J=true -D_FLOAT_BITS_=64 -o predg064 predg.cc
predgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -o predgmp predg.cc
predg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -D_FLOAT_BITS_=32 -o predg32mp predg.cc
predg64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -D_FLOAT_BITS_=64 -o predg64mp predg.cc
predg128mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -D_FLOAT_BITS_=128 -o predg128mp predg.cc
predg0mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=true -o predg0mp predg.cc
predg032mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=true -D_FLOAT_BITS_=32 -o predg032mp predg.cc
predg064mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=true -D_FLOAT_BITS_=64 -o predg064mp predg.cc
qredg:
	${CXX} ${CXXFLAGS} -static -DP0J=false -o qredg qredg.cc
qredg32:
	${CXX} ${CXXFLAGS} -static -DP0J=false -D_FLOAT_BITS_=32 -o qredg32 qredg.cc
qredg64:
	${CXX} ${CXXFLAGS} -static -DP0J=false -D_FLOAT_BITS_=64 -o qredg64 qredg.cc
qredg0:
	${CXX} ${CXXFLAGS} -static -DP0J=true -o qredg0 qredg.cc
qredg032:
	${CXX} ${CXXFLAGS} -static -DP0J=true -D_FLOAT_BITS_=32 -o qredg032 qredg.cc
qredg064:
	${CXX} ${CXXFLAGS} -static -DP0J=true -D_FLOAT_BITS_=64 -o qredg064 qredg.cc
qredgmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -o qredgmp qredg.cc
qredg32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -D_FLOAT_BITS_=32 -o qredg32mp qredg.cc
qredg64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -D_FLOAT_BITS_=64 -o qredg64mp qredg.cc
qredg128mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=false -D_FLOAT_BITS_=128 -o qredg128mp qredg.cc
qredg0mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=true -o qredg0mp qredg.cc
qredg032mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=true -D_FLOAT_BITS_=32 -o qredg032mp qredg.cc
qredg064mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -DP0J=true -D_FLOAT_BITS_=64 -o qredg064mp qredg.cc
tcont:
	${CXX} ${CXXFLAGS} -static -o tcont tcont.cc
tcont32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o tcont32 tcont.cc

