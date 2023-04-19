CXX=	clang++
#CXX=	eg++

# compiler flags.
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o ddpmopt ddpmopt32 ddpmoptmp ddpmopt32mp predg predg32 predgmp predg32mp qredg qredg32 qredgmp qredg32mp topt topt32 toptmp topt32mp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 ddpmoptmp ddpmopt32mp predg predg32 predgmp predg32mp qredg qredg32 qredgmp qredg32mp topt topt32 toptmp topt32mp

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
topt:
	${CXX} ${CXXFLAGS} -static -o topt topt.cc
topt32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o topt32 topt.cc
topt64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o topt64 topt.cc
toptmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o toptmp topt.cc
topt32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o topt32mp topt.cc
topt64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o topt64mp topt.cc
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
predg128mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=128 -o predg128mp predg.cc
predg256mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=256 -o predg256mp predg.cc
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
qredg128mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=128 -o qredg128mp qredg.cc
qredg256mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=256 -o qredg256mp qredg.cc

