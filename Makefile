#CXX=	clang++
CXX=	eg++

# compiler flags.
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
#LDFLAGS+=	-lc++ -L/usr/local/lib
LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o ddpmopt ddpmopt32 ddpmopt64 enlarge enlarge32 enlarge64

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 ddpmopt64 enlarge enlarge32 enlarge64
ddpmopt:
	${CXX} ${CXXFLAGS} -static -o ddpmopt ddpmopt.cc
ddpmopt32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmopt32 ddpmopt.cc
ddpmopt64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmopt64 ddpmopt.cc
enlarge:
	${CXX} ${CXXFLAGS} -static -o enlarge enlarge.cc
enlarge32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o enlarge32 enlarge.cc
enlarge64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o enlarge64 enlarge.cc

