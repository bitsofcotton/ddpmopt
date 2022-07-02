CXX=	clang++
#CXX=	/usr/local/bin/eg++

# compiler flags.
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-Oz -mtune=native -gfull
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o ddpmopt ddpmoptO0 ddpmopt64

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt
ddpmopt:
	${CXX} ${CXXFLAGS} -static -o ddpmopt tools.cc
ddpmoptO0:
	${CXX} ${CXXFLAGS} -static -O0 -o ddpmoptO0 tools.cc

