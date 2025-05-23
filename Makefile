CXX=	clang++
#CXX=	eg++

# compiler flags.
##CXXFLAGS+=	-O0 -mtune=generic -gfull
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g3
# This doesn't work, we need operator >>, operator << with ongoing stdlibc++.
#CXXFLAGS+=	-I/usr/local/include -mlong-double-128
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-pg
#CXXFLAGS+=	--analyze
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lgomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib
# Same as -mlong-double-128
#LDFLAGS+=	-lquadmath -lm

# lieonn.hh compile options
#CXXFLAGS+=	-D_ARCFOUR_
CXXFLAGS+=	-D_FEED_MUCH_

CLEANFILES= *.o ddpmopt ddpmopt32 ddpmopt64 ddpmoptmp ddpmopt32mp ddpmopt64mp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 ddpmopt64 ddpmoptmp ddpmopt32mp ddpmopt64mp

ddpmopt:
	${CXX} ${CXXFLAGS} -static -D_ARCFOUR_ -o ddpmopt ddpmopt.cc
ddpmopt32:
	${CXX} ${CXXFLAGS} -static -D_ARCFOUR_ -D_FLOAT_BITS_=32 -o ddpmopt32 ddpmopt.cc
ddpmopt64:
	${CXX} ${CXXFLAGS} -static -D_ARCFOUR_ -D_FLOAT_BITS_=64 -o ddpmopt64 ddpmopt.cc
ddpmoptmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_ARCFOUR_ -o ddpmoptmp ddpmopt.cc
ddpmopt32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_ARCFOUR_ -D_FLOAT_BITS_=32 -o ddpmopt32mp ddpmopt.cc
ddpmopt64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_ARCFOUR_ -D_FLOAT_BITS_=64 -o ddpmopt64mp ddpmopt.cc
ddpmopt128mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_ARCFOUR_ -D_FLOAT_BITS_=128 -o ddpmopt128mp ddpmopt.cc

