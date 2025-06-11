CXX=	clang++
#CXX=	eg++
#CXX=	c++

# compiler flags.
##CXXFLAGS+=	-O0 -mtune=generic -gfull
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g3
# This doesn't work, we need operator >>, operator << with ongoing stdlibc++.
#CXXFLAGS+=	-I/usr/local/include -mlong-double-128
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-O2 -g3
#CXXFLAGS+=	-pg
#CXXFLAGS+=	--analyze
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lgomp -fopenmp
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-std=gnu++98
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib
# Same as -mlong-double-128
#LDFLAGS+=	-lquadmath -lm

# lieonn.hh compile options
CXXFLAGS+=	-D_ARCFOUR_

# N.B. sed -e s/static\ inline//g | sed -e s/inline//g
#CXXFLAGS+=     -D_OLDCPP_ -ftemplate-depth-99

CLEANFILES= *.o ddpmopt ddpmopt32 ddpmoptp ddpmoptmp ddpmopt32mp ddpmoptpmp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 ddpmoptp ddpmoptmp ddpmopt32mp ddpmoptpmp

ddpmopt:
	${CXX} ${CXXFLAGS} -static -o ddpmopt ddpmopt.cc
ddpmopt32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmopt32 ddpmopt.cc
ddpmopt64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmopt64 ddpmopt.cc
ddpmoptp:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -D_PERSISTENT_ -o ddpmoptp ddpmopt.cc
ddpmoptmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o ddpmoptmp ddpmopt.cc
ddpmopt32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o ddpmopt32mp ddpmopt.cc
ddpmopt64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o ddpmopt64mp ddpmopt.cc
ddpmoptpmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -D_PERSISTENT_ -o ddpmoptpmp ddpmopt.cc

