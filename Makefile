CXX=	clang++
#CXX=	eg++
#CXX=	c++

# compiler flags.
#CXXFLAGS+=	-O0 -mtune=generic -gfull
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g3
# This doesn't work, we need operator >>, operator << with ongoing stdlibc++.
#CXXFLAGS+=	-I/usr/local/include -mlong-double-128
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O2 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-O2 -g3
#CXXFLAGS+=	-mno-sse2 -mno-sse -mno-3dnow -mno-mmx -msoft-float -fno-omit-frame-pointer
#CXXFLAGS+=	-pg
#CXXFLAGS+=	--analyze
#CXXFLAGS+=      -D_LIBCPP_HARDENING_MODE_DEBUG
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lgomp -fopenmp
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-std=gnu++98
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib
# Same as -mlong-double-128
#LDFLAGS+=	-lquadmath -lm

# lieonn.hh compile options
#CXXFLAGS+=	-D_P_PRNG_=11
# N.B. sed -e s/static\ inline//g | sed -e s/inline//g
#CXXFLAGS+=     -D_OLDCPP_ -ftemplate-depth-99
# N.B. _SIMPLEALLOC_=align needs env VM_LIEONN=(mem usage MB).

CLEANFILES= *.o ddpmopta ddpmoptpa ddpmoptmp ddpmoptpmp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopta ddpmoptpa ddpmoptmp ddpmoptpmp

ddpmopt:
	${CXX} ${CXXFLAGS} -static -o ddpmopt ddpmopt.cc
ddpmopt32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmopt32 ddpmopt.cc
ddpmopt64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmopt64 ddpmopt.cc
ddpmoptp:
	${CXX} ${CXXFLAGS} -static -D_PERSISTENT_ -o ddpmoptp ddpmopt.cc
ddpmopta:
	${CXX} ${CXXFLAGS} -static -D_SIMPLEALLOC_=64 -o ddpmopta ddpmopt.cc
ddpmopt32a:
	${CXX} ${CXXFLAGS} -static -D_SIMPLEALLOC_=64 -D_FLOAT_BITS_=32 -o ddpmopt32a ddpmopt.cc
ddpmoptpa:
	${CXX} ${CXXFLAGS} -static -D_SIMPLEALLOC_=64 -D_PERSISTENT_ -o ddpmoptpa ddpmopt.cc
ddpmoptmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o ddpmoptmp ddpmopt.cc
ddpmopt32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o ddpmopt32mp ddpmopt.cc
ddpmopt64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o ddpmopt64mp ddpmopt.cc
ddpmoptpmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_PERSISTENT_ -o ddpmoptpmp ddpmopt.cc

