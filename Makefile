CXX=	clang++
#CXX=	eg++

# compiler flags.
#CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
CXXFLAGS+=	-O2 -mtune=native -gfull
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o ddpmopt ddpmopt32 ddpmopt64 ddpmoptp ddpmoptp32 ddpmoptp64 ddpmoptq ddpmoptq32 ddpmoptq64

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 ddpmopt64 ddpmoptp ddpmoptp32 ddpmoptp64 ddpmoptq ddpmoptq32 ddpmoptq64
ddpmopt:
	${CXX} ${CXXFLAGS} -static -o ddpmopt ddpmopt.cc
ddpmopt32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmopt32 ddpmopt.cc
ddpmopt64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmopt64 ddpmopt.cc
ddpmoptp:
	${CXX} ${CXXFLAGS} -static -o ddpmoptp ddpmoptp.cc
ddpmoptp32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmoptp32 ddpmoptp.cc
ddpmoptp64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmoptp64 ddpmoptp.cc
ddpmoptq:
	${CXX} ${CXXFLAGS} -static -o ddpmoptq ddpmoptq.cc
ddpmoptq32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmoptq32 ddpmoptq.cc
ddpmoptq64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmoptq64 ddpmoptq.cc

