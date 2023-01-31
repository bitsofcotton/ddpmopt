CXX=	clang++
#CXX=	eg++

# compiler flags.
#CXXFLAGS+=	-Ofast -mtune=native -gfull
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o ddpmopt ddpmopt32 ddpmoptp ddpmoptp32 ddpmoptq ddpmoptq32 ddpmoptmp ddpmopt32mp ddpmoptpmp ddpmoptp32mp ddpmoptqmp ddpmoptq32mp ddpmoptpm ddpmoptpm32 ddpmoptpmmp ddpmoptpm32mp ddpmoptqm ddpmoptqm32 ddpmoptqmmp ddpmoptqm32mp

clean:
	@rm -rf ${CLEANFILES}

all:	ddpmopt ddpmopt32 ddpmoptp ddpmoptp32 ddpmoptq ddpmoptq32 ddpmoptmp ddpmopt32mp ddpmoptpmp ddpmoptp32mp ddpmoptqmp ddpmoptq32mp ddpmoptpm ddpmoptpm32 ddpmoptpmmp ddpmoptpm32mp ddpmoptqm ddpmoptqm32 ddpmoptqmmp ddpmoptqm32mp
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
ddpmoptpm:
	${CXX} ${CXXFLAGS} -static -D_CONDORCET_JURY_ -o ddpmoptpm ddpmoptp.cc
ddpmoptpm32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -D_CONDORCET_JURY_ -o ddpmoptpm32 ddpmoptp.cc
ddpmoptpm64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -D_CONDORCET_JURY_ -o ddpmoptpm64 ddpmoptp.cc
ddpmoptq:
	${CXX} ${CXXFLAGS} -static -o ddpmoptq ddpmoptq.cc
ddpmoptq32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o ddpmoptq32 ddpmoptq.cc
ddpmoptq64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o ddpmoptq64 ddpmoptq.cc
ddpmoptqm:
	${CXX} ${CXXFLAGS} -static -D_CONDORCET_JURY_ -o ddpmoptqm ddpmoptq.cc
ddpmoptqm32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -D_CONDORCET_JURY_ -o ddpmoptqm32 ddpmoptq.cc
ddpmoptqm64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -D_CONDORCET_JURY_ -o ddpmoptqm64 ddpmoptq.cc
ddpmoptmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o ddpmoptmp ddpmopt.cc
ddpmopt32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o ddpmopt32mp ddpmopt.cc
ddpmopt64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o ddpmopt64mp ddpmopt.cc
ddpmoptpmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o ddpmoptpmp ddpmoptp.cc
ddpmoptp32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o ddpmoptp32mp ddpmoptp.cc
ddpmoptp64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o ddpmoptp64mp ddpmoptp.cc
ddpmoptpmmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONDORCET_JURY_ -o ddpmoptpmmp ddpmoptp.cc
ddpmoptpm32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -D_CONDORCET_JURY_ -o ddpmoptpm32mp ddpmoptp.cc
ddpmoptpm64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -D_CONDORCET_JURY_ -o ddpmoptpm64mp ddpmoptp.cc
ddpmoptqmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -o ddpmoptqmp ddpmoptq.cc
ddpmoptq32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -o ddpmoptq32mp ddpmoptq.cc
ddpmoptq64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -o ddpmoptq64mp ddpmoptq.cc
ddpmoptqmmp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_CONDORCET_JURY_ -o ddpmoptqmmp ddpmoptq.cc
ddpmoptqm32mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=32 -D_CONDORCET_JURY_ -o ddpmoptqm32mp ddpmoptq.cc
ddpmoptqm64mp:
	${CXX} ${CXXFLAGS} ${MPFLAGS} -D_FLOAT_BITS_=64 -D_CONDORCET_JURY_ -o ddpmoptqm64mp ddpmoptq.cc

