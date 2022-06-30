CXX=	clang++
#CXX=	/usr/local/bin/eg++

# compiler flags.
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-Oz -mtune=native -gfull
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
CXXFLAGS+=	-std=c++11
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib

CLEANFILES= *.o tools

clean:
	@rm -rf ${CLEANFILES}

all:	tools

