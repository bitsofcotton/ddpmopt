CXX=	clang++
#CXX=	/usr/local/bin/clang++
#CXX=	x86_64-unknown-openbsd7.1-eg++

# compiler flags.
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-fopenmp -L/usr/local/lib -lomp
#CXXFLAGS+=	-pg
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g2
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-Oz -mno-sse2 -mno-sse -mno-3dnow -mno-mmx -msoft-float -gfull
LDFLAGS+=	-lc++
#LDFLAGS+=	-lestdc++

clean:
	@rm -rf p0 p0-32 p0r p0r32
all:	p0 p0-32 p0r p0r32
p0:
	${CXX} ${CXXFLAGS} -static -O0 -o p0 p0.cc
p0-32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -O0 -o p0-32 p0.cc
p0-64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -O0 -o p0-64 p0.cc
p0r:
	${CXX} ${CXXFLAGS} -static -o p0r p0r.cc
p0r32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o p0r32 p0r.cc
p0r64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o p0r64 p0r.cc

