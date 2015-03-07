#Profiling CFLAGS
# march=core2 - this was before by awfy
MATH_FLAGS=-ffast-math -funsafe-math-optimizations -fassociative-math -ffinite-math-only -fno-signed-zeros
CFLAGS=-march=native -msse4.1 -std=c++11 -O3 -W -Wall -Wextra -DDEBUG -DEXPBACKOFF -g -pedantic $(MATH_FLAGS) -funroll-all-loops -fvariable-expansion-in-unroller -fno-omit-frame-pointer -Iinclude
RELEASE_CFLAGS=-march=corei7-avx -msse4.2 -std=c++11 -O3 -W -Wall -Wextra -DNDBGPRINT $(MATH_FLAGS) -ftree-vectorize -funroll-all-loops -fvariable-expansion-in-unroller -fomit-frame-pointer -Iinclude

all: reference

reference: reference.cpp; g++ $(RELEASE_CFLAGS) -I/usr/local/lib/tbb/include  -L/usr/local/lib/tbb/lib/intel64/gcc4.4 -fopenmp -o$@ reference.cpp -lpthread -ltbb
#	g++ -fcilkplus -fopenmp $(RELEASE_CFLAGS) -o$@ reference.cpp -lpthread 

clean:
	rm reference
