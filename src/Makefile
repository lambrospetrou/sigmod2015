#Profiling CFLAGS
# march=core2 - this was before by awfy
MATH_FLAGS=-ffast-math -funsafe-math-optimizations -fassociative-math -ffinite-math-only -fno-signed-zeros
CFLAGS=-march=native -msse4.1 -std=c++11 -O3 -W -Wall -Wextra -DDEBUG -DEXPBACKOFF -g -pedantic $(MATH_FLAGS) -funroll-all-loops -fvariable-expansion-in-unroller -fno-omit-frame-pointer -Iinclude
RELEASE_CFLAGS=-march=native -msse4.1 -std=c++11 -O3 -W -Wall -Wextra -Werror -DNDEBUG -DNDBGPRINT -DEXPBACKOFF $(MATH_FLAGS) -funroll-all-loops -fvariable-expansion-in-unroller -Iinclude

all: reference

reference: reference.cpp 
	g++ $(RELEASE_CFLAGS) -Llib/tbb/lib/intel64/gcc4.4 -ltbb  -o$@ reference.cpp -lpthread 

clean:
	rm reference
