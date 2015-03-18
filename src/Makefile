#Profiling CFLAGS
# march=core2 - this was before by awfy
# march=corei7-avx
MATH_FLAGS=-ffast-math -funsafe-math-optimizations -fassociative-math -ffinite-math-only -fno-signed-zeros -funsafe-loop-optimizations -ftree-loop-if-convert-stores
RELEASE_CFLAGS=-march=native -std=c++11 -O3 -W -Wall -Wextra $(MATH_FLAGS) -fno-builtin -ftree-vectorize -funroll-all-loops -fvariable-expansion-in-unroller -fomit-frame-pointer -freorder-blocks-and-partition  -Iinclude 

#RELEASE_MIN=-march=corei7-avx -msse4.2 -std=c++11 -O3 -W -Wall -Wextra -funroll-all-loops -fomit-frame-pointer -Iinclude 

all: reference
    reference: reference.cpp; g++ -L./include/asm  $(RELEASE_CFLAGS) -fopenmp -o$@ reference.cpp -lpthread -laelf64o -fabi-version=0 -fopt-info-vec #-fopt-info-vec-missed 
#reference: reference.cpp; g++ $(RELEASE_CFLAGS) -I/usr/local/lib/tbb/include  -L/usr/local/lib/tbb/lib/intel64/gcc4.4 -fopenmp -o$@ reference.cpp -lpthread -ltbb -g

clean:
	rm reference
