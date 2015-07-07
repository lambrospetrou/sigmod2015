#Profiling CFLAGS
# march=core2 
# march=corei7-avx
MATH_FLAGS=-ffast-math -funsafe-math-optimizations -fassociative-math -ffinite-math-only -fno-signed-zeros -funsafe-loop-optimizations -ftree-loop-if-convert-stores
RELEASE_CFLAGS=-march=native -std=c++11 -Ofast -O3 -W -Wall -Wextra -Wunused $(MATH_FLAGS) -fno-builtin -ftree-vectorize -funroll-all-loops -fvariable-expansion-in-unroller -fomit-frame-pointer -freorder-blocks-and-partition -Iinclude 

all: reference

# just compile always anyway ;)
reference: reference.cpp; g++ -g -L./include/asm  $(RELEASE_CFLAGS) -o$@ reference.cpp -lpthread -laelf64o -fabi-version=0 -fopt-info-vec #-fopt-info-vec-missed 

clean:
	rm reference
