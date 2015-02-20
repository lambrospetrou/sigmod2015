all: reference

reference: reference.cpp BoundedSRSWQueue.hpp
	g++ -std=c++11 -O3 -march=native -Wall -Wextra -Werror -o$@ BoundedSRSWQueue.hpp LockFreeBoundedSRSWQueue.hpp reference.cpp -lpthread

clean:
	rm reference
