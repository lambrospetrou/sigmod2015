all: reference

reference: reference.cpp lputils.cpp
	g++ -std=c++11 -O3 -march=native -Wall -Wextra -Werror -o$@ reference.cpp lputils.cpp -lpthread

clean:
	rm reference
