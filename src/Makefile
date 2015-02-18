all: reference

reference: reference.cpp lputils.cpp
	g++ -std=c++11 -g -O3 -march=native -Wall -Wextra -Werror -o$@ reference.cpp lputils.cpp

clean:
	rm reference
