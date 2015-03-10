#include <cstdio>

#define OP_AB(a, b, op) a##op##b;

int main() {

    printf("hello: %d\n", OP_AB(1,2,+));
    return 0;
}
