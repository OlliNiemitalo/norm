#include <stdio.h>

int main() {
  // Set d to -nan
  double d = -0.0/0.0;
  printf("%d\n", d == d);
  return 0;
}

// g++ test.cpp -ffast-math -march=native -O3 -Wno-unused-result