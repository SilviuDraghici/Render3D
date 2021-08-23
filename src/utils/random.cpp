#include "random.h"
#include <random>


thread_local static int rng_x = 123456789;
thread_local static int rng_y = 362436069;
thread_local static int rng_z = 521288629;
thread_local static int rng_w = 88675123;

double xor128(void) {
  //https://en.wikipedia.org/wiki/Xorshift
  int t;
  t = rng_x ^ (rng_x << 11);
  rng_x = rng_y;
  rng_y = rng_z;
  rng_z = rng_w;
  return (rng_w = rng_w ^ (rng_w >> 19) ^ (t ^ (t >> 8))) / 2147483647.0;
}