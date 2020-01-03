#include "utils.h"

#include <math.h>

double dot(struct point *u, struct point *v) {
   // Computes the dot product of 3D vectors u and v.
   // The function assumes the w components of both vectors
   // are 1.
   return ((u->x * v->x) + (u->y * v->y) + (u->z * v->z));
}

void normalize(struct point *v)
{
   // Normalizes a vector to unit length.
   // Note that this assumes the w component of v is one, so make
   // sure you don't have homogeneous vectors with w != 1
   // floating around or things will go south.
   double l;
   l = v->x * v->x;
   l += (v->y * v->y);
   l += (v->z * v->z);
   l = 1.0 / sqrt(l);
   v->x *= l;
   v->y *= l;
   v->z *= l;
}