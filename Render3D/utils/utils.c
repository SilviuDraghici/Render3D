#include "utils.h"

#include "ray.h"

#include <math.h>

double dot(struct point *u, struct point *v) {
   // Computes the dot product of 3D vectors u and v.
   // The function assumes the w components of both vectors
   // are 1.
   return ((u->x * v->x) + (u->y * v->y) + (u->z * v->z));
}

struct point *cross(struct point *u, struct point *v)
{
   // Allocates and returns a vector with the cross product u x v.
   // The function assumes the w components of both vectors
   // are 1.
   struct point *cp;
   cp = (struct point *)calloc(1, sizeof(struct point));

   cp->x = (u->y * v->z) - (v->y * u->z);
   cp->y = (v->x * u->z) - (u->x * v->z);
   cp->z = (u->x * v->y) - (v->x * u->y);
   cp->w = 1;
   return (cp);
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

void solveQuadratic(struct ray *ray, double *l1, double *l2) {
   struct point a_sub_c;
   double A, B, C, D;
   a_sub_c = ray->p0;
   A = dot(&(ray->d), &(ray->d));
   B = dot(&a_sub_c, &(ray->d));
   C = dot(&a_sub_c, &a_sub_c) - 1;
   D = B * B - A * C;

   if (D < 0) {
      return;
   }

   if (D > 0) {
      *l1 = -B / A - sqrt(D) / A;
      *l2 = -B / A + sqrt(D) / A;
   }
}