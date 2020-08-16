#include "../utils/objects.h"

#ifndef RAYTRACER_H
#define RAYTRACER_H
void rayTraceMain(int argc, char *argv[]);
void rayTrace(Ray *ray, int depth, color *col, Object *Os);
void findFirstHit(Ray *ray, double *lambda, Object *Os, Object **obj, point *p, point *n, double *a, double *b);
void rtShade(Object *obj, point *p, point *n, Ray *ray, int depth, double a, double b, color *col);
void rt_brandished_trace(Ray *ray, Object *obj, color *col, int depth);
#endif