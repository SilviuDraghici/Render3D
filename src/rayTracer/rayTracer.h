#include "../utils/objects.h"

#ifndef RAYTRACER_H
#define RAYTRACER_H
void rayTraceMain(int argc, char *argv[]);
void rayTrace(struct ray *ray, int depth, struct color *col, Object *Os);
void findFirstHit(struct ray *ray, double *lambda, Object *Os, Object **obj, struct point *p, struct point *n, double *a, double *b);
void rtShade(Object *obj, struct point *p, struct point *n, struct ray *ray, int depth, double a, double b, struct color *col);
void rt_brandished_trace(struct ray *ray, Object *obj, struct color *col, int depth);
#endif