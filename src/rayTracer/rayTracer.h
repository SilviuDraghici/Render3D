#ifndef RAYTRACER_H
#define RAYTRACER_H
void rayTraceMain(int argc, char *argv[]);
void rayTrace(struct ray *ray, int depth, struct color *col, struct object *Os);
void findFirstHit(struct ray *ray, double *lambda, struct object *Os, struct object **obj, struct point *p, struct point *n, double *a, double *b);
void rtShade(struct object *obj, struct point *p, struct point *n, struct ray *ray, int depth, double a, double b, struct color *col);
void rt_brandished_trace(struct ray *ray, struct object *obj, struct color *col, int depth);
#endif