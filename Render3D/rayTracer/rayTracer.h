#ifndef RAYTRACER_H
#define RAYTRACER_H
void rayTraceMain(int argc, char *argv[]);
void rayTrace(struct ray *ray, int depth, struct color *col, struct object3D *Os);
void findFirstHit(struct ray *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point *p, struct point *n, double *a, double *b);
void rtShade(struct object3D *obj, struct point *p, struct point *n, struct ray *ray, int depth, double a, double b, struct color *col);
#endif