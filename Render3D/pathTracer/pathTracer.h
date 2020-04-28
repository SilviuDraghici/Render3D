#ifndef PATHTRACER_H
#define PATHTRACER_H
void pathTraceMain(int argc, char *argv[]);
void PathTrace(struct ray *ray, int depth, struct color *col, struct object *Os, struct object *explicit_l);
void normalizeLightWeights(struct object *object_list);
#endif