#include "../utils/objects.h"

#ifndef PATHTRACER_H
#define PATHTRACER_H
void pathTraceMain(int argc, char *argv[]);
void PathTrace(struct ray *ray, int depth, struct color *col, Object *Os, Object *explicit_l);
void normalizeLightWeights(Object *object_list);
#endif