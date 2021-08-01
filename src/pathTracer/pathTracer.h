#include "../utils/objects.h"

#include <list>

#ifndef PATHTRACER_H
#define PATHTRACER_H
void pathTraceMain(int argc, char *argv[]);
void PathTrace(Ray *ray, int depth, color *col, Object *Os, Object *explicit_l);
void normalizeLightWeights(std::list<Object *>& object_list);
#endif