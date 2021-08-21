#pragma once

#include "../utils/objects.h"

void normalsDisplayMain(int argc, char *argv[]);
void rayTrace(Ray *ray, color *col, Object *Os);
void findFirstHit(Ray *ray, double *lambda, Object *Os, Object **obj, point *p, point *n, double *a, double *b);

