#ifndef MAPPING_H
#define MAPPING_H


#include "color.h"
#include "objects.h"
#include "scene.h"

void loadTexture(Object *o, const std::string& filename, int type, Scene *scene);

void textureMap(Object *obj, double a, double b, color *col);
void normalMap(Object *obj, double a, double b, struct point *n);

void alphaMap(Object *obj, double a, double b, double *set_alpha, double obj_alpha);
#endif