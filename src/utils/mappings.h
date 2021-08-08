#ifndef MAPPING_H
#define MAPPING_H

#include <list>

#include "color.h"
#include "objects.h"
#include "scene.h"


textureNode* loadTexture(const std::string& filename, int type, std::list<textureNode*>& texture_list);
void setTexture(Object *o, const std::string& filename, int type, std::list<textureNode*>& texture_list);

void textureMap(Object *obj, double a, double b, color *col);
void normalMap(Object *obj, double a, double b, struct point *n);

void alphaMap(Object *obj, double a, double b, double *set_alpha, double obj_alpha);
#endif