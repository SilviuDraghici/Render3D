#ifndef MAPPING_H
#define MAPPING_H

#include <list>

#include "color.h"
#include "objects.h"
#include "scene.h"


textureNode* load_texture(const std::string& filename, int type, std::list<textureNode*>& texture_list);
void set_texture(Object *o, const std::string& filename, int type, std::list<textureNode*>& texture_list);

void textureMap(Object *obj, double a, double b, color *col);
void normalMap(Object *obj, double a, double b, struct point *n);

void alphaMap(Object *obj, double a, double b, double *set_alpha, double obj_alpha);
#endif