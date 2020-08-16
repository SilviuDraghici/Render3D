#ifndef MAPPING_H
#define MAPPING_H

#include "color.h"
#include "objects.h"

struct textureNode {
    char name[1024];
    int type;
    struct image *im;
    struct textureNode *next;
};

extern struct textureNode *texture_list;

void loadTexture(Object *o, const char *filename, int type, struct textureNode **t_list);

void textureMap(Object *obj, double a, double b, color *col);
void normalMap(Object *obj, double a, double b, struct point *n);

void alphaMap(Object *obj, double a, double b, double *set_alpha, double obj_alpha);
#endif