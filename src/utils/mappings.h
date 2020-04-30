#ifndef MAPPING_H
#define MAPPING_H

#include "color.h"

struct textureNode {
    char name[1024];
    int type;
    struct image *im;
    struct textureNode *next;
};

extern struct textureNode *texture_list;

void loadTexture(struct object *o, const char *filename, int type, struct textureNode **t_list);

void textureMap(struct object *obj, double a, double b, struct color *col);
void normalMap(struct object *obj, double a, double b, struct point *n);

void alphaMap(struct object *obj, double a, double b, double *set_alpha, double obj_alpha);
#endif