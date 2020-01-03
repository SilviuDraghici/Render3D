#ifndef MAPPING_H
#define MAPPING_H

struct textureNode {
   char name[1024];
   int type;
   struct image *im;
   struct textureNode *next;
};

void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list);
void texMap(struct image *img, double a, double b, double *R, double *G, double *B);
void texMapN(struct image *img, double a, double b, double *R, double *G, double *B);
void alphaMap(struct image *img, double a, double b, double *R, double *G, double *B);
#endif