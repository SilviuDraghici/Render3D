#include "utils.h"
#include "mappings.h"
#include "objects.h"
#include "imageProcessor.h"

struct textureNode *texture_list;

void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list) {
   // Load a texture or normal map image from file and assign it to the
   // specified object.
   // type:   1  ->  Texture map  (RGB, .ppm)
   //         2  ->  Normal map   (RGB, .ppm)
   //         3  ->  Alpha map    (grayscale, .pgm)
   // Stores loaded images in a linked list to avoid replication
   struct image *im;
   struct textureNode *p;

   if (o != NULL) {
      // Check current linked list
      p = *(t_list);
      while (p != NULL) {
         if (strcmp(&p->name[0], filename) == 0) {
            // Found image already on the list
            if (type == 1)
               o->texImg = p->im;
            else if (type == 2)
               o->normalMap = p->im;
            else
               o->alphaMap = p->im;
            return;
         }
         p = p->next;
      }

      // Load this texture image
      if (type == 1 || type == 2)
         im = readPPMimage(filename);
      else if (type == 3)
         im = readPGMimage(filename);

      // Insert it into the texture list
      if (im != NULL) {
         p = (struct textureNode *)calloc(1, sizeof(struct textureNode));
         strcpy(&p->name[0], filename);
         p->type = type;
         p->im = im;
         p->next = NULL;
         // Insert into linked list
         if ((*(t_list)) == NULL)
            *(t_list) = p;
         else {
            p->next = (*(t_list))->next;
            (*(t_list))->next = p;
         }
         // Assign to object
         if (type == 1)
            o->texImg = im;
         else if (type == 2)
            o->normalMap = im;
         else
            o->alphaMap = im;
      }

   }  // end if (o != NULL)
}

void texMap(struct image *img, double a, double b, double *R, double *G, double *B) {
   /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
   a given object.

  The colour is returned in R, G, B. Uses bi-linear interpolation
  to determine texture colour.
 */

   a = MIN(1, MAX(0, a));
   b = MIN(1, MAX(0, b));
   a = a * (img->sx - 1);
   b = b * (img->sy - 1);

   int x1 = (int)floor(a);
   int y1 = (int)floor(b);
   int x2 = MIN(img->sy - 1, (int)ceil(a));
   int y2 = MIN(img->sy - 1, (int)ceil(b));
   //printf("a: %f b: %f\n", a, b);
   //printf("x: %d y: %d\n", x, y);
   double *rgbIm = (double *)img->rgbdata;
   *(R) = ((double)rgbIm[3 * (y1 * img->sx + x1) + 0]);
   *(G) = ((double)rgbIm[3 * (y1 * img->sx + x1) + 1]);
   *(B) = ((double)rgbIm[3 * (y1 * img->sx + x1) + 2]);
   //printf("r: %f g: %f b: %f\n", *R, *G, *B);
   //here
   return;
}

void texMapN(struct image *img, double a, double b, double *R, double *G, double *B) {
   /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
  a given object.

  The colour is returned in R, G, B. Uses nearest neighbour
  to determine texture colour.
 */

   a = MIN(1, MAX(0, a));
   b = MIN(1, MAX(0, b));
   //printf("a: %f b: %f\n", a, b);
   int x = (int)floor(a * (img->sx - 1));
   int y = (int)floor(b * (img->sy - 1));
   //printf("x: %d y: %d\n", x, y);
   double *rgbIm = (double *)img->rgbdata;
   *(R) = ((double)rgbIm[3 * (y * img->sx + x)]);
   *(G) = ((double)rgbIm[3 * (y * img->sx + x) + 1]);
   *(B) = ((double)rgbIm[3 * (y * img->sx + x) + 2]);
   //printf("r: %f g: %f b: %f\n", *R, *G, *B);
   //here
   return;
}

void alphaMap(struct image *img, double a, double b, double *alpha) {
   // Just like texture map but returns the alpha value at a,b,
   // notice that alpha maps are single layer grayscale images, hence
   // the separate function.

   *(alpha) = 1;  // Returns 1 which means fully opaque. Replace
   return;        // with your code if implementing alpha maps.
}