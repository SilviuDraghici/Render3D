#include "mappings.h"
#include "imageProcessor.h"
#include "objects.h"
#include "utils.h"

struct textureNode *texture_list;

inline void uvMap(struct image *img, double u, double v, struct color *col) {
    u = MIN(1, MAX(0, u));
    v = MIN(1, MAX(0, v));
    u = u * (img->sx - 1);
    v = v * (img->sy - 1);

    int x1 = (int)floor(u);
    int y1 = (int)floor(v);
    int x2 = MIN(img->sy - 1, (int)ceil(u));
    int y2 = MIN(img->sy - 1, (int)ceil(v));
    //printf("u: %f v: %f\n", u, v);
    //printf("x: %d y: %d\n", x, y);
    double *rgbIm = (double *)img->rgbdata;
    col->R = ((double)rgbIm[3 * (y1 * img->sx + x1) + 0]);
    col->G = ((double)rgbIm[3 * (y1 * img->sx + x1) + 1]);
    col->B = ((double)rgbIm[3 * (y1 * img->sx + x1) + 2]);
    //printf("r: %f g: %f b: %f\n", *R, *G, *B);
}

void loadTexture(struct object *o, const char *filename, int type, struct textureNode **t_list) {
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

void textureMap(struct object *obj, double a, double b, struct color *col) {
    /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
   a given object.

  The colour is returned in R, G, B. Uses bi-linear interpolation
  to determine texture colour.
 */
    if (obj->texImg == NULL)  // Not textured, use object colour
    {
        *col = obj->col;
    } else {
        // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
        // for the object. Note that we will use textures also for Photon Mapping.
        uvMap(obj->texImg, a, b, col);
    }
    return;
}

void normalMap(struct object *obj, double a, double b, struct point *n) {
    if (obj->normalMap != NULL) {
        //printf("normal mapping \n");
        struct color delta;
        uvMap(obj->normalMap, a, b, &delta);
        n->x -= (2 * delta.R - 1);
        n->y -= (2 * delta.G - 1);
        n->z -= (2 * delta.B - 1);
        normalize(n);
        n->w = 1;
    }
}

void alphaMap(struct object *obj, double a, double b, double *set_alpha, double obj_alpha) {
    // Just like texture map but returns the alpha value at a,b,
    // notice that alpha maps are single layer grayscale images, hence
    // the separate function.
    
    if (obj->alphaMap == NULL) {
        *set_alpha = obj_alpha;
    } else {
        struct image *img = obj->alphaMap;
        double xc1,xc2;

        a = MIN(1, MAX(0, a));
        b = MIN(1, MAX(0, b));
        a = a*(img->sx - 1);
        b = b*(img->sy - 1);

        int x1 = (int)floor(a);
        int y1 = (int)floor(b);
        int x2 = MIN(img->sx -1, (int)ceil(a));
        int y2 = MIN(img->sy -1, (int)ceil(b));

        double *rgbIm=(double *)img->rgbdata;
        double ax1, ax2, ay1, ay2;
        ax1 = (x2 - a)/(x2 - x1);
        ax2 = (a - x1)/(x2 - x1);
        xc1 = ax1*((double)rgbIm[(y1*img->sx + x1)])  + ax2*((double)rgbIm[(y1*img->sx + x2)]);
        xc2 = ax1*((double)rgbIm[(y2*img->sx + x1)])  + ax2*((double)rgbIm[(y2*img->sx + x2)]);


        //printf("a: %f b: %f\n", a, b);
        //printf("x: %d y: %d\n", x, y);
        ay1 = (y2 - b)/(y2 - y1);
        ay2 = (b - y1)/(y2 - y1);
        //printf("r: %f g: %f b: %f\n", *R, *G, *B);
        //here

        *set_alpha = ay1*xc1 + ay2*xc2;
    }
}