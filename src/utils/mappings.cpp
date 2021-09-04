#include "mappings.h"

#include <algorithm>

#include "imageProcessor.h"
#include "objects.h"
#include "utils.h"

auto get_extension(std::string_view fn) {
  return fn.substr(fn.find_last_of(".") + 1);
}

inline void uvMap(struct image *img, double u, double v, color *col) {
    int x = (unsigned int)(u * img->sx) % img->sx;
    int y = (unsigned int)(v * img->sy) % img->sy;

    double *rgbIm = (double *)img->rgbdata;
    col->R = ((double)rgbIm[3 * (y * img->sx + x) + 0]);
    col->G = ((double)rgbIm[3 * (y * img->sx + x) + 1]);
    col->B = ((double)rgbIm[3 * (y * img->sx + x) + 2]);
}

void setTexture(Object *o, const std::string& filename, int type, std::list<textureNode*>& texture_list) {
    // Load a texture or normal map image from file and assign it to the
    // specified object.
    // type:   1  ->  Texture map  (RGB, .ppm)
    //         2  ->  Normal map   (RGB, .ppm)
    //         3  ->  Alpha map    (grayscale, .pgm)
    // Stores loaded images in a linked list to avoid replication
    image *im;
    textureNode *texture_node = loadTexture(filename, type, texture_list);

    // Found image already on the list
    if (type == 1)
        o->texImg = texture_node->im;
    else if (type == 2)
        o->normalMap = texture_node->im;
    else
        o->alphaMap = texture_node->im;
}

textureNode* loadTexture(const std::string& filename, int type, std::list<textureNode*>& texture_list){
    image *im;
    textureNode* texture_node = *std::find(texture_list.begin(), texture_list.end(), filename);;
    if (texture_node != *texture_list.end()) {
        // Found image already on the list
        return texture_node;
    }
    // Load this texture image
    if (type == 1 || type == 2){
        auto extension = get_extension(filename);
        if (extension == "ppm") {
            im = readPPMimage(filename.c_str());
        } else if (extension == "png") {
            im = readPNGimage(filename.c_str());
        }
    } else if (type == 3)
        im = readPGMimage(filename.c_str());

    // Insert it into the texture list
    if (im != NULL) {
        texture_node = new textureNode;
        texture_node->name = filename;
        texture_node->im = im;
        // Insert into linked list
        texture_list.push_front(texture_node);
    } else {
        texture_node = NULL;
    }
    return texture_node;
}

void textureMap(Object *obj, double a, double b, color *col) {
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
        if (obj->isLightSource){
            col->R *= obj->col.R;
            col->G *= obj->col.G;
            col->B *= obj->col.B;
        }
    }
    return;
}

void normalMap(Object *obj, double a, double b, struct point *n) {
    if (obj->normalMap != NULL) {
        //printf("normal mapping \n");
        color delta;
        uvMap(obj->normalMap, a, b, &delta);
        n->x -= (2 * delta.R - 1);
        n->y -= (2 * delta.G - 1);
        n->z -= (2 * delta.B - 1);
        normalize(n);
        n->w = 1;
    }
}

void alphaMap(Object *obj, double a, double b, double *set_alpha, double obj_alpha) {
    // Just like texture map but returns the alpha value at a,b,
    // notice that alpha maps are single layer grayscale images, hence
    // the separate function.

    if (obj->alphaMap == NULL) {
        *set_alpha = obj_alpha;
    } else {
        struct image *img = obj->alphaMap;
        double xc1, xc2;

        a = MIN(1, MAX(0, a));
        b = MIN(1, MAX(0, b));
        a = a * (img->sx - 1);
        b = b * (img->sy - 1);

        int x1 = (int)floor(a);
        int y1 = (int)floor(b);
        int x2 = MIN(img->sx - 1, (int)ceil(a));
        int y2 = MIN(img->sy - 1, (int)ceil(b));

        double *rgbIm = (double *)img->rgbdata;
        double ax1, ax2, ay1, ay2;
        ax1 = (x2 - a) / (x2 - x1);
        ax2 = (a - x1) / (x2 - x1);
        xc1 = ax1 * ((double)rgbIm[(y1 * img->sx + x1)]) + ax2 * ((double)rgbIm[(y1 * img->sx + x2)]);
        xc2 = ax1 * ((double)rgbIm[(y2 * img->sx + x1)]) + ax2 * ((double)rgbIm[(y2 * img->sx + x2)]);

        //printf("a: %f b: %f\n", a, b);
        //printf("x: %d y: %d\n", x, y);
        ay1 = (y2 - b) / (y2 - y1);
        ay2 = (b - y1) / (y2 - y1);
        //printf("r: %f g: %f b: %f\n", *R, *G, *B);
        //here

        *set_alpha = ay1 * xc1 + ay2 * xc2;
    }
}