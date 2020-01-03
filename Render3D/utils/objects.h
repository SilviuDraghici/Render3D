#include "imageProcessor.h"
#include "affinetransforms.h"
#include "color.h"

#ifndef OBJECTS_H
#define OBJECTS_H
struct object3D {
   struct albedosPhong alb;  // Object's albedos for Phong model
   struct color col;         // Object's colour in RGB
   struct transform T;       // T holds the transformation applied to this object.
   struct transform Tinv;    // Tinv holds the inverse transformation

   char label;

   // Below we set up space for a pointer to the intersection function for this object.
   // Note that the intersection function must compute the lambda at the intersection, the
   // intersection point p, the normal at that point n, and the texture coordinates (a,b).
   // The texture coordinates are not used unless texImg!=NULL and a textureMap function
   // has been provided
   void (*intersect)(struct object3D *obj, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);

   // Texture mapping function. Takes normalized texture coordinates (a,b) and returns the
   // texture colour at that point using bi-linear interpolation
   void (*textureMap)(struct image *img, double a, double b, double *R, double *G, double *B);

   // Functions to return coordinates on the surface of the object. One takes as input the a and b
   // parameters for the parametric function of the object and returns the (x,y,z) coordinates
   // on the object surface. The second returns a uniformly random-sampled point on the surface.
   // These are needed for Photon Mapping.
   void (*surfaceCoords)(struct object3D *obj, double a, double b, double *x, double *y, double *z);
   void (*randomPoint)(struct object3D *obj, double *x, double *y, double *z);
   struct image *texImg;     // Pointer to structure holding the texture for this object
   struct image *photonMap;  // Photon map for this object
   struct image *normalMap;  // Normal map for this object
   struct image *alphaMap;   // Alpha map for the object

   // Material properties
   double alpha;       // Opacity - if less than 1 this is a semi transparent object and refraction rays
                       // should be implemented
   double r_index;     // Index of refraction
   double shinyness;   // Exponent for the Phong specular component
   int frontAndBack;   // Flag to indicate that both sides of the object
                       // should be lit.
   int isLightSource;  // Flag to indicate if this is an area light source
   int isCSG;          // Object is part of a CSG composite object. Links to components via CSGnext
   int photonMapped;   // This object accumulates photons under photon mapping
   int normalMapped;   // This object has an associated normal map
   int alphaMapped;    // This object has an associated alpha map

   struct object3D *next;     // Pointer to next entry in object linked list
};

struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
void planeIntersect(struct object3D *plane, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void planeSample(struct object3D *plane, double *x, double *y, double *z);

struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newCyl(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newCone(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny);
struct object3D *newBox(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny);

void insertObject(struct object3D *o, struct object3D **list);

void normalTransform(struct point *n_orig, struct point *n_transformed, struct object3D *obj);

#endif