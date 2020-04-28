#include "imageProcessor.h"
#include "affinetransforms.h"
#include "color.h"

#ifndef OBJECTS_H
#define OBJECTS_H

extern struct object *object_list;
extern struct pointLS *light_list;

struct RT_properties
{
   double ambient; // Ambient light albedo
   double diffuse; // Diffuse component albedo
   double specular; // Specular component albedo
   double global; // Global component albedo

   double alpha; // Opacity - if less than 1 this is a semi transparent object
   double shinyness;   // Exponent for the Phong specular component
};

struct PT_properties{
    double diffuse; // Diffuse component proportion
    double reflect; // Purely reflecting component proportion
    double refract; // Refracting component proportion
                   // NOTE: diffuse+reflect+refract=1.0

    double LSweight;   // If this is an area light source, keeps track of its weight (volume)
};


struct object {
   char label[20]; // for debugging

   struct RT_properties rt;  // Object's albedos for Phong model (for ray tracing)
   struct PT_properties pt; // Object's surface properties for Path tracing

   struct color col;         // Object's colour in RGB
   struct matrix T;       // T holds the transformation applied to this object.
   struct matrix Tinv;    // Tinv holds the inverse transformation

   // Below we set up space for a pointer to the intersection function for this object.
   // Note that the intersection function must compute the lambda at the intersection, the
   // intersection point p, the normal at that point n, and the texture coordinates (a,b).
   // The texture coordinates are not used unless texImg!=NULL and a textureMap function
   // has been provided
   void (*intersect)(struct object *obj, struct ray *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);

   // Functions to return coordinates on the surface of the object. One takes as input the a and b
   // parameters for the parametric function of the object and returns the (x,y,z) coordinates
   // on the object surface. The second returns a uniformly random-sampled point on the surface.
   // These are needed for Photon Mapping.
   void (*surfaceCoords)(struct object *obj, double a, double b, double *x, double *y, double *z);
   void (*randomPoint)(struct object *obj, double *x, double *y, double *z);
   struct image *texImg;     // Pointer to structure holding the texture for this object
   struct image *photonMap;  // Photon map for this object
   struct image *normalMap;  // Normal map for this object
   struct image *alphaMap;   // Alpha map for the object

   // Material properties
   double refl_sig;
   double r_index;     // Index of refraction
   
   int frontAndBack;   // Flag to indicate that both sides of the object
                       // should be lit.
   int isLightSource;  // Flag to indicate if this is an area light source

   struct object *next;     // Pointer to next entry in object linked list
};

/* The structure below defines a point light source */
struct pointLS
{
   struct color col;     // Light source colour
   struct point p0;      // Light source location
   struct pointLS *next; // Pointer to next light in the scene
};

extern struct object *object_list;
extern struct pointLS *light_list;

struct object *newPlane(double r, double g, double b);
void planeIntersect(struct object *plane, struct ray *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void planeCoordinates(struct object *plane, double a, double b, double *x, double *y, double *z);
void planeSample(struct object *plane, double *x, double *y, double *z);

struct object *newSphere(double r, double g, double b);
void sphereIntersect(struct object *sphere, struct ray *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);
void sphereCoordinates(struct object *plane, double a, double b, double *x, double *y, double *z);
void sphereSample(struct object *plane, double *x, double *y, double *z);

struct object *newCyl(double r, double g, double b);


struct object *newCone(double r, double g, double b);


struct object *newBox(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny);

void set_rayTrace_properties(struct object *o, double ambient, double diffuse, double specular, double global, double alpha, double shiny);
void set_pathTrace_properties(struct object *o, double diffuse, double reflect, double refract);

void insertObject(struct object *o, struct object **list);

void normalTransform(struct point *n_orig, struct point *n_transformed, struct object *obj);

struct pointLS *newPLS(struct point *p0, double r, double g, double b);
void insertPLS(struct pointLS *l, struct pointLS **list);
#endif