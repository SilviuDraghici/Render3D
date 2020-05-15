#include "affinetransforms.h"
#include "color.h"
#include "imageProcessor.h"

#ifndef OBJECTS_H
#define OBJECTS_H

struct RT_properties {
    double ambient;   // Ambient light albedo
    double diffuse;   // Diffuse component albedo
    double specular;  // Specular component albedo
    double global;    // Global component albedo

    double alpha;      // Opacity - if less than 1 this is a semi transparent object
    double shinyness;  // Exponent for the Phong specular component
};

struct PT_properties {
    double diffuse;  // Diffuse component proportion
    double reflect;  // Purely reflecting component proportion
    double refract;  // Refracting component proportion
                     // NOTE: diffuse+reflect+refract=1.0

    double LSweight;  // If this is an area light source, keeps track of its weight (volume)
};

class Object{
    public:
    
    char label[20];  // for debugging

    struct RT_properties rt;  // Object's albedos for Phong model (for ray tracing)
    struct PT_properties pt;  // Object's surface properties for Path tracing

    struct color col;    // Object's colour in RGB
    struct matrix T;     // T holds the transformation applied to this object.
    struct matrix Tinv;  // Tinv holds the inverse transformation

    struct image *texImg;     // Pointer to structure holding the texture for this object
    struct image *photonMap;  // Photon map for this object
    struct image *normalMap;  // Normal map for this object
    struct image *alphaMap;   // Alpha map for the object

    // Material properties
    double refl_sig;
    double r_index;  // Index of refraction

    int frontAndBack;   // Flag to indicate that both sides of the object
                        // should be lit.
    int isLightSource;  // Flag to indicate if this is an area light source

    Object(double r, double g, double b);

    void set_rayTrace_properties(double ambient, double diffuse, double specular, double global, double alpha, double shiny);
    void set_pathTrace_properties(double diffuse, double reflect, double refract);

    virtual void intersect(struct ray *r, double *lambda, struct point *p, struct point *n, double *a, double *b) = 0;
    virtual void surfaceCoordinates(double a, double b, double *x, double *y, double *z) = 0;
    virtual void randomPoint(double *x, double *y, double *z) = 0;

    Object *next = NULL;
};

class Plane : public Object{
    public:
    Plane(double r, double g, double b);
    void intersect(struct ray *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
    void surfaceCoordinates(double a, double b, double *x, double *y, double *z);
    void randomPoint(double *x, double *y, double *z);
};

class Sphere : public Object{
    public:
    using Object::Object;
    void intersect(struct ray *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
    void surfaceCoordinates(double a, double b, double *x, double *y, double *z);
    void randomPoint(double *x, double *y, double *z);
};

class Box : public Object{
    public:
    using Object::Object;
    void intersect(struct ray *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
    void surfaceCoordinates(double a, double b, double *x, double *y, double *z);
    void randomPoint(double *x, double *y, double *z);
};

/* The structure below defines a point light source */
struct pointLS {
    struct color col;      // Light source colour
    struct point p0;       // Light source location
    struct pointLS *next;  // Pointer to next light in the scene
};

void insertObject(Object *o, Object **list);

void normalTransform(struct point *n_orig, struct point *n_transformed, Object *obj);

struct pointLS *newPLS(struct point *p0, double r, double g, double b);
void insertPLS(struct pointLS *l, struct pointLS **list);
#endif