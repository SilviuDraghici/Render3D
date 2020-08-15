#include <string.h>

#include <iostream>

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

    double alpha;  // Opacity - if less than 1 this is a semi transparent object
    double shinyness;  // Exponent for the Phong specular component
};

struct PT_properties {
    double diffuse;  // Diffuse component proportion
    double reflect;  // Purely reflecting component proportion
    double refract;  // Refracting component proportion
                     // NOTE: diffuse+reflect+refract=1.0

    double LSweight;  // If this is an area light source, keeps track of its
                      // weight (volume)
};

class Bounds {
   public:
    point min, max;
    Bounds(){
        min.x = min.y = min.z = INFINITY;
        max.x = max.y = max.z = -INFINITY;
    }
    Axis longestAxis() const;
    double offset(const point &p, Axis dim) const;
    double surface_area() const;
};

void union_bounds(Bounds &a, point &b, Bounds &union_box);
void union_bounds(Bounds &a, Bounds &b, Bounds &union_box);

class Primitive {
   public:
    virtual double intersect(struct ray *r, double lambda) = 0;
    virtual bool isprim() = 0;
    virtual double min_x() const = 0;
    virtual double min_y() const = 0;
    virtual double min_z() const = 0;
    virtual double max_x() const = 0;
    virtual double max_y() const = 0;
    virtual double max_z() const = 0;
    virtual Axis longestAxis() const;
};

class Object : public Primitive{
   public:
    char label[20];  // for debugging

    struct RT_properties
        rt;  // Object's albedos for Phong model (for ray tracing)
    struct PT_properties pt;  // Object's surface properties for Path tracing

    struct color col;    // Object's colour in RGB
    struct matrix T;     // T holds the transformation applied to this object.
    struct matrix Tinv;  // Tinv holds the inverse transformation

    struct image
        *texImg;  // Pointer to structure holding the texture for this object
    struct image *normalMap;  // Normal map for this object
    struct image *alphaMap;   // Alpha map for the object

    // Material properties
    double refl_sig;
    double r_index;  // Index of refraction

    int frontAndBack;   // Flag to indicate that both sides of the object
                        // should be lit.
    int isLightSource;  // Flag to indicate if this is an area light source

    Bounds w_bound; // the bounds for this object in world coordinates

    Object(double r, double g, double b);
    void set_color(double r, double g, double b);
    void set_rayTrace_properties(double ambient, double diffuse,
                                 double specular, double global, double alpha,
                                 double shiny);
    void set_pathTrace_properties(double diffuse, double reflect,
                                  double refract);

    virtual void set_canonical_bounds();
    void invert_and_bound(); // calculates inverse matrix and sets w_bound

    virtual void intersect(struct ray *r, double *lambda, struct point *p,
                           struct point *n, double *a, double *b) = 0;

    virtual void surfaceCoordinates(double a, double b, double *x, double *y,
                                    double *z);
    virtual void randomPoint(double *x, double *y, double *z);

    double intersect(struct ray *r, double lambda);
    bool isprim();
    double min_x() const;
    double min_y() const;
    double min_z() const;
    double max_x() const;
    double max_y() const;
    double max_z() const;

    Object *next = NULL;
};

class Plane : public Object {
   public:
    Plane(double r, double g, double b);
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
    void set_canonical_bounds();
    void surfaceCoordinates(double a, double b, double *x, double *y,
                            double *z);
    void randomPoint(double *x, double *y, double *z);
};

class Sphere : public Object {
   public:
    using Object::Object;
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
    void set_canonical_bounds();
    void surfaceCoordinates(double a, double b, double *x, double *y,
                            double *z);
    void randomPoint(double *x, double *y, double *z);
};

class Box : public Object {
   public:
    using Object::Object;
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
};

class Triangle : public Object {
   public:
    Triangle(double r, double g, double b);
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);

    void setPoints(double x1, double y1, double z1, double x2, double y2,
                   double z2, double x3, double y3, double z3);

    void setPoints(point p1, point p2, point p3);

   private:
    point p1, p2, p3;
};

class Polygon : public Object {
   public:
    Polygon(double r, double g, double b);
    void intersect(struct ray *r, double *lambda, struct point *p,
                   struct point *n, double *a, double *b);
    void setNumPoints(int num);
    void addPoint(point point);
    void addPoint(double x, double y, double z);

   private:
    int numPoints = 3;
    int currPoint = 0;  // current point used for adding points
    point *p;           // pointer to array of vertices
    point *e;  // pointer to array of edge vectors used for intersection test
    void calculate_edge_vectors();
};

/* The structure below defines a point light source */
struct pointLS {
    struct color col;      // Light source colour
    struct point p0;       // Light source location
    struct pointLS *next;  // Pointer to next light in the scene
};

void normalTransform(struct point *n_orig, struct point *n_transformed,
                     Object *obj);

struct pointLS *newPLS(struct point *p0, double r, double g, double b);
void insertPLS(struct pointLS *l, struct pointLS **list);
#endif