/*
   utils.h
   Basecode by F. J. Estrada
  
   Uses Tom F. El-Maraghi's code for computing inverse
   matrices.

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

#include "RayTracer.h"
#include "svdDynamic.h"
#include "modes.h"

#ifndef __utils_header
#define __utils_header

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define THR 0.0001
// Functions to apply transformations to objects.
// If you add any transformations to the list below, document them carefully

// Matrix manipulation - Mind the fact that these functions will change the matrix T - for hierarchical objects you will
// need to be careful to make local copies for manipulation where needed.
void invert(double *T, double *Tinv);

void printmatrix(double mat[4][4]);

// Vector management
inline void normalize(struct point *v)
{
   // Normalizes a vector to unit length.
   // Note that this assumes the w component of v is one, so make
   // sure you don't have homogeneous vectors with w != 1
   // floating around or things will go south.
   double l;
   l = v->x * v->x;
   l += (v->y * v->y);
   l += (v->z * v->z);
   l = 1.0 / sqrt(l);
   v->x *= l;
   v->y *= l;
   v->z *= l;
}

inline double dot(struct point *u, struct point *v)
{
   // Computes the dot product of 3D vectors u and v.
   // The function assumes the w components of both vectors
   // are 1.
   return ((u->x * v->x) + (u->y * v->y) + (u->z * v->z));
}

inline struct point *cross(struct point *u, struct point *v)
{
   // Allocates and returns a vector with the cross product u x v.
   // The function assumes the w components of both vectors
   // are 1.
   struct point *cp;
   cp = (struct point *)calloc(1, sizeof(struct point));

   cp->x = (u->y * v->z) - (v->y * u->z);
   cp->y = (v->x * u->z) - (u->x * v->z);
   cp->z = (u->x * v->y) - (v->x * u->y);
   cp->w = 1;
   return (cp);
}

inline double length(struct point *a)
{
   // Compute and return the length of a vector
   return (sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z)));
}

// Functions to instantiate primitives
struct point *newPoint(double px, double py, double pz);
struct pointLS *newPLS(struct point *p0, double r, double g, double b);

// Ray management inlines
inline void rayPosition(struct ray3D *ray, double lambda, struct point *pos)
{
   // Compute and return 3D position corresponding to a given lambda
   // for the ray.
   *pos = ray->p0 + (ray->d * lambda);
}

struct point path(double lambda);

inline void initRay(struct ray3D *ray, struct point *p0, struct point *d)
{
   // Initializes the given ray3D struct with the the position
   // and direction vectors. Note that this function DOES NOT normalize
   // d to be a unit vector.

   memcpy(&ray->p0, p0, sizeof(struct point));
   memcpy(&ray->d, d, sizeof(struct point));
}

// Ray and normal transformations to enable the use of canonical intersection tests with transformed objects
void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj);
void normalTransform(struct point *n_orig, struct point *n_transformed, struct object3D *obj);

void rayReflect(struct ray3D *ray_orig, struct point *p, struct point *n, struct ray3D *ray_reflected);

// Functions to create new objects, one for each type of object implemented.
// You'll need to add code for these functions in utils.c
struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newCyl(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newCone(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny);
struct object3D *newMesh(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny, struct triangle *triangles, struct point *vertices);
struct object3D *newBox(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny);

void importMesh(struct triangle **tg_list, struct point **v);

void newTree(struct object3D **object_list, struct transform hierarchy, struct color *col, double distFromC, double numBranches, double maxRotation, double depth);
void newBranch(struct object3D **object_list, struct transform hierarchy, struct color *col, double numBranches, double maxRotation, double curr_depth, double depth);
void newFlower(struct object3D **object_list, struct color *petalCol, struct transform hierarchy);

// Functions to obtain surface coordinates on objects
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void sphereCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void cylCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void planeSample(struct object3D *plane, double *x, double *y, double *z);
void sphereSample(struct object3D *plane, double *x, double *y, double *z);
void cylSample(struct object3D *plane, double *x, double *y, double *z);

// Functions to compute intersections for objects.
// You'll need to add code for these in utils.c
void planeIntersect(struct object3D *plane, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void sphereIntersect(struct object3D *sphere, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void coneIntersect(struct object3D *cone, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);
void meshIntersect(struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);
void triangleIntersect(struct triangle *triangle, struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);
void boxIntersect(struct object3D *box, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);

// Functions to texture-map objects
// You will need to add code for these if you implement texture mapping.
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list);
void texMap(struct image *img, double a, double b, double *R, double *G, double *B);
void texMapN(struct image *img, double a, double b, double *R, double *G, double *B);
void alphaMap(struct image *img, double a, double b, double *R, double *G, double *B);

// Functions to insert objects and lights into their respective lists
void insertObject(struct object3D *o, struct object3D **list);
void insertPLS(struct pointLS *l, struct pointLS **list);
void addAreaLight(double sx, double sy, double nx, double ny, double nz,
                  double tx, double ty, double tz, int lx, int ly,
                  double r, double g, double b, struct object3D **o_list, struct pointLS **l_list);

// Function to set up the camera and viewing coordinate frame.
// You will have to add code to this function's body in utils.c
struct view *setupView(struct point *e, struct point *g, struct point *up, double f, double wl, double wt, double wsize);

// Image management output. Note that you will need to free() any images you
// allocate with newImage() using deleteImage().
struct image *readPPMimage(const char *filename);
struct image *readPGMimage(const char *filename);
struct image *newImage(int size_x, int size_y);
void imageOutput(struct image *im, const char *filename);
void deleteImage(struct image *im);

// Cleanup: Release memory allocated to objects and light sources. Note that you will
// need to do your own clean-up wherever you have requested ray positions,
// since each call to the ray position function returns a newly allocated point3D
// structure.
void cleanup(struct object3D *o_list, struct pointLS *l_list, struct textureNode *t_list);

#endif
