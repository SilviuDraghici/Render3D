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
inline void matMult(double A[4][4], double B[4][4])
{
   // Performs matrix multiplication B=A*B (notice the result is left
   // in B). This is so that we can compose transformations by
   // premultiplying a given transformation matrix by one of our
   // simple transformation matrices (rotate, translate, scale, etc).
   // Note the indexing convention is [row][col]

   double C[4][4];
   int i, j, k;

   memset(C, 0, 16 * sizeof(double));
   for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
         C[i][j] = (A[i][0] * B[0][j]) + (A[i][1] * B[1][j]) + (A[i][2] * B[2][j]) + (A[i][3] * B[3][j]);

   memcpy(B, C, 16 * sizeof(double));
}

inline void matVecMult(double A[4][4], struct point3D *pt)
{
   // Matrix vector multiplication pt=A*pt, notice that the result
   // is left in pt. This is useful for performing transformations
   // on points and rays
   struct point3D pr;

   pr.px = (A[0][0] * pt->px) + (A[0][1] * pt->py) + (A[0][2] * pt->pz) + (A[0][3] * pt->pw);
   pr.py = (A[1][0] * pt->px) + (A[1][1] * pt->py) + (A[1][2] * pt->pz) + (A[1][3] * pt->pw);
   pr.pz = (A[2][0] * pt->px) + (A[2][1] * pt->py) + (A[2][2] * pt->pz) + (A[2][3] * pt->pw);
   pr.pw = (A[3][0] * pt->px) + (A[3][1] * pt->py) + (A[3][2] * pt->pz) + (A[3][3] * pt->pw);
   memcpy(pt, &pr, 4 * sizeof(double));
}

// Matrix manipulation - Mind the fact that these functions will change the matrix T - for hierarchical objects you will
// need to be careful to make local copies for manipulation where needed.
void invert(double *T, double *Tinv);
void RotateXMat(double T[4][4], double theta);                      // X-axis rotation by theta radians
void RotateYMat(double T[4][4], double theta);                      // Y-axis rotation by theta radians
void RotateZMat(double T[4][4], double theta);                      // Z-axis rotation by theta radians
void TranslateMat(double T[4][4], double tx, double ty, double tz); // 3D translation
void ScaleMat(double T[4][4], double sx, double sy, double sz);     // 3D non-uniform scaling

// Functions for geometric manipulation of *INDIVIDUAL* objects - be careful with composite objects: You will have
// to apply any transforms to each of the components.
void RotateX(struct object3D *o, double theta);                      // X-axis rotation by theta radians
void RotateY(struct object3D *o, double theta);                      // Y-axis rotation by theta radians
void RotateZ(struct object3D *o, double theta);                      // Z-axis rotation by theta radians
void Translate(struct object3D *o, double tx, double ty, double tz); // 3D translation
void Scale(struct object3D *o, double sx, double sy, double sz);     // 3D non-uniform scaling
void printmatrix(double mat[4][4]);

// Vector management
inline void normalize(struct point3D *v)
{
   // Normalizes a vector to unit length.
   // Note that this assumes the w component of v is one, so make
   // sure you don't have homogeneous vectors with w != 1
   // floating around or things will go south.
   double l;
   l = v->px * v->px;
   l += (v->py * v->py);
   l += (v->pz * v->pz);
   l = 1.0 / sqrt(l);
   v->px *= l;
   v->py *= l;
   v->pz *= l;
}

inline double dot(struct point3D *u, struct point3D *v)
{
   // Computes the dot product of 3D vectors u and v.
   // The function assumes the w components of both vectors
   // are 1.
   return ((u->px * v->px) + (u->py * v->py) + (u->pz * v->pz));
}

inline struct point3D *cross(struct point3D *u, struct point3D *v)
{
   // Allocates and returns a vector with the cross product u x v.
   // The function assumes the w components of both vectors
   // are 1.
   struct point3D *cp;
   cp = (struct point3D *)calloc(1, sizeof(struct point3D));

   cp->px = (u->py * v->pz) - (v->py * u->pz);
   cp->py = (v->px * u->pz) - (u->px * v->pz);
   cp->pz = (u->px * v->py) - (v->px * u->py);
   cp->pw = 1;
   return (cp);
}

inline void addVectors(struct point3D *a, struct point3D *b)
{
   // Performs the vector addition b=a+b. Note the result
   // is left in b. This function assumes the w components
   // of both vectors are set to 1.
   b->px = b->px + a->px;
   b->py = b->py + a->py;
   b->pz = b->pz + a->pz;
   b->pw = 1; // Mind the homogeneous coordinate!
}

inline void subVectors(struct point3D *a, struct point3D *b)
{
   // Performs the vector subtraction b=b-a. Note the result
   // is left in b. This function assumes the w components
   // of both vectors are set to 1.
   b->px = b->px - a->px;
   b->py = b->py - a->py;
   b->pz = b->pz - a->pz;
   b->pw = 1; // Mind the homogeneous coordinate!
}

inline double length(struct point3D *a)
{
   // Compute and return the length of a vector
   return (sqrt((a->px * a->px) + (a->py * a->py) + (a->pz * a->pz)));
}

// Functions to instantiate primitives
struct point3D *newPoint(double px, double py, double pz);
struct pointLS *newPLS(struct point3D *p0, double r, double g, double b);

// Ray management inlines
inline void rayPosition(struct ray3D *ray, double lambda, struct point3D *pos)
{
   // Compute and return 3D position corresponding to a given lambda
   // for the ray.
   pos->px = ray->p0.px + (lambda * ray->d.px);
   pos->py = ray->p0.py + (lambda * ray->d.py);
   pos->pz = ray->p0.pz + (lambda * ray->d.pz);
   pos->pw = 1;
}

struct point3D path(double lambda);

inline void initRay(struct ray3D *ray, struct point3D *p0, struct point3D *d)
{
   // Initializes the given ray3D struct with the the position
   // and direction vectors. Note that this function DOES NOT normalize
   // d to be a unit vector.

   memcpy(&ray->p0, p0, sizeof(struct point3D));
   memcpy(&ray->d, d, sizeof(struct point3D));
   ray->rayPos = &rayPosition;
}

// Ray and normal transformations to enable the use of canonical intersection tests with transformed objects
void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj);
void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj);

void rayReflect(struct ray3D *ray_orig, struct point3D *p, struct point3D *n, struct ray3D *ray_reflected);

// Functions to create new objects, one for each type of object implemented.
// You'll need to add code for these functions in utils.c
struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newCyl(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny);
struct object3D *newCone(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny);
struct object3D *newMesh(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny, struct triangle *triangles, struct point3D *vertices);
struct object3D *newBox(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny);

void importMesh(struct triangle **tg_list, struct point3D **v);

void newTree(struct object3D **object_list, double hierarchyMat[4][4], struct color *col, double distFromC, double numBranches, double maxRotation, double depth);
void newBranch(struct object3D **object_list, double hierarchyMat[4][4], struct color *col, double numBranches, double maxRotation, double curr_depth, double depth);
void newFlower(struct object3D **object_list, struct color *petalCol, double hierarchyMat[4][4]);

// Functions to obtain surface coordinates on objects
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void sphereCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void cylCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void planeSample(struct object3D *plane, double *x, double *y, double *z);
void sphereSample(struct object3D *plane, double *x, double *y, double *z);
void cylSample(struct object3D *plane, double *x, double *y, double *z);

// Functions to compute intersections for objects.
// You'll need to add code for these in utils.c
void planeIntersect(struct object3D *plane, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void sphereIntersect(struct object3D *sphere, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void coneIntersect(struct object3D *cone, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void meshIntersect(struct object3D *mesh, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void triangleIntersect(struct triangle *triangle, struct object3D *mesh, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);
void boxIntersect(struct object3D *box, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);

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
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize);

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
