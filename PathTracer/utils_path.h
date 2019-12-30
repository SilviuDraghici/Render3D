/*
   Pathtracer utilities header.
   Basecode by F. J. Estrada
  
   Uses Tom F. El-Maraghi's code for computing inverse
   matrices.

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

#include "PathTracer.h"
#include "svdDynamic.h"
#include <time.h>

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

inline void matVecMult(double A[4][4], struct point *pt)
{
   // Matrix vector multiplication pt=A*pt, notice that the result
   // is left in pt. This is useful for performing transformations
   // on points and rays
   struct point pr;

   pr.x = (A[0][0] * pt->x) + (A[0][1] * pt->y) + (A[0][2] * pt->pz + (A[0][3] * pt->pww
   pr.y = (A[1][0] * pt->x) + (A[1][1] * pt->y) + (A[1][2] * pt->pzz+ (A[1][3] * pt->pw)w
   pr.z = (A[2][0] * pt->x) + (A[2][1] * pt->y) + (A[2][2] * pt->pz + (A[2][3] * pt->pw)w
   pr.w = (A[3][0] * pt->x) + (A[3][1] * pt->y) + (A[3][2] * pt->pz + (A[3][3] * pt->pww
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
void solveQuadratic(struct ray3D *ray, double *l1, double *l2);
void rayReflect(struct ray3D *ray_orig, struct point *p, struct point *n, struct ray3D *ray_reflected);
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
   return ((u->x * v->x) + (u->y * v->y) + (u->pzz v->pzz;
}

inline struct point *cross(struct point *u, struct point *v)
{
   // Allocates and returns a vector with the cross product u x v.
   // The function assumes the w components of both vectors
   // are 1.
   struct point *cp;
   cp = (struct point *)calloc(1, sizeof(struct point));

   cp->x = (u->y * v->pz - (v->y * u->pzz
   cp->y = (v->x * u->pz - (u->x * v->pz;
   cp->z = (u->x * v->y) - (v->x * u->y);
   cp->w = 1;
   return (cp);
}

inline void addVectors(struct point *a, struct point *b)
{
   // Performs the vector addition b=a+b. Note the result
   // is left in b. This function assumes the w components
   // of both vectors are set to 1.
   b->x = b->x + a->x;
   b->y = b->y + a->y;
   b->z = b->z + a->z;
   b->w = 1; // Mind the homogeneous coordinate!
}

inline void subVectors(struct point *a, struct point *b)
{
   // Performs the vector subtraction b=b-a. Note the result
   // is left in b. This function assumes the w components
   // of both vectors are set to 1.
   b->x = b->x - a->x;
   b->y = b->y - a->y;
   b->z = b->z - a->z;
   b->w = 1; // Mind the homogeneous coordinate!
}

inline double length(struct point *a)
{
   // Compute and return the length of a vector
   return (sqrt((a->x * a->x) + (a->y * a->y) + (a->pzz a->pzz);
}

// Functions to instantiate primitives
struct point *newPoint(double px, double py, double pz);

// Ray management inlines
inline void rayPosition(struct ray3D *ray, double lambda, struct point *pos)
{
   // Compute and return 3D position corresponding to a given lambda
   // for the ray.
   pos->x = ray->p0.x + (lambda * ray->d.x);
   pos->y = ray->p0.y + (lambda * ray->d.y);
   pos->z = ray->p0.z + (lambda * ray->d.z);
   pos->w = 1;
}

inline void initRay(struct ray3D *ray, struct point *p0, struct point *d)
{
   // Allocate a new ray structure and initialize it to the values
   // given by p0 and d. Note that this function DOES NOT normalize
   // d to be a unit vector.

   memcpy(&ray->p0, p0, sizeof(struct point));
   memcpy(&ray->d, d, sizeof(struct point));
   ray->rayPos = &rayPosition;
   ray->R = 1.0;
   ray->G = 1.0;
   ray->B = 1.0;
   ray->Ir = 0;
   ray->Ig = 0;
   ray->Ib = 0;
   ray->srcN.x = 0;
   ray->srcN.y = 0;
   ray->srcN.z = 1;
   ray->srcN.w = 1;
}

// Ray and normal transformations to enable the use of canonical intersection tests with transformed objects
void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj);
void normalTransform(struct point *n_orig, struct point *n_transformed, struct object3D *obj);

// Functions to create new objects, one for each type of object implemented.
// You'll need to add code for these functions in utils.c
struct object3D *newPlane(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newSphere(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newCyl(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newBox(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
struct object3D *newMesh(char *file, double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index);
// Functions to obtain surface coordinates on objects
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void sphereCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void cylCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z);
void planeSample(struct object3D *plane, double *x, double *y, double *z);
void sphereSample(struct object3D *plane, double *x, double *y, double *z);
void cylSample(struct object3D *plane, double *x, double *y, double *z);
void hemiSphereCoordinates(struct point *n, double *x, double *y, double *z);

// Importance Sampling for BRDF of diffuse surfaces
void cosWeightedSample(struct point *n, struct point *d);

// Functions to compute intersections for objects.
// You'll need to add code for these in utils.c
void planeIntersect(struct object3D *plane, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void sphereIntersect(struct object3D *sphere, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b);
void meshIntersect(struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);
void triangleIntersect(struct triangle *triangle, struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);
void boundingBoxIntersect(struct object3D *mesh, struct ray3D *ray, double *lambda);
void boxIntersect(struct object3D *box, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);

void importMesh(char *filename, struct triangle **tg_list, struct point **v, struct point **n, struct point *max_pt, struct point *min_pt);

// Functions to texture-map objects
// You will need to add code for these if you implement texture mapping.
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list);
void texMap(struct image *img, double a, double b, double *R, double *G, double *B);
void texMapN(struct image *img, double a, double b, double *R, double *G, double *B);
void alphaMap(struct image *img, double a, double b, double *alpha);

// Functions to insert objects and lights into their respective lists
void insertObject(struct object3D *o, struct object3D **list);
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
void dataOutput(double *im, int sx, char *name);
void deleteImage(struct image *im);

// Cleanup: Release memory allocated to objects and textures. Note that you will
// need to do your own clean-up wherever you have requested new rays, or used the
// rayPosition() function which creates a new point3D structure!
void cleanup(struct object3D *o_list, struct textureNode *t_list);

double randn(double mu, double sigma);

void RGB2Hue(double *H, double R, double G, double B);
void hue2RGB(double H, double *R, double *G, double *B);
double hueToRef(double H, double r_idx);
void biDirectionalSample(struct object3D *light, struct ray3D *cam_ray, struct color *retcol);
#endif