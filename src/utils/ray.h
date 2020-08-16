#include "color.h"
#include "utils.h"
#include "objects.h"
#include "buildscene.h"

#ifndef RAY_H
#define RAY_H

struct ray_pathTrace_elements {
    // this stuct contains ray elements used only in Path Tracing
    color ray_col;  // Keeps track of the ray's colour as it
                           // bounces through the scene

    color expl_col;  // Accumulators of brightness for explicit
                            // light sampling

    int isLightRay = 0;
};

/* The structure below defines a ray in 3D homogeneous coordinates,
   the point corresponds to the representation r(\lambda)=p+(\lambda)d */
struct Ray {
    point p0;  // Ray origin (at lambda=0)
    point d;   // Ray direction

    ray_pathTrace_elements pt;
};

void rayTransform(Ray *ray_orig, Ray *ray_transformed, Object *obj);
void rayPosition(Ray *ray, double lambda, point *pos);
void rayReflect(Ray *ray_orig, point *p, point *n, Ray *ray_reflected);
void rayRefract(Ray *ray_orig, Object *obj, point *p, point *n, Ray *ray_refracted, double *s, double *R_Shlick);
void findFirstHit(Scene *scene, Ray *ray, double *lambda, Object *Os, Object **obj, point *p, point *n, double *a, double *b);

#endif