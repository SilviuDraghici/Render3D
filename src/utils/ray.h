#include "color.h"
#include "utils.h"
#include "objects.h"

#ifndef RAY_H
#define RAY_H

struct ray_pathTrace_elements {
    // this stuct contains ray elements used only in Path Tracing
    struct color ray_col;  // Keeps track of the ray's colour as it
                           // bounces through the scene

    struct color expl_col;  // Accumulators of brightness for explicit
                            // light sampling

    int isLightRay = 0;
};

/* The structure below defines a ray in 3D homogeneous coordinates,
   the point corresponds to the representation r(\lambda)=p+(\lambda)d */
struct ray {
    struct point p0;  // Ray origin (at lambda=0)
    struct point d;   // Ray direction

    struct ray_pathTrace_elements pt;
};

void rayTransform(struct ray *ray_orig, struct ray *ray_transformed, Object *obj);
void rayPosition(struct ray *ray, double lambda, struct point *pos);
void rayReflect(struct ray *ray_orig, struct point *p, struct point *n, struct ray *ray_reflected);
void rayRefract(struct ray *ray_orig, Object *obj, struct point *p, struct point *n, struct ray *ray_refracted, double *s, double *R_Shlick);
void findFirstHit(struct ray *ray, double *lambda, Object *Os, Object **obj, struct point *p, struct point *n, double *a, double *b);

#endif