#include "utils.h" 

#ifndef RAY_H
#define RAY_H
/* The structure below defines a ray in 3D homogeneous coordinates,
   the point corresponds to the representation r(\lambda)=p+(\lambda)d */
struct ray {
   struct point p0;  // Ray origin (at lambda=0)
   struct point d;   // Ray direction

   /* You may add data here to keep track of any values associated */
   /* with this ray when implementing advanced raytracing features */
};

void rayTransform(struct ray *ray_orig, struct ray *ray_transformed, struct object *obj);
void rayPosition(struct ray *ray, double lambda, struct point *pos);
void rayReflect(struct ray *ray_orig, struct point *p, struct point *n, struct ray *ray_reflected);
#endif