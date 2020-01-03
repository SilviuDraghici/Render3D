#include "ray.h"

#include "affineTransforms.h"
#include "objects.h"

void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj) {
   // Transforms a ray using the inverse transform for the specified object. This is so that we can
   // use the intersection test for the canonical object. Note that this has to be done carefully!

   //copy original ray to new ray
   memcpy(ray_transformed, ray_orig, sizeof(struct ray3D));

   //inverse tranform ray origin
   ray_transformed->p0 = obj->Tinv * ray_transformed->p0;

   //inverse tranform ray direction without translation
   ray_transformed->d.w = 0;
   ray_transformed->d = obj->Tinv * ray_transformed->d;
   ray_transformed->d.w = 1;
}

void rayPosition(struct ray3D *ray, double lambda, struct point *pos) {
   // Compute and return 3D position corresponding to a given lambda
   // for the ray.
   *pos = ray->p0 + (ray->d * lambda);
}