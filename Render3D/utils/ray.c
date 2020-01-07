#include "ray.h"

#include "affineTransforms.h"
#include "objects.h"

void rayTransform(struct ray *ray_orig, struct ray *ray_transformed, struct object *obj) {
   // Transforms a ray using the inverse transform for the specified object. This is so that we can
   // use the intersection test for the canonical object. Note that this has to be done carefully!

   //copy original ray to new ray
   memcpy(ray_transformed, ray_orig, sizeof(struct ray));

   //inverse tranform ray origin
   ray_transformed->p0 = obj->Tinv * ray_transformed->p0;

   //inverse tranform ray direction without translation
   ray_transformed->d.w = 0;
   ray_transformed->d = obj->Tinv * ray_transformed->d;
   ray_transformed->d.w = 1;
}

void rayPosition(struct ray *ray, double lambda, struct point *pos) {
   // Compute and return 3D position corresponding to a given lambda
   // for the ray.
   *pos = ray->p0 + (ray->d * lambda);
}

void rayReflect(struct ray *ray_orig, struct point *p, struct point *n, struct ray *ray_reflected) {
   //this function assumes n is unit length!

   //reflection starts at point of intersection
   memcpy(&(ray_reflected->p0), p, sizeof(struct point));

   //r=d−2(d⋅n)n
   double ddotn = dot(&(ray_orig->d), n);
   ray_reflected->d = ray_orig->d - *n * 2 * ddotn;
   //normalize(&ray_reflected->d);
}