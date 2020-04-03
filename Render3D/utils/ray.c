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

void rayRefract(struct ray *ray_orig, struct object *obj, struct point *p, struct point *n, struct ray *ray_refracted) {

   double r_index = obj->r_index;
   double r, n1 = 1, n2 = 1, theta = -1;
   double c = dot(n, &(ray_orig->d));
   double R_Shlick, dice;
   // moving from outside object to inside the object
   if (c > 0)
   {
      *n *= -1;
      n1 = r_index;
   }
   else
   {
      n2 = r_index;
      c *= -1;
   }
   theta = c;
   r = n1 / n2;
   //theta = c;

   memcpy(ray_refracted, ray_orig, sizeof(struct ray));
   ray_refracted->p0 = *p - *n * THR;
   double s = 1 - (r * r) * (1 - (c * c));
   if (s > 0)
   {
      ray_refracted->d = ray_orig->d * r + *n * (r * c - sqrt(s));

      normalize(&ray_refracted->d);
      
      /*
      // Use Shlick's to figure out amount of reflected and refracted light
      double R0 = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2));
      // R(theta)=R0+((1-R0)*(1-cos(theta))^5)
      R_Shlick = MIN(1, R0 + (1 - R0) * pow(1 - theta, 5));
      dice = drand48();
      // randomly choose if refract or reflect on each iteration
      if (dice < (1 - R_Shlick))
      {
         dice = drand48();

         max_col = MAX(MAX(refractRay.R, refractRay.G), MAX(refractRay.R, refractRay.B));

         if (sqrt(dice) < max_col)
         {
            return PathTrace(&refractRay, depth + 1, col, obj, explt);
         }
      }
      */
   }
}