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

void rayRefract(struct ray *ray_orig, struct object *obj, struct point *p, struct point *n, struct ray *ray_refracted, double *s, double *R_Shlick) {
    double r_index = obj->r_index;
    double r, n1 = 1, n2 = 1, theta = -1;
    double c = dot(n, &(ray_orig->d));
    // moving from outside object to inside the object
    if (c > 0) {
        *n *= -1;
        n1 = r_index;
    } else {
        n2 = r_index;
        c *= -1;
    }
    theta = c;
    r = n1 / n2;
    //theta = c;

    memcpy(ray_refracted, ray_orig, sizeof(struct ray));
    ray_refracted->p0 = *p - *n * THR;
    *s = 1 - (r * r) * (1 - (c * c));

    // Use Shlick's to figure out amount of reflected and refracted light
    double R0 = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2));
    // R(theta)=R0+((1-R0)*(1-cos(theta))^5)
    *R_Shlick = MIN(1, R0 + (1 - R0) * pow(1 - theta, 5));
    ray_refracted->d = ray_orig->d * r + *n * (r * c - sqrt(*s));
}

void findFirstHit(struct ray *ray, double *lambda, struct object *Os, struct object **obj, struct point *p, struct point *n, double *a, double *b) {
    // Find the closest intersection between the ray and any objects in the scene.
    // Inputs:
    //   *ray    -  A pointer to the ray being traced
    //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
    //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
    //              projection
    // Outputs:
    //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
    //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
    //              this ray (this is required so you can do the shading)
    //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
    //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
    //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

    struct object *curr_obj = object_list;
    double curr_l, curr_a, curr_b;
    struct point curr_p, curr_n;
    *lambda = INFINITY;

    while (curr_obj != NULL) {
        curr_obj->intersect(curr_obj, ray, &curr_l, &curr_p, &curr_n, &curr_a, &curr_b);
        if (THR < curr_l && curr_l < *lambda) {
            *lambda = curr_l;
            *obj = curr_obj;
            *p = curr_p;
            *n = curr_n;
            *a = curr_a;
            *b = curr_b;
        }
        curr_obj = curr_obj->next;
    }
}
