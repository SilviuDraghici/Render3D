#include "objects.h"

#include <stdio.h>
#include <stdlib.h>

#include "mappings.h"
#include "ray.h"
#include "utils.h"

Object *object_list;
struct pointLS *light_list;

Object::Object(double r = 1, double g = 1, double b = 1) {
    col.R = r;
    col.G = g;
    col.B = b;
    rt.alpha = 1;
    rt.shinyness = 2;
    pt.LSweight = 1;
    r_index = 1;
    refl_sig = 0;
    texImg = NULL;
    photonMap = NULL;
    normalMap = NULL;
    alphaMap = NULL;
    frontAndBack = 0;
    isLightSource = 0;
    T = I();
    next = NULL;
}

void Object::set_rayTrace_properties(double ambient, double diffuse, double specular, double global, double alpha, double shiny) {
    rt.ambient = ambient;
    rt.diffuse = diffuse;
    rt.specular = specular;
    rt.global = global;
    rt.alpha = alpha;
    rt.shinyness = shiny;

    // create a similar set of params incase the object is drawn in pathtrace mode
    pt.diffuse = alpha * (ambient + diffuse);
    pt.reflect = alpha * (1 - diffuse);
    pt.refract = 1 - alpha;

    // ensure sum to 1
    double sum = pt.diffuse + pt.reflect + pt.refract;
    pt.diffuse /= sum;
    pt.reflect /= sum;
    pt.refract /= sum;
}

void Object::set_pathTrace_properties(double diffuse, double reflect, double refract) {
    pt.diffuse = diffuse;
    pt.reflect = reflect;
    pt.refract = refract;

    // ensure sum to 1
    double sum = pt.diffuse + pt.reflect + pt.refract;
    pt.diffuse /= sum;
    pt.reflect /= sum;
    pt.refract /= sum;

    // create a similar set of parameters for ray tracing
    rt.ambient = 0.1 * diffuse;
    rt.diffuse = 0.9 * diffuse;
    rt.specular = reflect;
    rt.global = reflect;
    rt.alpha = 1 - refract;
    rt.shinyness = reflect * 10;
}

Plane::Plane(double r = 1, double g = 1, double b = 1) : Object(r, g, b){
    frontAndBack = 1;
}

void Plane::intersect(struct ray *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical plane.

    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);
    *lambda = -1;
    struct point norm;
    // normal of canonical plane
    norm.x = 0;
    norm.y = 0;
    norm.z = 1;
    struct point p1;

    double l;
    double d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0) {
        p1.x = -ray_transformed.p0.x;
        p1.y = -ray_transformed.p0.y;
        p1.z = -ray_transformed.p0.z;

        l = dot(&(p1), &norm) / d_dot_n;
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, p);
        double x = p->x;
        double y = p->y;
        if (fabs(p->z) < THR && fabs(p->x) <= 1 + THR && fabs(p->y) <= 1 + THR) {
            *lambda = l;
            rayPosition(ray, l, p);
            //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
            normalTransform(&norm, n, this);
            //printf("nt: (%f, %f, %f)\n", n->px, n->py, n->pz);
        }

        *a = (x + 1) / 2;
        *b = (-y + 1) / 2;
    }
}

void Plane::surfaceCoordinates(double a, double b, double *x, double *y, double *z){ 
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
    // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
    // bottom of the plane.

    struct point p;
    p.x = -2 * a + 1;
    p.y = -2 * b + 1;
    p.z = 0;
    p.w = 1;

    p = T * p;
    //matVecMult(plane->T, &p);

    *x = p.x;
    *y = p.y;
    *z = p.z;
}

void Plane::randomPoint(double *x, double *y, double *z){
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
    // Sapling should be uniform, meaning there should be an equal change of gedtting
    // any spot on the plane

    double a = drand48();
    double b = drand48();
    surfaceCoordinates(a, b, x, y, z);
}

void Sphere::intersect(struct ray *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical sphere.

    struct ray ray_transformed;
    double l1 = -1, l2 = -1;
    rayTransform(ray, &ray_transformed, this);
    solveQuadratic(&ray_transformed, &l1, &l2);

    *lambda = -1;
    if (MAX(l1, l2) >= THR) {
        if (MIN(l1, l2) <= THR) {
            *lambda = MAX(l1, l2);
        } else {
            *lambda = MIN(l1, l2);
        }
    }

    double x, y, z;
    if (*lambda != -1) {
        rayPosition(&ray_transformed, *lambda, p);
        n->x = p->x;
        n->y = p->y;
        n->z = p->z;
        x = p->x;
        y = p->y;
        z = p->z;
        n->w = 1;

        normalize(n);

        normalTransform(n, n, this);

        rayPosition(ray, *lambda, p);
        *a = 0.5 + (atan2(z, x)) / (2 * PI);
        *b = 0.5 - (asin(y)) / (PI);
    }
}

void Sphere::surfaceCoordinates(double a, double b, double *x, double *y, double *z){
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
    // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
    // bottom of the plane.

    struct point p;
    p.x = cos(a) * sin(b);
    p.y = sin(a) * sin(b);
    p.z = cos(b);
    p.w = 1;

    //matVecMult(sphere->T, &p);
    p = T * p;

    *x = p.x;
    *y = p.y;
    *z = p.z;
}

void Sphere::randomPoint(double *x, double *y, double *z){
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
    // Sapling should be uniform, meaning there should be an equal change of gedtting
    // any spot on the plane

    double a = drand48() * 2 * PI;
    double b = drand48() * 2 - 1;
    b = acos(b);
    surfaceCoordinates(a, b, x, y, z);
}

void insertObject(Object *o, Object **list) {
    if (o == NULL)
        return;
    // Inserts an object into the object list.
    if (*(list) == NULL) {
        *(list) = o;
        (*(list))->next = NULL;
    } else {
        o->next = (*(list))->next;
        (*(list))->next = o;
    }
}

inline void normalTransform(struct point *n, struct point *n_transformed, Object *obj) {
    // Computes the normal at an affinely transformed point given the original normal and the
    // object's inverse transformation. From the notes:
    // n_transformed=A^-T*n normalized.
    double x = obj->Tinv.T[0][0] * n->x + obj->Tinv.T[1][0] * n->y + obj->Tinv.T[2][0] * n->z;
    double y = obj->Tinv.T[0][1] * n->x + obj->Tinv.T[1][1] * n->y + obj->Tinv.T[2][1] * n->z;
    double z = obj->Tinv.T[0][2] * n->x + obj->Tinv.T[1][2] * n->y + obj->Tinv.T[2][2] * n->z;
    n_transformed->x = x;
    n_transformed->y = y;
    n_transformed->z = z;
    n_transformed->w = 1;
    normalize(n_transformed);
}

struct pointLS *newPLS(struct point *p0, double r, double g, double b) {
    // Allocate a new point light sourse structure. Initialize the light
    // source to the specified RGB colour
    // Note that this is a point light source in that it is a single point
    // in space, if you also want a uniform direction for light over the
    // scene (a so-called directional light) you need to place the
    // light source really far away.

    struct pointLS *ls = (struct pointLS *)calloc(1, sizeof(struct pointLS));
    if (!ls)
        fprintf(stderr, "Out of memory allocating light source!\n");
    else {
        memcpy(&ls->p0, p0, sizeof(struct point));  // Copy light source location

        ls->col.R = r;  // Store light source colour and
        ls->col.G = g;  // intensity
        ls->col.B = b;
    }
    return (ls);
}

void insertPLS(struct pointLS *l, struct pointLS **list) {
    if (l == NULL)
        return;
    // Inserts a light source into the list of light sources
    if (*(list) == NULL) {
        *(list) = l;
        (*(list))->next = NULL;
    } else {
        l->next = (*(list))->next;
        (*(list))->next = l;
    }
}
