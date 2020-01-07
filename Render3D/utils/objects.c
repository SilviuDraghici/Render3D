#include <stdio.h>
#include <stdlib.h>

#include "objects.h"

#include "mappings.h"
#include "ray.h"
#include "utils.h"

struct object *object_list;
struct pointLS *light_list;

struct object *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
    // Intialize a new plane with the specified parameters:
    // ra, rd, rs, rg - Albedos for the components of the Phong model
    // r, g, b, - Colour for this plane
    // alpha - Transparency, must be set to 1 unless you are doing refraction
    // r_index - Refraction index if you are doing refraction.
    // shiny - Exponent for the specular component of the Phong model
    //
    // The plane is defined by the following vertices (CCW)
    // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
    // With normal vector (0,0,1) (i.e. parallel to the XY plane)

    struct object *plane = (struct object *)calloc(1, sizeof(struct object));

    if (!plane)
        fprintf(stderr, "Unable to allocate new plane, out of memory!\n");
    else {
        plane->alb.ra = ra;
        plane->alb.rd = rd;
        plane->alb.rs = rs;
        plane->alb.rg = rg;
        plane->col.R = r;
        plane->col.G = g;
        plane->col.B = b;
        plane->alpha = alpha;
        plane->r_index = r_index;
        plane->shinyness = shiny;
        plane->intersect = &planeIntersect;
        plane->surfaceCoords = &planeCoordinates;
        plane->randomPoint = &planeSample;
        plane->texImg = NULL;
        plane->photonMap = NULL;
        plane->normalMap = NULL;
        plane->frontAndBack = 1;
        plane->photonMapped = 0;
        plane->normalMapped = 0;
        plane->isCSG = 0;
        plane->isLightSource = 0;
        plane->next = NULL;
    }
    return (plane);
}

void planeIntersect(struct object *plane, struct ray *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical plane.

    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, plane);
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
            normalTransform(&norm, n, plane);
            //printf("nt: (%f, %f, %f)\n", n->px, n->py, n->pz);
        }

        *a = (x + 1) / 2;
        *b = (-y + 1) / 2;
    }
}

void planeCoordinates(struct object *plane, double a, double b, double *x, double *y, double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
    // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
    // bottom of the plane.

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void planeSample(struct object *plane, double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
    // Sapling should be uniform, meaning there should be an equal change of gedtting
    // any spot on the plane

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

struct object *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
    // Intialize a new sphere with the specified parameters:
    // ra, rd, rs, rg - Albedos for the components of the Phong model
    // r, g, b, - Colour for this plane
    // alpha - Transparency, must be set to 1 unless you are doing refraction
    // r_index - Refraction index if you are doing refraction.
    // shiny -Exponent for the specular component of the Phong model
    //
    // This is assumed to represent a unit sphere centered at the origin.
    //

    struct object *sphere = (struct object *)calloc(1, sizeof(struct object));

    if (!sphere)
        fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
    else {
        sphere->alb.ra = ra;
        sphere->alb.rd = rd;
        sphere->alb.rs = rs;
        sphere->alb.rg = rg;
        sphere->col.R = r;
        sphere->col.G = g;
        sphere->col.B = b;
        sphere->alpha = alpha;
        sphere->r_index = r_index;
        sphere->shinyness = shiny;
        sphere->intersect = &sphereIntersect;
        sphere->surfaceCoords = &sphereCoordinates;
        sphere->randomPoint = &sphereSample;
        sphere->texImg = NULL;
        sphere->photonMap = NULL;
        sphere->normalMap = NULL;
        sphere->frontAndBack = 0;
        sphere->photonMapped = 0;
        sphere->normalMapped = 0;
        sphere->isCSG = 0;
        sphere->isLightSource = 0;
        sphere->next = NULL;
    }
    return (sphere);
}

void sphereIntersect(struct object *sphere, struct ray *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical sphere.

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
    struct ray ray_transformed;
    double l1 = -1, l2 = -1;
    rayTransform(ray, &ray_transformed, sphere);
    solveQuadratic(&ray_transformed, &l1, &l2);

    *lambda = -1;
    if (MAX(l1, l2) >= THR) {
        if (MIN(l1, l2) <= THR) {
            *lambda = MAX(l1, l2);
        } else {
            *lambda = MIN(l1, l2);
        }
    }
#ifdef DEBUG
    if (0) {
        printf("l1: %f, l2: %f\n", l1, l2);

        printf("sphere tray d: %f %f %f\n", ray_transformed.d.px, ray_transformed.d.py, ray_transformed.d.pz);
        printf("sphere tray p0: %f %f %f\n", ray_transformed.p0.px, ray_transformed.p0.py, ray_transformed.p0.pz);

        printf("sphere ray d: %f %f %f\n", ray->d.px, ray->d.py, ray->d.pz);
        printf("sphere ray p0: %f %f %f\n", ray->p0.px, ray->p0.py, ray->p0.pz);
    }
#endif  // DEBUG

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

        normalTransform(n, n, sphere);

        rayPosition(ray, *lambda, p);
    }
    *a = 0.5 + (atan2(z, x)) / (2 * PI);
    *b = 0.5 - (asin(y)) / (PI);
}

void sphereCoordinates(struct object *plane, double a, double b, double *x, double *y, double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
    // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
    // bottom of the plane.

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void sphereSample(struct object *plane, double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
    // Sapling should be uniform, meaning there should be an equal change of gedtting
    // any spot on the plane

    /////////////////////////////////
    // TO DO: Complete this function.
    /////////////////////////////////
}

void insertObject(struct object *o, struct object **list) {
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

inline void normalTransform(struct point *n, struct point *n_transformed, struct object *obj) {
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
