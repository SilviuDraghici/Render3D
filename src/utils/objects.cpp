#include "objects.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <string>

#include "affineTransforms.h"
#include "mappings.h"
#include "ray.h"
#include "utils.h"

void union_bounds(Bounds &a, point &b, Bounds &union_box) {
    union_box.min.x = MIN(a.min.x, b.x);
    union_box.min.y = MIN(a.min.y, b.y);
    union_box.min.z = MIN(a.min.z, b.z);
    union_box.max.x = MAX(a.max.x, b.x);
    union_box.max.y = MAX(a.max.y, b.y);
    union_box.max.z = MAX(a.max.z, b.z);
}

void union_bounds(Bounds &a, Bounds &b, Bounds &union_box) {
    //assigns the bounds of union_box to fit a and b
    union_box.min.x = MIN(a.min.x, b.min.x);
    union_box.min.y = MIN(a.min.y, b.min.y);
    union_box.min.z = MIN(a.min.z, b.min.z);
    union_box.max.x = MAX(a.max.x, b.max.x);
    union_box.max.y = MAX(a.max.y, b.max.y);
    union_box.max.z = MAX(a.max.z, b.max.z);
}

Axis Primitive::longestAxis() const {
    double x = max_x() - min_x();
    double y = max_y() - min_y();
    double z = max_z() - min_z();
    if (x > y && x > z) {
        return Axis::X;
    } else if (y > z) {
        return Axis::Y;
    } else {
        return Axis::Z;
    }
}

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
    normalMap = NULL;
    alphaMap = NULL;
    frontAndBack = 0;
    isLightSource = 0;
    T = I();
    next = NULL;
}

void Object::set_color(double r, double g, double b) {
    col.R = r;
    col.G = g;
    col.B = b;
}

void Object::set_rayTrace_properties(double ambient, double diffuse,
                                     double specular, double global,
                                     double alpha, double shiny) {
    rt.ambient = ambient;
    rt.diffuse = diffuse;
    rt.specular = specular;
    rt.global = global;
    rt.alpha = alpha;
    rt.shinyness = shiny;

    // create a similar set of params incase the object is drawn in pathtrace
    // mode
    pt.diffuse = alpha * (ambient + diffuse);
    pt.reflect = alpha * (1 - diffuse);
    pt.refract = 1 - alpha;

    // ensure sum to 1
    double sum = pt.diffuse + pt.reflect + pt.refract;
    pt.diffuse /= sum;
    pt.reflect /= sum;
    pt.refract /= sum;
}

void Object::set_pathTrace_properties(double diffuse, double reflect,
                                      double refract) {
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

void Object::surfaceCoordinates(double a, double b, double *x, double *y,
                                double *z) {
    fprintf(stderr, "Object::surfaceCoordinates\n");
}

void Object::set_canonical_bounds() {
    std::cout << "Canonical Bounds not set!\n";
}

void Object::randomPoint(double *x, double *y, double *z) {
    point r = 0;
    r = T * r;
    *x = r.x;
    *y = r.y;
    *z = r.z;
    fprintf(stderr, "Object::randomPoint\n");
}

double Object::intersect(struct ray *r, double lambda) {
    lambda = INFINITY;
    double dummy_doub;
    point dummy_point;
    intersect(r, &lambda, &dummy_point, &dummy_point, &dummy_doub, &dummy_doub);
    return lambda;
}

void Object::invert_and_bound() {
    // calculates inverse matrix and sets w_bound
    invert(&T.T[0][0], &Tinv.T[0][0]);

    set_canonical_bounds();
    //std::cout << label << "----------------------\n";
    //std::cout << "c bound:"<< w_bound.min << " " << w_bound.max << "\n";

    //convert canonical bounds to world coordinates
    point a, b, c, d, e, f, g, h;
    a = T * point(w_bound.min.x, w_bound.min.y , w_bound.min.z);
    b = T * point(w_bound.min.x, w_bound.min.y , w_bound.max.z);
    c = T * point(w_bound.min.x, w_bound.max.y , w_bound.min.z);
    d = T * point(w_bound.min.x, w_bound.max.y , w_bound.max.z);
    
    e = T * point(w_bound.max.x, w_bound.max.y , w_bound.max.z);
    f = T * point(w_bound.max.x, w_bound.max.y , w_bound.min.z);
    g = T * point(w_bound.max.x, w_bound.min.y , w_bound.max.z);
    h = T * point(w_bound.max.x, w_bound.min.y , w_bound.min.z);

    //printmatrix(T);
    //std::cout << "\n" << T << "\n";
    //std::cout << "w bound:"<< a << " " << b << "\n";

    // find world_bounds
    w_bound.min.x = MIN(MIN(MIN(MIN(a.x, b.x), c.x), d.x), MIN(MIN(MIN(e.x, f.x), g.x), h.x));
    w_bound.min.y = MIN(MIN(MIN(MIN(a.y, b.y), c.y), d.y), MIN(MIN(MIN(e.y, f.y), g.y), h.y));
    w_bound.min.z = MIN(MIN(MIN(MIN(a.z, b.z), c.z), d.z), MIN(MIN(MIN(e.z, f.z), g.z), h.z));
    w_bound.max.x = MAX(MAX(MAX(MAX(a.x, b.x), c.x), d.x), MAX(MAX(MAX(e.x, f.x), g.x), h.x));
    w_bound.max.y = MAX(MAX(MAX(MAX(a.y, b.y), c.y), d.y), MAX(MAX(MAX(e.y, f.y), g.y), h.y));
    w_bound.max.z = MAX(MAX(MAX(MAX(a.z, b.z), c.z), d.z), MAX(MAX(MAX(e.z, f.z), g.z), h.z));

}

bool Object::isprim() {
    return true;
}

double Object::min_x() const { return w_bound.min.x; }
double Object::min_y() const { return w_bound.min.y; }
double Object::min_z() const { return w_bound.min.z; }
double Object::max_x() const { return w_bound.max.x; }
double Object::max_y() const { return w_bound.max.y; }
double Object::max_z() const { return w_bound.max.z; }

///////////////////////////////Plane//////////////////////////////////////

Plane::Plane(double r = 1, double g = 1, double b = 1) : Object(r, g, b) {
    frontAndBack = 1;
}

void Plane::intersect(struct ray *ray, double *lambda, struct point *p,
                      struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical plane.

    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);
    *lambda = INFINITY;
    struct point norm;
    // normal of canonical plane
    norm.x = 0;
    norm.y = 0;
    norm.z = 1;
    struct point p1;

    double l;
    double d_dot_n = dot(&(ray_transformed.d), &norm);
    if (d_dot_n != 0) {
        p1 = -ray_transformed.p0;

        l = dot(&(p1), &norm) / d_dot_n;
        if (l < THR) return;
        // Check if the intersection point is inside the plane
        rayPosition(&ray_transformed, l, p);
        double x = p->x;
        double y = p->y;
        if (fabs(p->z) < THR && fabs(p->x) <= 1 + THR &&
            fabs(p->y) <= 1 + THR) {
            *lambda = l;
            rayPosition(ray, l, p);
            // printf("n: (%f, %f, %f)\n", n->x, n->y, n->z);
            normalTransform(&norm, n, this);
            // printf("nt: (%f, %f, %f)\n", n->x, n->y, n->z);
        }

        *a = (x + 1) / 2;
        *b = (-y + 1) / 2;
    }
}

void Plane::set_canonical_bounds() {
    w_bound.min = point(-1, -1, 0);
    w_bound.max = point(1, 1, 0);
}

void Plane::surfaceCoordinates(double a, double b, double *x, double *y,
                               double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b in [0,1]. 'a' controls displacement from the left side of
    // the plane, 'b' controls displacement from the bottom of the plane.

    struct point p;
    p.x = -2 * a + 1;
    p.y = -2 * b + 1;
    p.z = 0;
    p.w = 1;

    p = T * p;
    // matVecMult(plane->T, &p);

    *x = p.x;
    *y = p.y;
    *z = p.z;
}

void Plane::randomPoint(double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the
    // plane Sapling should be uniform, meaning there should be an equal chance
    // of getting any spot on the plane

    double a = drand48();
    double b = drand48();
    surfaceCoordinates(a, b, x, y, z);
}

void Sphere::intersect(struct ray *ray, double *lambda, struct point *p,
                       struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical sphere.

    struct ray ray_transformed;
    double l1 = -1, l2 = -1;
    rayTransform(ray, &ray_transformed, this);
    solveQuadratic(&ray_transformed, &l1, &l2);

    *lambda = INFINITY;
    if (MAX(l1, l2) >= THR) {
        if (MIN(l1, l2) <= THR) {
            *lambda = MAX(l1, l2);
        } else {
            *lambda = MIN(l1, l2);
        }
    }

    double x, y, z;
    if (*lambda != INFINITY) {
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

        // use the lambda to set the position along the untransformed ray
        rayPosition(ray, *lambda, p);
        *a = 0.5 + (atan2(z, x)) / (2 * PI);
        *b = 0.5 - (asin(y)) / (PI);
    }
}

void Sphere::set_canonical_bounds() {
    w_bound.min = point(-1, -1, -1);
    w_bound.max = point(1, 1, 1);
}

void Sphere::surfaceCoordinates(double a, double b, double *x, double *y,
                                double *z) {
    // Return in (x,y,z) the coordinates of a point on the plane given by the 2
    // parameters a,b in [0,1]. 'a' controls displacement from the left side of
    // the plane, 'b' controls displacement from the bottom of the plane.

    struct point p;
    p.x = cos(a) * sin(b);
    p.y = sin(a) * sin(b);
    p.z = cos(b);
    p.w = 1;

    // matVecMult(sphere->T, &p);
    p = T * p;

    *x = p.x;
    *y = p.y;
    *z = p.z;
}

void Sphere::randomPoint(double *x, double *y, double *z) {
    // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the
    // plane Sapling should be uniform, meaning there should be an equal change
    // of gedtting any spot on the plane

    double a = drand48() * 2 * PI;
    double b = drand48() * 2 - 1;
    b = acos(b);
    surfaceCoordinates(a, b, x, y, z);
}

void Box::intersect(struct ray *ray, double *lambda, struct point *p,
                    struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified canonical box.

    *lambda = INFINITY;

    // current intersection lambda
    double b_lambda;

    // use the objects inverse transform to get a ray in the cannonical space
    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    // y-z plane box face at x = -0.5
    b_lambda = (-0.5 - ray_transformed.p0.x) / ray_transformed.d.x;
    if (THR < b_lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->y && p->y < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at x = -0.5
            *lambda = b_lambda;
            *n = point(-1, 0, 0);
        }
    }

    // y-z plane box face at x = 0.5
    b_lambda = (0.5 - ray_transformed.p0.x) / ray_transformed.d.x;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->y && p->y < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at x = -0.5
            *lambda = b_lambda;
            *n = point(1, 0, 0);
        }
    }

    // x-z plane box face at y = -0.5
    b_lambda = (-0.5 - ray_transformed.p0.y) / ray_transformed.d.y;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at y = -0.5
            *lambda = b_lambda;
            *n = point(0, -1, 0);
        }
    }

    // x-z plane box face at y = 0.5
    b_lambda = (0.5 - ray_transformed.p0.y) / ray_transformed.d.y;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->z && p->z < 0.5)) {
            // hit the face at y = -0.5
            *lambda = b_lambda;
            *n = point(0, 1, 0);
        }
    }

    // x-y plane box face at z = -0.5
    b_lambda = (-0.5 - ray_transformed.p0.z) / ray_transformed.d.z;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->y && p->y < 0.5)) {
            // hit the face at z = -0.5
            *lambda = b_lambda;
            *n = point(0, 0, -1);
        }
    }

    // x-y plane box face at z = 0.5
    b_lambda = (0.5 - ray_transformed.p0.z) / ray_transformed.d.z;
    if (THR < b_lambda && b_lambda < *lambda) {
        rayPosition(&ray_transformed, b_lambda, p);
        if ((-0.5 < p->x && p->x < 0.5) && (-0.5 < p->y && p->y < 0.5)) {
            // hit the face at z = -0.5
            *lambda = b_lambda;
            *n = point(0, 0, 1);
        }
    }

    // printf("intersect point: (%f, %f, %f)\n", p->x, p->y, p->z);

    // if there is an intersection, update the normal and intersection point
    if (*lambda != INFINITY) {
        normalTransform(n, n, this);
        rayPosition(ray, *lambda, p);
    } else {  // return no intersection flag
        *lambda = -1;
    }
}

Triangle::Triangle(double r = 1, double g = 1, double b = 1) : Object(r, g, b) {
    frontAndBack = 1;
    // equilateral triangle on x-y plane unit circle centered at 0
    p1.x = -0.866, p1.y = -0.5;
    p2.x = 0.866, p2.y = -0.5;
    p3.x = 0.0, p3.y = 1.0;
}

void Triangle::intersect(struct ray *ray, double *lambda, struct point *p,
                         struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified triangle.

    // default no intersection
    *lambda = -1;

    // use the objects inverse transform to get a ray in the cannonical space
    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    point e12 = p2 - p1;
    point e23 = p3 - p2;

    point normal = cross(&e12, &e23);
    // std::cout << "normal: " << normal << std::endl;

    // ray plane intersection calculation
    point r = p1 - ray_transformed.p0;
    double d_dot_n = dot(&ray_transformed.d, &normal);
    double r_dot_n = dot(&r, &normal);
    double t_lambda = r_dot_n / d_dot_n;

    if (THR < t_lambda) {  // check that intersection with plane is positive
        rayPosition(&ray_transformed, t_lambda, p);
        // e12 has already been calculated
        point v1i = *p - p1;

        point crossp = cross(&e12, &v1i);
        if (dot(&crossp, &normal) >= 0) {  // check within first edge
            // e23 has already been calculated
            point v2i = *p - p2;
            crossp = cross(&e23, &v2i);
            if (dot(&crossp, &normal) >= 0) {  // check within second edge
                point e31 = p1 - p3;
                point v3i = *p - p3;
                crossp = cross(&e31, &v3i);
                if (dot(&crossp, &normal) >= 0) {  // check within third edge
                    // intersection is within triangle
                    *lambda = t_lambda;
                    normalTransform(&normal, n, this);
                    normalize(n);
                    rayPosition(ray, *lambda, p);
                }
            }
        }
    }
}

void Triangle::setPoints(double x1, double y1, double z1, double x2, double y2,
                         double z2, double x3, double y3, double z3) {
    p1.x = x1, p1.y = y1, p1.z = z1;
    p2.x = x2, p2.y = y2, p2.z = z2;
    p3.x = x3, p3.y = y3, p3.z = z3;
}

void Triangle::setPoints(point p1, point p2, point p3) {
    this->p1 = p1, this->p2 = p2, this->p3 = p3;
}

Polygon::Polygon(double r = 1, double g = 1, double b = 1) : Object(r, g, b) {
    frontAndBack = 1;
    p = (point *)malloc(3 * sizeof(point));
    e = (point *)malloc(3 * sizeof(point));
    // equilateral triangle on x-y plane unit circle centered at 0
    p[0].x = -0.866, p[0].y = -0.5;
    p[1].x = 0.866, p[1].y = -0.5;
    p[2].x = 0.0, p[2].y = 1.0;
    calculate_edge_vectors();
}

void Polygon::intersect(struct ray *ray, double *lambda, struct point *pi,
                        struct point *n, double *a, double *b) {
    // Computes and returns the value of 'lambda' at the intersection
    // between the specified ray and the specified polygon.

    // default no intersection
    *lambda = -1;

    // use the objects inverse transform to get a ray in the cannonical space
    struct ray ray_transformed;
    rayTransform(ray, &ray_transformed, this);

    point e12 = p[1] - p[0];
    point e23 = p[2] - p[1];

    point normal = cross(&e[0], &e[1]);
    // std::cout << "normal: " << normal << std::endl;

    // ray plane intersection calculation
    point r = p[0] - ray_transformed.p0;
    double d_dot_n = dot(&ray_transformed.d, &normal);
    double r_dot_n = dot(&r, &normal);
    double p_lambda = r_dot_n / d_dot_n;

    if (THR < p_lambda) {  // check that intersection with plane is positive
        rayPosition(&ray_transformed, p_lambda, pi);
        point v, crossp;
        int in_edges = 0;
        for (int i = 0; i < numPoints; i++) {
            v = *pi - p[i];
            crossp = cross(&e[i], &v);
            if (dot(&crossp, &normal) < 0) break;
            in_edges++;
        }
        if (in_edges == numPoints) {
            // intersection is within polygon
            *lambda = p_lambda;
            normalTransform(&normal, n, this);
            normalize(n);
            rayPosition(ray, *lambda, pi);
        }
    }
}

void Polygon::setNumPoints(int num) {
    free(p);
    free(e);
    numPoints = num;
    currPoint = 0;
    p = (point *)malloc(numPoints * sizeof(point));
    e = (point *)malloc(numPoints * sizeof(point));
}

void Polygon::addPoint(point point) {
    p[currPoint] = point;
    currPoint++;
    if (currPoint == numPoints) {
        calculate_edge_vectors();
    }
}

void Polygon::addPoint(double x, double y, double z) {
    p[currPoint] = point(x, y, z);
    currPoint++;
    if (currPoint == numPoints) {
        calculate_edge_vectors();
    }
}

void Polygon::calculate_edge_vectors() {
    // calculate vectors for all edges
    for (int i = 0; i < numPoints - 1; i++) {
        e[i] = p[i + 1] - p[i];
    }
    // last edge goes from last point to first point
    e[numPoints - 1] = p[0] - p[numPoints - 1];
}

void normalTransform(struct point *n, struct point *n_transformed,
                     Object *obj) {
    // Computes the normal at an affinely transformed point given the original
    // normal and the object's inverse transformation. From the notes:
    // n_transformed=A^-T*n normalized.
    double x = obj->Tinv.T[0][0] * n->x + obj->Tinv.T[1][0] * n->y +
               obj->Tinv.T[2][0] * n->z;
    double y = obj->Tinv.T[0][1] * n->x + obj->Tinv.T[1][1] * n->y +
               obj->Tinv.T[2][1] * n->z;
    double z = obj->Tinv.T[0][2] * n->x + obj->Tinv.T[1][2] * n->y +
               obj->Tinv.T[2][2] * n->z;
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
        memcpy(&ls->p0, p0,
               sizeof(struct point));  // Copy light source location

        ls->col.R = r;  // Store light source colour and
        ls->col.G = g;  // intensity
        ls->col.B = b;
    }
    return (ls);
}

void insertPLS(struct pointLS *l, struct pointLS **list) {
    if (l == NULL) return;
    // Inserts a light source into the list of light sources
    if (*(list) == NULL) {
        *(list) = l;
        (*(list))->next = NULL;
    } else {
        l->next = (*(list))->next;
        (*(list))->next = l;
    }
}
