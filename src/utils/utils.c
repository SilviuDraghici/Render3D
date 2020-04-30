#include "utils.h"

#include <math.h>

#include "ray.h"

struct matrix I() {
    struct matrix i;
    memset(&i.T[0][0], 0, 16 * sizeof(double));
    i.T[0][0] = 1;
    i.T[1][1] = 1;
    i.T[2][2] = 1;
    i.T[3][3] = 1;
    return i;
}

double dot(struct point *u, struct point *v) {
    // Computes the dot product of 3D vectors u and v.
    // The function assumes the w components of both vectors
    // are 1.
    return ((u->x * v->x) + (u->y * v->y) + (u->z * v->z));
}

struct point *cross(struct point *u, struct point *v) {
    // Allocates and returns a vector with the cross product u x v.
    // The function assumes the w components of both vectors
    // are 1.
    struct point *cp;
    cp = (struct point *)calloc(1, sizeof(struct point));

    cp->x = (u->y * v->z) - (v->y * u->z);
    cp->y = (v->x * u->z) - (u->x * v->z);
    cp->z = (u->x * v->y) - (v->x * u->y);
    cp->w = 1;
    return (cp);
}

void normalize(struct point *v) {
    // Normalizes a vector to unit length.
    // Note that this assumes the w component of v is one, so make
    // sure you don't have homogeneous vectors with w != 1
    // floating around or things will go south.
    double l;
    l = v->x * v->x;
    l += (v->y * v->y);
    l += (v->z * v->z);
    l = 1.0 / sqrt(l);
    v->x *= l;
    v->y *= l;
    v->z *= l;
}

void solveQuadratic(struct ray *ray, double *l1, double *l2) {
    struct point a_sub_c;
    double A, B, C, D;
    a_sub_c = ray->p0;
    A = dot(&(ray->d), &(ray->d));
    B = dot(&a_sub_c, &(ray->d));
    C = dot(&a_sub_c, &a_sub_c) - 1;
    D = B * B - A * C;

    if (D < 0) {
        return;
    }

    if (D > 0) {
        *l1 = -B / A - sqrt(D) / A;
        *l2 = -B / A + sqrt(D) / A;
    }
}

void hemiSphereCoordinates(struct point *n, struct point *d) {
    // returns a random vector from the hemisphere around n
    double a = drand48() * 2 * PI;
    double b = drand48() * 2 - 1;
    b = acos(b);

    d->x = cos(a) * sin(b);
    d->y = sin(a) * sin(b);
    d->z = cos(b);
    d->w = 1;
    while (dot(n, d) < 0) {
        a = drand48() * 2 * PI;
        b = drand48() * 2 - 1;
        b = acos(b);
        d->x = cos(a) * sin(b);
        d->y = sin(a) * sin(b);
        d->z = cos(b);
    }
}

double rand_normal_dist(double mu, double sigma) {
    // https://en.wikipedia.org/wiki/Marsaglia_polar_method
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1) {
        call = !call;
        return (mu + sigma * X2);
    }

    do {
        U1 = -1 + drand48() * 2;
        U2 = -1 + drand48() * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * X1);
}