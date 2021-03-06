#include "utils.h"

#include <iomanip>
#include <math.h>
#include <string>

#include "ray.h"

std::ostream &operator<<(std::ostream &strm, const struct point &a) {
    return strm << "(" << a.x << ", " << a.y << ", " << a.z << ", " << a.w << ")";
}

std::ostream &operator<<(std::ostream &strm, const matrix &m){
    int wid = 0;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            wid = MAX(wid, std::to_string(m.T[i][j]).length());
        }
    }
    strm << std::setw(wid) << std::to_string(m.T[0][0]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[0][1]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[0][2]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[0][3]) << "\n";
    strm << std::setw(wid) << std::to_string(m.T[1][0]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[1][1]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[1][2]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[1][3]) << "\n";
    strm << std::setw(wid) << std::to_string(m.T[2][0]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[2][1]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[2][2]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[2][3]) << "\n";
    strm << std::setw(wid) << std::to_string(m.T[3][0]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[3][1]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[3][2]) << " ";
    strm << std::setw(wid) << std::to_string(m.T[3][3]) << "\n";
    return strm;
}

matrix I() {
    matrix i;
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

struct point cross(struct point *u, struct point *v) {
    // Allocates and returns a vector with the cross product u x v.
    // The function assumes the w components of both vectors
    // are 1.
    struct point cp;

    cp.x = (u->y * v->z) - (v->y * u->z);
    cp.y = (v->x * u->z) - (u->x * v->z);
    cp.z = (u->x * v->y) - (v->x * u->y);
    cp.w = 1;
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

void solveQuadratic(struct Ray *ray, double *l1, double *l2) {
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
    double a = xor128() * 2 * PI;
    double b = xor128() * 2 - 1;
    b = acos(b);

    d->x = cos(a) * sin(b);
    d->y = sin(a) * sin(b);
    d->z = cos(b);
    d->w = 1;
    while (dot(n, d) < 0) {
        a = xor128() * 2 * PI;
        b = xor128() * 2 - 1;
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
        U1 = -1 + xor128() * 2;
        U2 = -1 + xor128() * 2;
        W = pow(U1, 2) + pow(U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt((-2 * log(W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * X1);
}