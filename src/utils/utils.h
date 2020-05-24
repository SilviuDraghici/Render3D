#include <string.h>
#include <iostream>

#ifndef UTILS_H
#define UTILS_H

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define THR 0.0001
#define PI 3.14159265354

/* The structure below defines a point in 3D homogeneous coordinates        */
struct point {
    double x;
    double y;
    double z;
    double w;
    point(double px = 0, double py = 0, double pz = 0) {
        x = px;
        y = py;
        z = pz;
        w = 1;
    }
    point &operator=(const point &pt) {
        x = pt.x;
        y = pt.y;
        z = pt.z;
        return *this;
    }

    point &operator*=(double mult) {
        x *= mult;
        y *= mult;
        z *= mult;
        return *this;
    }

    point operator+(const point &a) const {
        return point(x + a.x, y + a.y, z + a.z);
    }
    point operator-(const point &a) const {
        return point(x - a.x, y - a.y, z - a.z);
    }

    point operator-() const { return point(-x, -y, -z); }

    point operator*(double scalar) const {
        return point(x * scalar, y * scalar, z * scalar);
    }

    //std::ostream &operator<<(std::ostream &strm) {
      //  return strm << "(" << x << ", " << y << ", " << z << ", " << w << ")";
    //}
};

std::ostream &operator<<(std::ostream &strm, const struct point &a);

struct matrix {
    // this struct defines a 4x4 matrix used for
    // 3d affine transforms in homogeneous coordinates
    double T[4][4];
    matrix() {
        memset(&T[0][0], 0, 16 * sizeof(double));
        T[0][0] = 1;
        T[1][1] = 1;
        T[2][2] = 1;
        T[3][3] = 1;
    }

    matrix operator*(const matrix &b) const {
        struct matrix result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result.T[i][j] = (T[i][0] * b.T[0][j]) + (T[i][1] * b.T[1][j]) +
                                 (T[i][2] * b.T[2][j]) + (T[i][3] * b.T[3][j]);

        return result;
    }
    matrix &operator*=(const matrix &b) {
        *this = b * *this;
        return *this;
    }
    point operator*(const point &p) const {
        struct point pr;
        pr.x = (T[0][0] * p.x) + (T[0][1] * p.y) + (T[0][2] * p.z) +
               (T[0][3] * p.w);
        pr.y = (T[1][0] * p.x) + (T[1][1] * p.y) + (T[1][2] * p.z) +
               (T[1][3] * p.w);
        pr.z = (T[2][0] * p.x) + (T[2][1] * p.y) + (T[2][2] * p.z) +
               (T[2][3] * p.w);
        return pr;
    }
    point operator*(const point *p) const {
        struct point pr;
        pr.x = (T[0][0] * p->x) + (T[0][1] * p->y) + (T[0][2] * p->z) +
               (T[0][3] * p->w);
        pr.y = (T[1][0] * p->x) + (T[1][1] * p->y) + (T[1][2] * p->z) +
               (T[1][3] * p->w);
        pr.z = (T[2][0] * p->x) + (T[2][1] * p->y) + (T[2][2] * p->z) +
               (T[2][3] * p->w);
        return pr;
    }
};

struct matrix I();

double dot(struct point *u, struct point *v);
struct point cross(struct point *u, struct point *v);

void normalize(struct point *v);

void solveQuadratic(struct ray *ray, double *l1, double *l2);

void hemiSphereCoordinates(struct point *n, struct point *d);

double rand_normal_dist(double mu, double sigma);
#endif