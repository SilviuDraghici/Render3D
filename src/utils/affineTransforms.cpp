#include "affineTransforms.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

inline double &d4x4(double *T, int r, int c) {
    return *(T + 4 * r + c);
}

inline double determinant_3x3(double A[3][3]) {
    //        ┌─     ─┐
    //        │ a b c │
    //  A  =  │ d e f │
    //        │ g h i │
    //        └─     ─┘
    // |A| = a(ei − fh) − b(di − fg) + c(dh − eg)

    double a = A[0][0], b = A[0][1], c = A[0][2];
    double d = A[1][0], e = A[1][1], f = A[1][2];
    double g = A[2][0], h = A[2][1], i = A[2][2];

    double eifh = (e * i) - (f * h);
    double difg = (d * i) - (f * g);
    double dheg = (d * h) - (e * g);

    return (a * eifh) - (b * difg) + (c * dheg);
}

matrix Sc(double Xscale, double Yscale, double Zscale) {
    // Returns tranform for left multiplying
    // to a transform or point to
    // scale the object
    matrix scale;
    scale.T[0][0] = Xscale;
    scale.T[1][1] = Yscale;
    scale.T[2][2] = Zscale;
    return scale;
}

matrix Sc(double uniform_scale) {
    // Returns tranform for left multiplying
    // to a transform or point to
    // scale the object
    matrix scale;
    scale.T[0][0] = uniform_scale;
    scale.T[1][1] = uniform_scale;
    scale.T[2][2] = uniform_scale;
    return scale;
}

matrix Tr(double Xtranslate, double Ytranslate, double Ztranslate) {
    // Returns tranform for left multiplying
    // to a transform or point to
    // tranlate the object
    matrix translate;
    translate.T[0][3] = Xtranslate;
    translate.T[1][3] = Ytranslate;
    translate.T[2][3] = Ztranslate;
    return translate;
}
matrix RotX(double theta) {
    // Returns tranform for left multiplying
    // to a transform or point to
    // rotate theta *RADIANS* around the
    // X axis.
    matrix rotateX;
    rotateX.T[1][1] = cos(theta);
    rotateX.T[1][2] = -sin(theta);
    rotateX.T[2][1] = sin(theta);
    rotateX.T[2][2] = cos(theta);
    return rotateX;
}
matrix RotY(double theta) {
    // Returns tranform for left multiplying
    // to a transform or point to
    // rotate theta *RADIANS* around the
    // Y axis.
    matrix rotateY;
    rotateY.T[0][0] = cos(theta);
    rotateY.T[0][2] = sin(theta);
    rotateY.T[2][0] = -sin(theta);
    rotateY.T[2][2] = cos(theta);
    return rotateY;
}
matrix RotZ(double theta) {
    // Returns tranform for left multiplying
    // to a transform or point to
    // rotate theta *RADIANS* around the
    // Z axis.
    matrix rotateZ;
    rotateZ.T[0][0] = cos(theta);
    rotateZ.T[0][1] = -sin(theta);
    rotateZ.T[1][0] = sin(theta);
    rotateZ.T[1][1] = cos(theta);
    return rotateZ;
}

double eye4x4[4][4] = {{1.0, 0.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0, 0.0},
                       {0.0, 0.0, 1.0, 0.0},
                       {0.0, 0.0, 0.0, 1.0}};

void invert(double *T, double *Tinv) {
    // Computes the inverse of transformation matrix T.
    // the result is returned in Tinv.

    // because all the affine transforms have the form:
    //        ┌─       ─┐
    //        │ a b c d │
    //  T  =  │ e f g h │
    //        │ i j k l │
    //        │ 0 0 0 1 │
    //        └─       ─┘
    // the determinant of the top left 3x3 is equal to T's determinant
    double A[3][3] = {{d4x4(T, 0, 0), d4x4(T, 0, 1), d4x4(T, 0, 2)},
                      {d4x4(T, 1, 0), d4x4(T, 1, 1), d4x4(T, 1, 2)},
                      {d4x4(T, 2, 0), d4x4(T, 2, 1), d4x4(T, 2, 2)}};

    double determinant = determinant_3x3(A);

    if (determinant == 0) {
        fprintf(stderr, "Error: Matrix not invertible for this object, returning identity\n");
        memcpy(Tinv, eye4x4, 16 * sizeof(double));
        return;
    }

    //builds the sub-matrices to build the cofactor matrix 
    int Tr, Tc;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            Tr = Tc = 0;
            for (int i = 0; i < 3; i++) {
                if (i == r) { Tr++; }
                
                Tc = 0;
                for (int j = 0; j < 3; j++) {
                    if (j == c) { Tc++; }
                    A[i][j] = d4x4(T, Tr, Tc);
                    Tc++;
                }

                Tr++;
            }
            //build the cofactor matrix
            d4x4(Tinv, r, c) = pow(-1, (r + 1) + (c + 1)) * determinant_3x3(A);
        }
    }

    double inv_det = 1 / determinant;

    //transpose(Tinv) = 1/det(T) * cofactor(T)
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            d4x4(Tinv, r, c) *= inv_det;
        }
    }
    
    //Takes the transpose of transpose(Tinv)
    double temp;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < r; c++) {
            temp = d4x4(Tinv, r, c);
            d4x4(Tinv, r, c) = d4x4(Tinv, c, r);
            d4x4(Tinv, c, r) = temp;
        }
    }
}

void printmatrix(matrix &matrix) {
    fprintf(stderr, "Matrix contains:\n");
    fprintf(stderr, "%f %f %f %f\n", matrix.T[0][0], matrix.T[0][1], matrix.T[0][2], matrix.T[0][3]);
    fprintf(stderr, "%f %f %f %f\n", matrix.T[1][0], matrix.T[1][1], matrix.T[1][2], matrix.T[1][3]);
    fprintf(stderr, "%f %f %f %f\n", matrix.T[2][0], matrix.T[2][1], matrix.T[2][2], matrix.T[2][3]);
    fprintf(stderr, "%f %f %f %f\n", matrix.T[3][0], matrix.T[3][1], matrix.T[3][2], matrix.T[3][3]);
}

//////////////////////////////////
// Importance sampling for BRDF
//////////////////////////////////
void cosWeightedSample(struct point *n, struct point *d) {
    // This function returns a randomly sampled direction over
    // a hemisphere whose pole is the normal direction n. The
    // sampled direction comes from a distribution weighted
    // by the cosine of the angle between n and d.
    // Use this for importance sampling for diffuse surfaces.

    double u1, r, theta, phi;
    double x, y, z, c;
    matrix R = I();
    struct point nz;

    // Random sample on hemisphere with cosine-weighted distribution
    u1 = drand48();
    r = sqrt(u1);
    theta = 2 * PI * drand48();
    nz.x = r * cos(theta);
    nz.y = r * sin(theta);
    nz.z = sqrt(1.0 - (nz.x * nz.x) - (nz.y * nz.y));
    nz.w = 1;

    // Rotation based on cylindrical coordinate conversion
    phi = acos(n->z);
    theta = atan2(n->y, n->x);

    //RotateYMat(R, phi);
    //RotateZMat(R, theta);
    R *= RotY(phi);
    R *= RotZ(theta);

    // Rotate d to align with normal
    *d = R * nz;

    return;
}