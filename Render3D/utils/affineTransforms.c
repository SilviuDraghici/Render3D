#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "affineTransforms.h"

#include "svdDynamic.c"
#include "svdDynamic.h"

struct matrix I() {
   struct matrix i;
   return i;
}

struct matrix Sc(double Xscale, double Yscale, double Zscale) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // scale the object
   struct matrix scale;
   scale.T[0][0] = Xscale;
   scale.T[1][1] = Yscale;
   scale.T[2][2] = Zscale;
   return scale;
}

struct matrix Sc(double uniform_scale) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // scale the object
   struct matrix scale;
   scale.T[0][0] = uniform_scale;
   scale.T[1][1] = uniform_scale;
   scale.T[2][2] = uniform_scale;
   return scale;
}

struct matrix Tr(double Xtranslate, double Ytranslate, double Ztranslate) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // tranlate the object
   struct matrix translate;
   translate.T[0][3] = Xtranslate;
   translate.T[1][3] = Ytranslate;
   translate.T[2][3] = Ztranslate;
   return translate;
}
struct matrix RotX(double theta) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // rotate theta *RADIANS* around the
   // X axis.
   struct matrix rotateX;
   rotateX.T[1][1] = cos(theta);
   rotateX.T[1][2] = -sin(theta);
   rotateX.T[2][1] = sin(theta);
   rotateX.T[2][2] = cos(theta);
   return rotateX;
}
struct matrix RotY(double theta) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // rotate theta *RADIANS* around the
   // Y axis.
   struct matrix rotateY;
   rotateY.T[0][0] = cos(theta);
   rotateY.T[0][2] = sin(theta);
   rotateY.T[2][0] = -sin(theta);
   rotateY.T[2][2] = cos(theta);
   return rotateY;
}
struct matrix RotZ(double theta) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // rotate theta *RADIANS* around the
   // Z axis.
   struct matrix rotateZ;
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

   double *U, *s, *V, *rv1;
   int singFlag, i;

   // Invert the affine transform
   U = NULL;
   s = NULL;
   V = NULL;
   rv1 = NULL;
   singFlag = 0;

   SVD(T, 4, 4, &U, &s, &V, &rv1);
   if (U == NULL || s == NULL || V == NULL) {
      fprintf(stderr, "Error: Matrix not invertible for this object, returning identity\n");
      memcpy(Tinv, eye4x4, 16 * sizeof(double));
      return;
   }

   // Check for singular matrices...
   for (i = 0; i < 4; i++)
      if (*(s + i) < 1e-9)
         singFlag = 1;
   if (singFlag) {
      fprintf(stderr, "Error: Transformation matrix is singular, returning identity\n");
      memcpy(Tinv, eye4x4, 16 * sizeof(double));
      return;
   }

   // Compute and store inverse matrix
   InvertMatrix(U, s, V, 4, Tinv);

   free(U);
   free(s);
   free(V);
}

void printmatrix(struct matrix matrix) {
   fprintf(stderr, "Matrix contains:\n");
   fprintf(stderr, "%f %f %f %f\n", matrix.T[0][0], matrix.T[0][1], matrix.T[0][2], matrix.T[0][3]);
   fprintf(stderr, "%f %f %f %f\n", matrix.T[1][0], matrix.T[1][1], matrix.T[1][2], matrix.T[1][3]);
   fprintf(stderr, "%f %f %f %f\n", matrix.T[2][0], matrix.T[2][1], matrix.T[2][2], matrix.T[2][3]);
   fprintf(stderr, "%f %f %f %f\n", matrix.T[3][0], matrix.T[3][1], matrix.T[3][2], matrix.T[3][3]);
}