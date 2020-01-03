#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "affineTransforms.h"

struct transform I() {
   struct transform i;
   return i;
}

struct transform Sc(double Xscale, double Yscale, double Zscale) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // scale the object
   struct transform scale;
   scale.T[0][0] = Xscale;
   scale.T[1][1] = Yscale;
   scale.T[2][2] = Zscale;
   return scale;
}

struct transform Sc(double uniform_scale) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // scale the object
   struct transform scale;
   scale.T[0][0] = uniform_scale;
   scale.T[1][1] = uniform_scale;
   scale.T[2][2] = uniform_scale;
   return scale;
}

struct transform Tr(double Xtranslate, double Ytranslate, double Ztranslate) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // tranlate the object
   struct transform translate;
   translate.T[0][3] = Xtranslate;
   translate.T[1][3] = Ytranslate;
   translate.T[2][3] = Ztranslate;
   return translate;
}
struct transform RotX(double theta) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // rotate theta *RADIANS* around the
   // X axis.
   struct transform rotateX;
   rotateX.T[1][1] = cos(theta);
   rotateX.T[1][2] = -sin(theta);
   rotateX.T[2][1] = sin(theta);
   rotateX.T[2][2] = cos(theta);
   return rotateX;
}
struct transform RotY(double theta) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // rotate theta *RADIANS* around the
   // Y axis.
   struct transform rotateY;
   rotateY.T[0][0] = cos(theta);
   rotateY.T[0][2] = sin(theta);
   rotateY.T[2][0] = -sin(theta);
   rotateY.T[2][2] = cos(theta);
   return rotateY;
}
struct transform RotZ(double theta) {
   // Returns tranform for left multiplying
   // to a transform or point to
   // rotate theta *RADIANS* around the
   // Z axis.
   struct transform rotateZ;
   rotateZ.T[0][0] = cos(theta);
   rotateZ.T[0][1] = -sin(theta);
   rotateZ.T[1][0] = sin(theta);
   rotateZ.T[1][1] = cos(theta);
   return rotateZ;
}

void printmatrix(struct transform matrix) {
   fprintf(stderr, "Matrix contains:\n");
   fprintf(stderr, "%f %f %f %f\n", matrix.T[0][0], matrix.T[0][1], matrix.T[0][2], matrix.T[0][3]);
   fprintf(stderr, "%f %f %f %f\n", matrix.T[1][0], matrix.T[1][1], matrix.T[1][2], matrix.T[1][3]);
   fprintf(stderr, "%f %f %f %f\n", matrix.T[2][0], matrix.T[2][1], matrix.T[2][2], matrix.T[2][3]);
   fprintf(stderr, "%f %f %f %f\n", matrix.T[3][0], matrix.T[3][1], matrix.T[3][2], matrix.T[3][3]);
}