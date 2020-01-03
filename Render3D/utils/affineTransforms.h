#include <string.h>
#include <math.h>

#ifndef AFFINETRANSFORMS_H
#define AFFINETRANSFORMS_H
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

   point operator+(const point &a) const {
      return point(x + a.x, y + a.y, z + a.z);
   }
   point operator-(const point &a) const {
      return point(x - a.x, y - a.y, z - a.z);
   }

   point operator-() const {
      return point(-x, -y, -z);
   }

   point operator*(double scalar) const {
      return point(x * scalar, y * scalar, z * scalar);
   }
};

struct transform {
   // this struct defines a 4x4 matrix used for
   // 3d affine transforms in homogeneous coordinates
   double T[4][4];
   transform() {
      memset(&T[0][0], 0, 16 * sizeof(double));
      T[0][0] = 1;
      T[1][1] = 1;
      T[2][2] = 1;
      T[3][3] = 1;
   }
   transform operator*(const transform &b) const {
      struct transform result;
      for (int i = 0; i < 4; i++)
         for (int j = 0; j < 4; j++)
            result.T[i][j] = (T[i][0] * b.T[0][j]) + (T[i][1] * b.T[1][j]) + (T[i][2] * b.T[2][j]) + (T[i][3] * b.T[3][j]);

      return result;
   }
   transform &operator*=(const transform &b) {
      *this = b * *this;
      return *this;
   }
   point operator*(const point &p) const {
      struct point pr;
      pr.x = (T[0][0] * p.x) + (T[0][1] * p.y) + (T[0][2] * p.z) + (T[0][3] * p.w);
      pr.y = (T[1][0] * p.x) + (T[1][1] * p.y) + (T[1][2] * p.z) + (T[1][3] * p.w);
      pr.z = (T[2][0] * p.x) + (T[2][1] * p.y) + (T[2][2] * p.z) + (T[2][3] * p.w);
      return pr;
   }
   point operator*(const point *p) const {
      struct point pr;
      pr.x = (T[0][0] * p->x) + (T[0][1] * p->y) + (T[0][2] * p->z) + (T[0][3] * p->w);
      pr.y = (T[1][0] * p->x) + (T[1][1] * p->y) + (T[1][2] * p->z) + (T[1][3] * p->w);
      pr.z = (T[2][0] * p->x) + (T[2][1] * p->y) + (T[2][2] * p->z) + (T[2][3] * p->w);
      return pr;
   }
};

//Functions for returning affine transforms
struct transform I();
struct transform Sc(double Xscale, double Yscale, double Zscale);
struct transform Sc(double uniform_scale);
struct transform Tr(double Xtranslate, double Ytranslate, double Ztranslate);
struct transform RotX(double theta);
struct transform RotY(double theta);
struct transform RotZ(double theta);

void printmatrix(struct transform matrix);
#endif