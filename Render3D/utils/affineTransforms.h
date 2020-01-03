#include <math.h>
#include "utils.h"

#ifndef AFFINETRANSFORMS_H
#define AFFINETRANSFORMS_H

//Functions for returning affine transforms
struct transform I();
struct transform Sc(double Xscale, double Yscale, double Zscale);
struct transform Sc(double uniform_scale);
struct transform Tr(double Xtranslate, double Ytranslate, double Ztranslate);
struct transform RotX(double theta);
struct transform RotY(double theta);
struct transform RotZ(double theta);

void invert(double *T, double *Tinv);

void printmatrix(struct transform matrix);
#endif