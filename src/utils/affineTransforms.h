#include <math.h>

#include "utils.h"

#ifndef AFFINETRANSFORMS_H
#define AFFINETRANSFORMS_H

//Functions for returning affine transforms
matrix Sc(double Xscale, double Yscale, double Zscale);
matrix Sc(double uniform_scale);
matrix Tr(double Xtranslate, double Ytranslate, double Ztranslate);
matrix RotX(double theta);
matrix RotY(double theta);
matrix RotZ(double theta);

void invert(double *T, double *Tinv);

void printmatrix(matrix &matrix);

void cosWeightedSample(struct point *n, struct point *d);
#endif