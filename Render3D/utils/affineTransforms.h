#include <math.h>
#include "utils.h"

#ifndef AFFINETRANSFORMS_H
#define AFFINETRANSFORMS_H

//Functions for returning affine transforms
struct matrix I();
struct matrix Sc(double Xscale, double Yscale, double Zscale);
struct matrix Sc(double uniform_scale);
struct matrix Tr(double Xtranslate, double Ytranslate, double Ztranslate);
struct matrix RotX(double theta);
struct matrix RotY(double theta);
struct matrix RotZ(double theta);

void invert(double *T, double *Tinv);

void printmatrix(struct matrix matrix);
#endif