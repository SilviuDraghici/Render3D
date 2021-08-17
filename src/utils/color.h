#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

#ifndef COLOR_H
#define COLOR_H
/*
   The structure below defines an RGB colour, values are
   in [0,1]
*/
struct color {
    double R;
    double G;
    double B;
    
    color(double r = 0, double g = 0, double b = 0) {
        R = r;
        G = g;
        B = b;
    }
    color &operator=(const color &col) {
        R = col.R;
        G = col.G;
        B = col.B;
        return *this;
    }
    color &operator+=(const color &col) {
        R += col.R;
        G += col.G;
        B += col.B;
        return *this;
    }
    color &operator*=(const color &col) {
        R *= col.R;
        G *= col.G;
        B *= col.B;
        return *this;
    }
    color &operator/=(double div) {
        R /= div;
        G /= div;
        B /= div;
        return *this;
    }
    color &operator=(const double val) {
        R = val;
        G = val;
        B = val;
        return *this;
    }
    color operator*(const double scalar) const {
        return color(R * scalar, G * scalar, B * scalar);
    }
    color operator*=(const double scalar) {
        R *= scalar;
        G *= scalar;
        B *= scalar;
        return *this;
    }

    color operator*(const color &a) const {
        return color{R * a.R, G * a.G, B * a.B};
    }

    color operator+(const color &a) const {
        return color(R + a.R, G + a.G, B + a.B);
    }

    color inverse(){
        return color(1.0 - R, 1.0 - G, 1.0 - B);
    }
};

inline std::ostream &operator<<(std::ostream &strm, const color &a){
    return strm << "(" << a.R << ", " << a.G << ", " << a.B << ")";
}

#endif