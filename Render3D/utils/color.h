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
struct color{
   double R;
   double G;
   double B;
   color(double r = 0, double g = 0, double b = 0)
   {
      R = r;
      G = g;
      B = b;
   }
   color &operator=(const color &col)
   {
      R = col.R;
      G = col.G;
      B = col.B;
      return *this;
   }
   color &operator+=(const color &col)
   {
      R += col.R;
      G += col.G;
      B += col.B;
      return *this;
   }
   color &operator=(const double val)
   {
      R = val;
      G = val;
      B = val;
      return *this;
   }
   color operator*(const double scalar) const
   {
      return color(R * scalar, G * scalar, B * scalar);
   }
   color operator*=(const double scalar)
   {
      R *= scalar;
      G *= scalar;
      B *= scalar;
      return *this;
   }

   color operator*(const color &a) const
   {
      return color{R * a.R, G * a.G, B * a.B};
   }

   color operator+(const color &a) const
   {
      return color(R + a.R, G + a.G, B + a.B);
   }
};

/*
   The structures below are used to define an object colour in terms of the
   components of the Phong illumination model. Note that we typically
   define colours together with objects, so you should not need to
   instantiate lone instances of the colour structure.

   Also, note that you can easily make your objects completely white
   (or completely monochromatic) by not being careful how the different
   components in the Phong model add up. Take a moment and think how
   you want your object to look before you set these values.
*/
struct albedosPhong
{
   double ra; // Ambient light albedo
   double rd; // Diffuse component albedo
   double rs; // Specular component albedo
   double rg; // Global component albedo
};

#endif