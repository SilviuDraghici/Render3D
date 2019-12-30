/*
   RayTracer header file.
   Basecode by F. J. Estrada

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __RayTracer_header
#define __RayTracer_header

#define PI 3.14159265354

/* The structure below is used to hold a single RGB image. 'rgbdata' points */
/* to an array of size sx*sy*[num layers] but note that your code must know */
/* the data type for the array. Within the raytracer, there will be double  */
/* floating point images, and also unsigned char images, and there will be  */
/* 3-layer images (your raytraced scene, texture and normal maps) as well   */
/* as 1-layer images (alpha maps)					      */
struct image {
   void *rgbdata;
   int sx;
   int sy;
};

/* The structure below is used to hold 1 image that will be used as either  */
/* a texture map, a normal map, or an alpha map. This allows us to store    */
/* texture data separately from objects and avoid duplication when multiple */
/* objects share the same texture/normal/alpha maps                         */
struct textureNode {
   char name[1024];
   int type;
   struct image *im;
   struct textureNode *next;
};

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
   transform operator*(const transform &b) const{
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
   point operator*(const point &p) const{
      struct point pr;
      pr.x = (T[0][0] * p.x) + (T[0][1] * p.y) + (T[0][2] * p.z) + (T[0][3] * p.w);
      pr.y = (T[1][0] * p.x) + (T[1][1] * p.y) + (T[1][2] * p.z) + (T[1][3] * p.w);
      pr.z = (T[2][0] * p.x) + (T[2][1] * p.y) + (T[2][2] * p.z) + (T[2][3] * p.w);
      return pr;
   }
   point operator*(const point *p) const{
      struct point pr;
      pr.x = (T[0][0] * p->x) + (T[0][1] * p->y) + (T[0][2] * p->z) + (T[0][3] * p->w);
      pr.y = (T[1][0] * p->x) + (T[1][1] * p->y) + (T[1][2] * p->z) + (T[1][3] * p->w);
      pr.z = (T[2][0] * p->x) + (T[2][1] * p->y) + (T[2][2] * p->z) + (T[2][3] * p->w);
      return pr;
   }
};

struct triangle {
   int v1;
   int v2;
   int v3;
   struct point n;
   struct triangle *next;
};

/* The structure below defines a ray in 3D homogeneous coordinates,
   the point corresponds to the representation r(\lambda)=p+(\lambda)d */
struct ray3D {
   struct point p0;  // Ray origin (at lambda=0)
   struct point d;   // Ray direction
   void (*rayPos)(struct ray3D *ray, double lambda, struct point *pos);
   // Function to return the
   // position along the ray
   // for a given lambda.

   /* You may add data here to keep track of any values associated */
   /* with this ray when implementing advanced raytracing features */
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
struct albedosPhong {
   double ra;  // Ambient light albedo
   double rd;  // Diffuse component albedo
   double rs;  // Specular component albedo
   double rg;  // Global component albedo
};

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
};

/*
   The structure below defines an Object within the World Coordinate Frame.
   For this ray tracer, we will use a simple linked list of objects (not
   a tree structure, so hierarchical transformations have to be fully
   specified). This is not really limiting and simplifies the data
   structure to hold the scene.

   Note that the Object3D structure is completely agnostic about the type
   of object defined, the difference is made by the intersect function.
   If you think carefully about it, in terms of raytracing objects are
   defined entirely by their surface, and the surface is defined in the
   intersect function for each object.

   Thus, to create additional objects, simply provide a suitable
   intersection function. The rest stays the same.
*/
struct object3D {
   struct albedosPhong alb;  // Object's albedos for Phong model
   struct color col;         // Object's colour in RGB
   struct transform T;       // T holds the transformation applied to this object.
   struct transform Tinv;    // Tinv holds the inverse transformation

   char label;

   // Below we set up space for a pointer to the intersection function for this object.
   // Note that the intersection function must compute the lambda at the intersection, the
   // intersection point p, the normal at that point n, and the texture coordinates (a,b).
   // The texture coordinates are not used unless texImg!=NULL and a textureMap function
   // has been provided
   void (*intersect)(struct object3D *obj, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b);

   // Texture mapping function. Takes normalized texture coordinates (a,b) and returns the
   // texture colour at that point using bi-linear interpolation
   void (*textureMap)(struct image *img, double a, double b, double *R, double *G, double *B);

   // Functions to return coordinates on the surface of the object. One takes as input the a and b
   // parameters for the parametric function of the object and returns the (x,y,z) coordinates
   // on the object surface. The second returns a uniformly random-sampled point on the surface.
   // These are needed for Photon Mapping.
   void (*surfaceCoords)(struct object3D *obj, double a, double b, double *x, double *y, double *z);
   void (*randomPoint)(struct object3D *obj, double *x, double *y, double *z);
   struct image *texImg;     // Pointer to structure holding the texture for this object
   struct image *photonMap;  // Photon map for this object
   struct image *normalMap;  // Normal map for this object
   struct image *alphaMap;   // Alpha map for the object

   // Material properties
   double alpha;       // Opacity - if less than 1 this is a semi transparent object and refraction rays
                       // should be implemented
   double r_index;     // Index of refraction
   double shinyness;   // Exponent for the Phong specular component
   int frontAndBack;   // Flag to indicate that both sides of the object
                       // should be lit.
   int isLightSource;  // Flag to indicate if this is an area light source
   int isCSG;          // Object is part of a CSG composite object. Links to components via CSGnext
   int photonMapped;   // This object accumulates photons under photon mapping
   int normalMapped;   // This object has an associated normal map
   int alphaMapped;    // This object has an associated alpha map

   struct object3D *CSGnext;  // For CSG objects, points to next component
   struct object3D *next;     // Pointer to next entry in object linked list

   // If needed for the advanced raytracer, you can modify this data structure to add any data/methods you
   // require.
   struct triangle *triangles;
   struct point *vertices;
};

/* The structure below defines a point light source */
struct pointLS {
   struct color col;      // Light source colour
   struct point p0;       // Light source location
   struct pointLS *next;  // Pointer to next light in the scene
};

/*
   The structure below is used to hold camera parameters. You will need
   to write code to initialize the camera position and orientation.
*/
struct view {
   struct point e;    // Location of the camera center
   struct point u;    // u vector
   struct point v;    // v vector
   struct point w;    // w vector
   double f;          // Focal length
   double wl;         // Left edge in camera coordinates
   double wt;         // Top edge in camera coordinates
   double wsize;      // Window size in distance units (not pixels!)
   double W2C[4][4];  // World2Camera conversion matrix
   double C2W[4][4];  // Camera2World conversion matrix
};

// Function definitions start here
int main(int argc, char *argv[]);  // Main raytracing function.

// Raytracing
void buildScene(void);                                                                // Scene set up. Defines objects and object transformations
void rayTrace(struct ray3D *ray, int depth, struct color *col, struct object3D *Os);  // RayTracing routine
void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point *p, struct point *n, double *a, double *b);
void rtShade(struct object3D *obj, struct point *p, struct point *n, struct ray3D *ray, int depth, double a, double b, struct color *col);

#endif
