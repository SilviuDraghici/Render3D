/*
   PathTracer header.
   Basecode by F. J. Estrada
  
   Uses Tom F. El-Maraghi's code for computing inverse
   matrices.

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#ifndef __PathTracer_header
#define __PathTracer_header

#define PI 3.14159265354
#define MAX_D 15

/* The structure below is used to hold a single RGB image. 'rgbdata' points */
/* to an array of size sx*sy*[num layers] but note that your code must know */
/* the data type for the array. Within the path tracer, there will be double  */
/* floating point images, and also unsigned char images, and there will be  */
/* 3-layer images (your path traced scene, texture and normal maps) as well   */
/* as 1-layer images (alpha maps)					      */
struct image
{
   void *rgbdata;
   int sx;
   int sy;
};

/* The structure below is used to hold 1 image that will be used as either  */
/* a texture map, a normal map, or an alpha map. This allows us to store    */
/* texture data separately from objects and avoid duplication when multiple */
/* objects share the same texture/normal/alpha maps                         */
struct textureNode
{
   char name[1024];
   int type;
   struct image *im;
   struct textureNode *next;
};

/* The structure below defines a point in 3D homogeneous coordinates        */
struct point3D
{
   double px;
   double py;
   double pz;
   double pw;
};

struct triangle
{
   int v1;
   int v2;
   int v3;
   struct point3D n;
   struct triangle *next;
};

/* The structure below defines a ray in 3D homogeneous coordinates,
   the point corresponds to the representation r(\lambda)=p+(\lambda)d */
struct ray3D
{
   struct point3D p0; // Ray origin (at lambda=0)
   struct point3D d;  // Ray direction
   void (*rayPos)(struct ray3D *ray, double lambda, struct point3D *pos);
   // Function to return the
   // position along the ray
   // for a given lambda.
   double H = 0.5;      // For monochromatic rays, HUE value (used to obtain colour, and as
                        //   a convenient substitute for wavelength) values in [0 1] go
                        //   from deep red to purple
   double R, G, B;      // Keeps track of the ray's colour as it
                        // bounces through the scene
   double Ir, Ig, Ib;   // Accumulators of brightness for explicit
                        // light sampling
   struct point3D srcN; // Source normal - needed for diffuse illumination
                        /* You may add data here to keep track of any values associated */
                        /* with this ray when implementing advanced raytracing features */

   int isLightRay = 0;
};

/*
   The structure below defines an RGB colour, values are
   in [0,1]
*/
struct color
{
   double R;
   double G;
   double B;
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
struct object3D
{
   struct color col; // Object's colour in RGB
   double T[4][4];       // T holds the transformation applied to this object.
   double Tinv[4][4];    // Tinv holds the inverse transformation
   char label[8];
   // Below we set up space for a pointer to the intersection function for this object.
   // Note that the intersection function must compute the lambda at the intersection, the
   // intersection point p, the normal at that point n, and the texture coordinates (a,b).
   // The texture coordinates are not used unless texImg!=NULL and a textureMap function
   // has been provided
   void (*intersect)(struct object3D *obj, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b);

   // Texture mapping function. Takes normalized texture coordinates (a,b) and returns the
   // texture colour at that point using bi-linear interpolation
   void (*textureMap)(struct image *img, double a, double b, double *R, double *G, double *B);

   // Functions to return coordinates on the surface of the object. One takes as input the a and b
   // parameters for the parametric function of the object and returns the (x,y,z) coordinates
   // on the object surface. The second returns a uniformly random-sampled point on the surface.
   // These are needed for Photon Mapping.
   void (*surfaceCoords)(struct object3D *obj, double a, double b, double *x, double *y, double *z);
   void (*randomPoint)(struct object3D *obj, double *x, double *y, double *z);

   struct image *texImg;       // Pointer to structure holding the texture for this object
   struct image *normalMap;    // Normal map for this object
   struct image *alphaMap;     // Alpha map for the object
   struct image *intersectMap; // Alpha map except where it's zero the intersection is ignored.

   // Material properties (no longer Phong model!)
   double diffPct; // Diffuse component proportion
   double reflPct; // Purely reflecting component proportion
   double tranPct; // Refracting component proportion
                   // NOTE: diffPct+reflPct+tranPct=1.0

   double refl_sig; // For reflective component, controls the amount of spread
   double r_index;  // Index of refraction

   int isLightSource; // Flag to indicate if this is an area light source
   double LSweight;   // If this is an area light source, keeps track of its weight (volume)

   int frontAndBack; // This flag indicates that normal direction should always
                     // point toward the viewer (e.g. as with planes)
   int isCSG;        // Object is part of a CSG composite object. Links to components via CSGnext

   struct object3D *CSGnext; // For CSG objects, points to next component
   struct object3D *next;    // Pointer to next entry in object linked list

   // You can modify this data structure to add any data/methods you
   // require to implement advanced features.
   struct triangle *triangles;
   struct point3D *vertices;
   struct point3D *normals;
   struct point3D *min;
   struct point3D *max;
};

/* The structure below defines a point light source */
struct pointLS
{
   struct color col; // Light source colour
   struct point3D p0;    // Light source location
   struct pointLS *next; // Pointer to next light in the scene
};

/*
   The structure below is used to hold camera parameters. You will need
   to write code to initialize the camera position and orientation.
*/
struct view
{
   struct point3D e; // Location of the camera center
   struct point3D u; // u vector
   struct point3D v; // v vector
   struct point3D w; // w vector
   double f;         // Focal length
   double wl;        // Left edge in camera coordinates
   double wt;        // Top edge in camera coordinates
   double wsize;     // Window size in distance units (not pixels!)
   double W2C[4][4]; // World2Camera conversion matrix
   double C2W[4][4]; // Camera2World conversion matrix
};

// Function definitions start here
int main(int argc, char *argv[]); // Main raytracing function.

// Raytracing
void buildScene(void); // Scene set up. Defines objects and object transformations
void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b);
void PathTrace(struct ray3D *ray, int depth, struct color *col, struct object3D *Os, struct object3D *explicit_l);

#endif
