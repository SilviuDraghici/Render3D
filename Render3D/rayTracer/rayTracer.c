#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rayTracer.h"

#include "../utils/affinetransforms.h"
#include "../utils/buildscene.h"
#include "../utils/camera.h"
#include "../utils/color.h"
#include "../utils/imageProcessor.h"
#include "../utils/mappings.h"
#include "../utils/ray.h"
#include "../utils/utils.h"

int MAX_DEPTH;

void rayTraceMain(int argc, char *argv[]) {
   if (argc < 5) {
      fprintf(stderr, "RayTracer: Can not parse input parameters\n");
      fprintf(stderr, "USAGE: Render3D 0 size rec_depth antialias output_name\n");
      fprintf(stderr, "   size = Image size (both along x and y)\n");
      fprintf(stderr, "   rec_depth = Recursion depth\n");
      fprintf(stderr, "   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
      fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
      exit(0);
   }

   int sx = atoi(argv[2]);
   MAX_DEPTH = atoi(argv[3]);
   strcpy(&output_name[0], argv[5]);

   unsigned char *rgbIm;
   // Allocate memory for the new image
   outImage = newImage(sx, sx);
   if (!outImage) {
      fprintf(stderr, "Unable to allocate memory for raytraced image\n");
      exit(0);
   } else {
      rgbIm = (unsigned char *)outImage->rgbdata;
   }

   // Camera center is at (0,0,-1)
   cam_pos = point(0, 0, -1);

   // To define the gaze vector, we choose a point 'pc' in the scene that
   // the camera is looking at, and do the vector subtraction pc-e.
   // Here we set up the camera to be looking at the origin.

   cam_gaze_point = point(0, 0, 0);
   cam_gaze = cam_gaze_point - cam_pos;

   cam_up = point(0, 1, 0);

   cam_focal = -1;

   buildScene();

   struct view *cam;  // Camera and view for this scene
   cam = setupView(&cam_pos, &cam_gaze, &cam_up, cam_focal, -2, 2, 4);

   if (cam == NULL) {
      fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
      //cleanup(object_list, light_list, texture_list);
      deleteImage(outImage);
      exit(0);
   }

   setPixelStep(cam, sx, sx);

   struct ray ray;
   double depth;
   color col;
   color background = 0;

   for (int j = 0; j < sx; j++) {  // For each of the pixels in the image
      for (int i = 0; i < sx; i++) {
         getRayFromPixel(&ray, cam, i, j);
         depth = 1;
         col = 0;
         rayTrace(&ray, depth, &col, NULL);
         if (col.R == -1) {
            col = background;
         }
         rgbIm[3 * (j * outImage->sx + i) + 0] = (unsigned char)(255 * col.R);
         rgbIm[3 * (j * outImage->sx + i) + 1] = (unsigned char)(255 * col.G);
         rgbIm[3 * (j * outImage->sx + i) + 2] = (unsigned char)(255 * col.B);
      }  // end for i
   }     // end for j

   fprintf(stderr, "Done!\n");

   // Output rendered image
   imageOutput(outImage, output_name);
}

void rayTrace(struct ray *ray, int depth, struct color *col, struct object *Os) {
   // Trace one ray through the scene.
   //
   // Parameters:
   //   *ray   -  A pointer to the ray being traced
   //   depth  -  Current recursion depth for recursive raytracing
   //   *col   - Pointer to an RGB colour structure so you can return the object colour
   //            at the intersection point of this ray with the closest scene object.
   //   *Os    - 'Object source' is a pointer to the object from which the ray
   //            originates so you can discard self-intersections due to numerical
   //            errors. NULL for rays originating from the center of projection.

   double lambda;              // Lambda at intersection
   double a, b;                // Texture coordinates
   struct object *obj = NULL;  // Pointer to object at intersection
   struct point p;             // Intersection point
   struct point n;             // Normal at intersection
   struct color I;             // Colour returned by shading function

   if (depth > MAX_DEPTH)  // Max recursion depth reached. Return invalid colour.
   {
      *col = -1;
      return;
   }
   findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
   if (lambda == -1 || obj == NULL) {
      *col = -1;
      return;
   }
   rtShade(obj, &p, &n, ray, depth, a, b, col);
}

void findFirstHit(struct ray *ray, double *lambda, struct object *Os, struct object **obj, struct point *p, struct point *n, double *a, double *b) {
   // Find the closest intersection between the ray and any objects in the scene.
   // Inputs:
   //   *ray    -  A pointer to the ray being traced
   //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
   //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
   //              projection
   // Outputs:
   //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
   //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
   //              this ray (this is required so you can do the shading)
   //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
   //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
   //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

   struct object *curr_obj = object_list;
   double curr_l, curr_a, curr_b;
   struct point curr_p, curr_n;
   *lambda = INFINITY;

   while (curr_obj != NULL) {
      curr_obj->intersect(curr_obj, ray, &curr_l, &curr_p, &curr_n, &curr_a, &curr_b);
      if (THR < curr_l && curr_l < *lambda) {
         *lambda = curr_l;
         *obj = curr_obj;
         *p = curr_p;
         *n = curr_n;
         *a = curr_a;
         *b = curr_b;
      }
      curr_obj = curr_obj->next;
   }
}

void rtShade(struct object *obj, struct point *p, struct point *n, struct ray *ray, int depth, double a, double b, struct color *col) {
   // This function implements the shading model as described in lecture. It takes
   // - A pointer to the first object intersected by the ray (to get the colour properties)
   // - The coordinates of the intersection point (in world coordinates)
   // - The normal at the point
   // - The ray (needed to determine the reflection direction to use for the global component, as well as for
   //   the Phong specular component)
   // - The current recursion depth
   // - The (a,b) texture coordinates (meaningless unless texture is enabled)
   //
   // Returns:
   // - The colour for this ray (using the col pointer)
   //

   struct color tmp_col;  // Accumulator for colour components
   struct color objcol;   // Colour for the object in R G and B

   // This will hold the colour as we process all the components of
   // the Phong illumination model
   tmp_col = 0;
   *col = 0;

   // set object color
   textureMap(obj, a, b, &objcol);

   // check for normal map
   normalMap(obj, a, b, n);

   //vector from intersection point to camera
   struct point c;
   c = -ray->d;
   normalize(&c);

   double ra = obj->alb.ra, rd = obj->alb.rd, rs = obj->alb.rs, rg = obj->alb.rg;

   struct color ambient;
   ambient = objcol * ra;

   struct color diffuse;

   struct color specular;

   struct color global;

   struct pointLS *light = light_list;

   while (light != NULL) {
      //ray from intersection point to light source
      struct ray pToLight;
      // pToLight.p0 = *p - n * THR;
      pToLight.p0 = *p;

      pToLight.d = light->p0 - *p;
      double light_lambda;
      //we don't care about these details of the hit
      struct object *dummyobj = NULL;
      struct point dummyp;
      double dummya;

      findFirstHit(&pToLight, &light_lambda, obj, &dummyobj, &dummyp, &dummyp, &dummya, &dummya);

      if (!(THR < light_lambda && light_lambda < 1)) {
         //vector from intersection point to light
         struct point s;
         s = light->p0 - *p;
         normalize(&s);

         double n_dot_s = MAX(0, dot(n, &s));
         if (obj->frontAndBack) {
            n_dot_s = MAX(0, abs(dot(n, &s)));
         }
         diffuse += objcol * rd * light->col * n_dot_s;

         //perfect light ray reflection
         struct ray light_reflection;
         normalize(&pToLight.d);
         pToLight.d = -pToLight.d;
         rayReflect(&pToLight, p, n, &light_reflection);
         double c_dot_m = dot(&c, &(light_reflection.d));

         specular += light->col * rs * pow(MAX(0, c_dot_m), obj->shinyness);

      }

      light = light->next;
   }
   //perfect camera ray reflection
   struct ray cam_reflection;
   rayReflect(ray, p, n, &cam_reflection);
   //printf("r(%f, %f, %f), n(%f, %f, %f), ref(%f, %f, %f)\n",ray->d.px,ray->d.py,ray->d.pz, n->px,n->py,n->pz,cam_reflection.d.px,cam_reflection.d.py,cam_reflection.d.pz);
   rayTrace(&cam_reflection, depth + 1, &global, obj);
   if (global.R == -1) {
      global = 0;
   } else {
      global *= rg;
   }

   col->R = MIN(1, ambient.R + diffuse.R + specular.R + global.R);
   col->G = MIN(1, ambient.G + diffuse.G + specular.G + global.G);
   col->B = MIN(1, ambient.B + diffuse.B + specular.B + global.B);
   return;
}
