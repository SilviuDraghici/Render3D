/*
   RayTracer code.
   Basecode by F. J. Estrada
  
   Uses Tom F. El-Maraghi's code for computing inverse
   matrices.

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

#include "utils.h" // <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

double iter_val = -1;
int num_frames = 24 * 5; //120 <- 5 frames

void buildScene(void)
{
#include "buildscene.c" // <-- Import the scene definition!
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct color *col)
{
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

   struct color tmp_col; // Accumulator for colour components
   struct color objcol;  // Colour for the object in R G and B

   // This will hold the colour as we process all the components of
   // the Phong illumination model
   tmp_col = 0;
   *col = 0;

   if (obj->texImg == NULL) // Not textured, use object colour
   {
      objcol = obj->col;
   }
   else
   {
      // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
      // for the object. Note that we will use textures also for Photon Mapping.
      
      //obj->textureMap(obj->texImg, a, b, &R, &G, &B);
   }

   //vector from intersection point to camera
   struct point3D c;
   c.px = -ray->d.px;
   c.py = -ray->d.py;
   c.pz = -ray->d.pz;
   c.pw = 1;
   normalize(&c);

   double ra = obj->alb.ra, rd = obj->alb.rd, rs = obj->alb.rs, rg = obj->alb.rg;
#ifdef SPECULAR
   ra = 0, rd = 0, rs = 0.5, rg = 0;
#endif
#ifdef DIFFUSE
   ra = 0.01, rs = 0, rg = 0;
#endif

   struct color ambient;
   ambient = objcol * ra;

   struct color diffuse;
   diffuse.R = 0, diffuse.G = 0, diffuse.B = 0;

   struct color specular;
   specular.R = 0, specular.G = 0, specular.B = 0;

   struct color global;

   struct pointLS *light = light_list;

   while (light != NULL)
   {

      //ray from intersection point to light source
      struct ray3D pToLight;
      // pToLight.p0.px = p->px - n->px*THR;
      // pToLight.p0.py = p->py - n->py*THR;
      // pToLight.p0.pz = p->pz - n->pz*THR;
      pToLight.p0.px = p->px;
      pToLight.p0.py = p->py;
      pToLight.p0.pz = p->pz;

      pToLight.p0.pw = 1;
      pToLight.d.px = light->p0.px - p->px;
      pToLight.d.py = light->p0.py - p->py;
      pToLight.d.pz = light->p0.pz - p->pz;
      pToLight.d.pw = 1;
      double light_lambda;
      //we don't care about these details of the hit
      struct object3D *dummyobj = NULL;
      struct point3D dummyp;
      double dummya;
#ifdef DEBUG
      if (0)
      {
         printf("light hit check:------\n");
         printf("startng from %c\n", obj->label);
         printf("p0: (%f, %f, %f)\n", p->px, p->py, p->pz);
      }
#endif
      findFirstHit(&pToLight, &light_lambda, obj, &dummyobj, &dummyp, &dummyp, &dummya, &dummya);
#ifdef DEBUG
      if (0)
      {
         if (dummyobj != NULL)
         {
            printf("intersected %c\n", dummyobj->label);
         }
         printf("light_lambda: %f\n", light_lambda);
      }
#endif

      if (!(THR < light_lambda && light_lambda < 1))
      {
         //vector from intersection point to light
         struct point3D s;
         s.px = light->p0.px - p->px;
         s.py = light->p0.py - p->py;
         s.pz = light->p0.pz - p->pz;
         s.pw = 1;
         normalize(&s);

         double n_dot_s = MAX(0, dot(n, &s));
         if (obj->frontAndBack)
         {
            n_dot_s = MAX(0, abs(dot(n, &s)));
         }
#ifdef DEBUG
         if (0)
         {
            printf("p: (%f, %f, %f)\n", p->px, p->py, p->pz);
            printf("l: (%f, %f, %f)\n", light->p0.px, light->p0.py, light->p0.pz);
            printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
            printf("s: (%f, %f, %f)\n", s.px, s.py, s.pz);
            printf("n_dot_s; %f\n", n_dot_s);
         }
#endif
         diffuse += objcol * rd * light->col * n_dot_s;

         //perfect light ray reflection
         struct ray3D light_reflection;
         normalize(&pToLight.d);
         pToLight.d.px = -pToLight.d.px;
         pToLight.d.py = -pToLight.d.py;
         pToLight.d.pz = -pToLight.d.pz;
         rayReflect(&pToLight, p, n, &light_reflection);
         double c_dot_m = dot(&c, &(light_reflection.d));

         specular += light->col * rs * pow(MAX(0, c_dot_m), obj->shinyness);

#ifdef REFLECTION

         col->R = MIN(1, (light_reflection.d.px + 1) * 0.5);
         col->G = MIN(1, (light_reflection.d.py + 1) * 0.5);
         col->B = MIN(1, (light_reflection.d.pz + 1) * 0.5);
//  fprintf(stderr, "R: %f, G: %f, B: %f\n", col->R, col->G, col->B);
#endif // DEBUG
      }

      light = light->next;
   }
   //perfect camera ray reflection
   struct ray3D cam_reflection;
   rayReflect(ray, p, n, &cam_reflection);
   //printf("r(%f, %f, %f), n(%f, %f, %f), ref(%f, %f, %f)\n",ray->d.px,ray->d.py,ray->d.pz, n->px,n->py,n->pz,cam_reflection.d.px,cam_reflection.d.py,cam_reflection.d.pz);
   rayTrace(&cam_reflection, depth + 1, &global, obj);
   if (global.R == -1)
   {
      global = 0;
   }
   else
   {
      global *= rg;
   }
#ifdef DEBUG
   if (1)
   {
      printf("obj r: %f  g: %f  b: %f\n", obj->col.R, obj->col.G, obj->col.B);
      printf("amb r: %f  g: %f  b: %f\n", ambient.R, ambient.G, ambient.B);
      printf("dif r: %f  g: %f  b: %f\n", diffuse.R, diffuse.G, diffuse.B);
      printf("fin r: %f  g: %f  b: %f\n", ambient.R + diffuse.R, ambient.G + diffuse.G, ambient.B + diffuse.B);
   }
#endif

#ifndef REFLECTION
   col->R = MIN(1, ambient.R + diffuse.R + specular.R + global.R);
   col->G = MIN(1, ambient.G + diffuse.G + specular.G + global.G);
   col->B = MIN(1, ambient.B + diffuse.B + specular.B + global.B);
#endif

#ifdef SIGNATURE
   *col = obj->col;
#endif

#ifdef NORMAL
   col->R = (n->px + 1) / 2;
   col->G = (n->py + 1) / 2;
   col->B = (n->pz + 1) / 2;
#endif
#ifdef REFLECTION_GLOBAL
   struct ray3D r;
   rayReflect(ray, p, n, &r);

   col->R = (r.d.px + 1) * 0.5;
   col->G = (r.d.py + 1) * 0.5;
   col->B = (r.d.pz + 1) * 0.5;
#endif // DEBUG

   return;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
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

   /////////////////////////////////////////////////////////////
   // TO DO: Implement this function. See the notes for
   // reference of what to do in here
   /////////////////////////////////////////////////////////////
   struct object3D *curr_obj = object_list;
   double curr_l, curr_a, curr_b;
   struct point3D curr_p, curr_n;
   *lambda = INFINITY;

   while (curr_obj != NULL)
   {
      curr_obj->intersect(curr_obj, ray, &curr_l, &curr_p, &curr_n, &curr_a, &curr_b);
#ifdef DEBUG
      if (0)
      {
         printf("checking %c\n", curr_obj->label);
         printf("fh lambda: %f\n", curr_l);
      }
#endif
      if (THR < curr_l && curr_l < *lambda)
      {
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

void rayTrace(struct ray3D *ray, int depth, struct color *col, struct object3D *Os)
{
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

   double lambda;               // Lambda at intersection
   double a, b;                 // Texture coordinates
   struct object3D *obj = NULL; // Pointer to object at intersection
   struct point3D p;            // Intersection point
   struct point3D n;            // Normal at intersection
   struct color I;          // Colour returned by shading function

   if (depth > MAX_DEPTH) // Max recursion depth reached. Return invalid colour.
   {
      col->R = -1;
      col->G = -1;
      col->B = -1;
      return;
   }
   ///////////////////////////////////////////////////////
   // TO DO: Complete this function. Refer to the notes
   // if you are unsure what to do here.
   ///////////////////////////////////////////////////////
   findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
   if (lambda == -1 || obj == NULL)
   {
      col->R = -1;
      col->G = -1;
      col->B = -1;
      return;
   }
   rtShade(obj, &p, &n, ray, depth, a, b, col);
}

int main(int argc, char *argv[])
{
   // Main function for the raytracer. Parses input parameters,
   // sets up the initial blank image, and calls the functions
   // that set up the scene and do the raytracing.
   struct image *im;       // Will hold the raytraced image
   struct view *cam;       // Camera and view for this scene
   int sx;                 // Size of the raytraced image
   int antialiasing;       // Flag to determine whether antialiaing is enabled or disabled
   char output_name[1024]; // Name of the output file for the raytraced .ppm image
   struct point3D e;       // Camera view parameters 'e', 'g', and 'up'
   struct point3D g;
   struct point3D up;
   double du, dv;               // Increase along u and v directions for pixel coordinates
   struct point3D pc, d;        // Point structures to keep the coordinates of a pixel and
                                // the direction or a ray
   struct ray3D ray;            // Structure to keep the ray from e to a pixel
   struct color col;        // Return colour for raytraced pixels
   struct color background; // Background colour
   int i, j;                    // Counters for pixel coordinates
   unsigned char *rgbIm;
   int depth;

   if (argc < 5)
   {
      fprintf(stderr, "RayTracer: Can not parse input parameters\n");
      fprintf(stderr, "USAGE: RayTracer size rec_depth antialias output_name\n");
      fprintf(stderr, "   size = Image size (both along x and y)\n");
      fprintf(stderr, "   rec_depth = Recursion depth\n");
      fprintf(stderr, "   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
      fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
      exit(0);
   }
   if (argc > 5)
   {
      iter_val = atof(argv[5]);
      printf("\nIteration: %f\n", iter_val);
   }
   sx = atoi(argv[1]);
   MAX_DEPTH = atoi(argv[2]);
   if (atoi(argv[3]) == 0)
      antialiasing = 0;
   else
      antialiasing = 1;
   strcpy(&output_name[0], argv[4]);

   //  fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
   //  fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
   //  if (!antialiasing) fprintf(stderr,"Antialising is off\n");
   //  else fprintf(stderr,"Antialising is on\n");
   //  fprintf(stderr,"Output file name: %s\n",output_name);

   object_list = NULL;
   light_list = NULL;
   texture_list = NULL;

   // Allocate memory for the new image
   im = newImage(sx, sx);
   if (!im)
   {
      fprintf(stderr, "Unable to allocate memory for raytraced image\n");
      exit(0);
   }
   else
      rgbIm = (unsigned char *)im->rgbdata;

   buildScene(); // Create a scene. This defines all the
                 // objects in the world of the raytracer

   // Mind the homogeneous coordinate w of all vectors below. DO NOT
   // forget to set it to 1, or you'll get junk out of the
   // geometric transformations later on.

   // Camera center is at (0,0,-1)
   e.px = 0;
   e.py = 0;
#ifdef CAMERA
   e.pz = -8;
#endif
#ifndef CAMERA
   e.pz = -1;
#endif
   e.pw = 1;

   // To define the gaze vector, we choose a point 'pc' in the scene that
   // the camera is looking at, and do the vector subtraction pc-e.
   // Here we set up the camera to be looking at the origin.
   g.px = 0 - e.px;
#ifdef CAMERA
   g.py = -3 - e.py;
#endif
#ifndef CAMERA
   g.py = -0 - e.py;
#endif
   g.pz = 0 - e.pz;
   g.pw = 1;
   // In this case, the camera is looking along the world Z axis, so
   // vector w should end up being [0, 0, -1]

   // Define the 'up' vector to be the Y axis
   up.px = 0;
   up.py = 1;
   up.pz = 0;
   up.pw = 1;

   // Set up view with given the above vectors, a 4x4 window,
   // and a focal length of -1 (why? where is the image plane?)
   // Note that the top-left corner of the window is at (-2, 2)
   // in camera coordinates.
   double focal = -1;
#ifdef CAMERA
   focal = -3;
#endif
   cam = setupView(&e, &g, &up, focal, -2, 2, 4);
   //cam=setupView(&e, &g, &up, -1, -2, 2, 4);

   if (cam == NULL)
   {
      fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
      cleanup(object_list, light_list, texture_list);
      deleteImage(im);
      exit(0);
   }

   // Set up background colour here
   background.R = 0;
   background.G = 0;
   background.B = 0;

   // Do the raytracing
   //////////////////////////////////////////////////////
   // TO DO: You will need code here to do the raytracing
   //        for each pixel in the image. Refer to the
   //        lecture notes, in particular, to the
   //        raytracing pseudocode, for details on what
   //        to do here. Make sure you undersand the
   //        overall procedure of raytracing for a single
   //        pixel.
   //////////////////////////////////////////////////////
   du = cam->wsize / (sx - 1);  // du and dv. In the notes in terms of wl and wr, wt and wb,
   dv = -cam->wsize / (sx - 1); // here we use wl, wt, and wsize. du=dv since the image is
                                // and dv is negative since y increases downward in pixel
                                // coordinates and upward in camera coordinates.

   //  fprintf(stderr,"View parameters:\n");
   //  fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
   //  fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
   //  printmatrix(cam->C2W);
   //  fprintf(stderr,"World to camera conversion matrix:\n");
   //  printmatrix(cam->W2C);
   //  fprintf(stderr,"\n");

   //fprintf(stderr,"Rendering\n");

#ifdef DEBUG
   for (j = 192; j < 320; j++)
   { // For each of the pixels in the image
      for (i = 192; i < 193; i++)
#else
   for (j = 0; j < sx; j++)
   { // For each of the pixels in the image
      for (i = 0; i < sx; i++)
#endif
      {
///////////////////////////////////////////////////////////////////
// TO DO - complete the code that should be in this loop to do the
//         raytracing!
///////////////////////////////////////////////////////////////////
#ifdef DEBUG
         printf("------------------------------------\n");
         printf("pixel: (%d %d)\n", i, j);
#endif
         pc.px = cam->wl + i * du;
         pc.py = cam->wt + j * dv;
         pc.pz = cam->f;
         pc.pw = 1;

         matVecMult(cam->C2W, &pc);

         ray.d.px = pc.px - cam->e.px;
         ray.d.py = pc.py - cam->e.py;
         ray.d.pz = pc.pz - cam->e.pz;
         ray.d.pw = 1;
         normalize(&ray.d);

         ray.p0.px = pc.px;
         ray.p0.py = pc.py;
         ray.p0.pz = pc.pz;
         ray.p0.pw = 1;
         depth = 1;
         col.R = 0;
         col.G = 0;
         col.B = 0;
         rayTrace(&ray, depth, &col, NULL);
         if (col.R == -1)
         {
            col = background;
         }
         rgbIm[3 * (j * im->sx + i)] = (unsigned char)(255 * col.R);
         rgbIm[3 * (j * im->sx + i) + 1] = (unsigned char)(255 * col.G);
         rgbIm[3 * (j * im->sx + i) + 2] = (unsigned char)(255 * col.B);
      } // end for i
   }    // end for j

   fprintf(stderr, "Done!\n");

   // Output rendered image
   imageOutput(im, output_name);

   // Exit section. Clean up and return.
   cleanup(object_list, light_list, texture_list); // Object, light, and texture lists
   deleteImage(im);                                // Rendered image
   free(cam);                                      // camera view
   exit(0);
}
