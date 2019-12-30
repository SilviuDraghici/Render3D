/*
   RayTracer code.
   Basecode by F. J. Estrada
  
   Uses Tom F. El-Maraghi's code for computing inverse
   matrices.

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

//  __  __           _         __                             _
// |  \/  |_   _ ___| |_ __ _ / _| __ _  __      ____ _ ___  | |__   ___ _ __ ___
// | |\/| | | | / __| __/ _` | |_ / _` | \ \ /\ / / _` / __| | '_ \ / _ \ '__/ _ \
// | |  | | |_| \__ \ || (_| |  _| (_| |  \ V  V / (_| \__ \ | | | |  __/ | |  __/
// |_|  |_|\__,_|___/\__\__,_|_|  \__,_|   \_/\_/ \__,_|___/ |_| |_|\___|_|  \___|
//
//     _    ____   _           _ _        _   _
//    / \  |  _ \ / \      ___(_) |_ __ _| |_(_) ___  _ __
//   / _ \ | |_) / _ \    / __| | __/ _` | __| |/ _ \| '_ \
//  / ___ \|  __/ ___ \  | (__| | || (_| | |_| | (_) | | | |
// /_/   \_\_| /_/   \_\  \___|_|\__\__,_|\__|_|\___/|_| |_|

#include "utils_path.h" // <-- This includes PathTracer.h
#define __USE_IS        // Use importance sampling for diffuse materials
#define __USE_ES        // Use explicit light sampling
//#define __USE_DISP      // Give rays random colors
//#define __USE_BD
//#define __DEBUG		  	// <-- Use this to turn on/off debugging output

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;

//array of light sources
struct object3D **light_list;
int num_lights;
int curr_light;

double pct;

struct textureNode *texture_list;
unsigned long int NUM_RAYS;
int MAX_DEPTH;

#include "buildScene.c" // Import scene definition

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point *p, struct point *n, double *a, double *b)
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
   struct point curr_p, curr_n;
   *lambda = INFINITY;

   while (curr_obj != NULL)
   {
      curr_obj->intersect(curr_obj, ray, &curr_l, &curr_p, &curr_n, &curr_a, &curr_b);

      int IntersectMapped = (curr_obj->intersectMap != NULL);
      double intMapVal;
      if (IntersectMapped)
      {
         alphaMap(curr_obj->intersectMap, curr_a, curr_b, &intMapVal);
      }

      if (THR < curr_l && curr_l < *lambda)
      {
         if (!IntersectMapped || IntersectMapped && intMapVal >= 0.5)
         {
            *lambda = curr_l;
            *obj = curr_obj;
            *p = curr_p;
            *n = curr_n;
            *a = curr_a;
            *b = curr_b;
         }
         else
         {
            //check if there sould be intersection with object's back face
            struct point op = ray->p0;
            ray->p0.x = curr_p.x + THR * ray->d.x;
            ray->p0.y = curr_p.y + THR * ray->d.y;
            ray->p0.z = curr_p.z + THR * ray->d.z;
            curr_obj->intersect(curr_obj, ray, &curr_l, &curr_p, &curr_n, &curr_a, &curr_b);
            alphaMap(curr_obj->intersectMap, curr_a, curr_b, &intMapVal);
            if (THR < curr_l && curr_l < *lambda && intMapVal >= 0.5)
            {
               *lambda = curr_l;
               *obj = curr_obj;
               *p = curr_p;
               *n = curr_n;
               *a = curr_a;
               *b = curr_b;
            }
            ray->p0 = op;
         }
      }
      curr_obj = curr_obj->next;
   }
}

void PathTrace(struct ray3D *ray, int depth, struct color *col, struct object3D *Os, struct object3D *explicit_l)
{
   // Trace one light path through the scene.
   //
   // Parameters:
   //   *ray   -  A pointer to the ray being traced
   //   depth  -  Current recursion depth for recursive raytracing
   //   *col   - Pointer to an RGB colour structure so you can return the object colour
   //            at the intersection point of this ray with the closest scene object.
   //   *Os    - 'Object source' is a pointer to the object from which the ray
   //            originates so you can discard self-intersections due to numerical
   //            errors. NULL for rays originating from the center of projection.
   NUM_RAYS++;
   double lambda;               // Lambda at intersection
   double a, b;                 // Texture coordinates
   struct object3D *obj = NULL; // Pointer to object at intersection
   struct point p;            // Intersection point
   struct point n;            // Normal at intersection
   struct point d;
   double R = 0, G = 0, B = 0; // Handy in case you need to keep track of some RGB colour value
   double diffPct, reflPct, tranPct;
   double R_Shlick = 1;
   double dice; // Handy to keep a random value
   double max_col;

   struct object3D *explt = explicit_l;

   if (depth > MAX_DEPTH) // Max recursion depth reached. Return black (no light coming into pixel from this path).
   {
#ifdef __DEBUG
      printf("max depth\n");
#endif
      col->R = ray->Ir; // These are accumulators, initialized at 0. Whenever we find a source of light these
      col->G = ray->Ig; // get incremented accordingly. At the end of the recursion, we return whatever light
      col->B = ray->Ib; // we accumulated into these three values.
      return;
   }

#ifdef __DEBUG
   printf("----------------------| Depth: %d |----------------------\n", depth);
#endif
   ///////////////////////////////////////////////////////
   // TO DO: Complete this function. Refer to the notes
   // if you are unsure what to do here.
   ///////////////////////////////////////////////////////
   findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);
   if (obj == NULL || lambda < THR)
   {
#ifdef __DEBUG
      printf("obj null or lambda < thr\n");
#endif
      col->R = ray->Ir;
      col->G = ray->Ig;
      col->B = ray->Ib;
      return;
   }

   if (obj->texImg == NULL)
   { // Not textured, use object colour
      R = obj->col.R;
      G = obj->col.G;
      B = obj->col.B;
   }
   else
   {
      obj->textureMap(obj->texImg, a, b, &R, &G, &B);
   }

   if (obj->normalMap != NULL)
   {
      double n_delta_x = 0, n_delta_y = 0, n_delta_z = 0;
      obj->textureMap(obj->normalMap, a, b, &n_delta_x, &n_delta_y,
                      &n_delta_z);
      n.x -= (2 * n_delta_x - 1);
      n.y -= (2 * n_delta_y - 1);
      n.z -= (2 * n_delta_z - 1);
      normalize(&n);
      n.w = 1;
   }

   if (obj->alphaMap == NULL)
   {
      diffPct = obj->diffPct;
      reflPct = obj->reflPct;
      tranPct = obj->tranPct;
   }
   else
   {
      alphaMap(obj->alphaMap, a, b, &diffPct);
      tranPct = 1 - diffPct;

      reflPct = obj->reflPct;
      //scale down so total is still 1
      diffPct *= (1 - reflPct);
      tranPct *= (1 - reflPct);
   }

   //all the other terms use this
   ray->R *= R;
   ray->G *= G;
   ray->B *= B;

   // if hit light source
   if (obj->isLightSource)
   {
#ifdef __DEBUG
      printf("hit light\n");
#endif

      col->R = ray->R + ray->Ir;
      col->G = ray->G + ray->Ig;
      col->B = ray->B + ray->Ib;
#ifdef __USE_ES
      if (ray->isLightRay == 0)
      {
         if (explt == obj)
         { // the ray cast is the same as explicit
            col->R = ray->Ir;
            col->G = ray->Ig;
            col->B = ray->Ib;
         }
      }
#endif
      return;
   }

   //make sure normal is on side of incoming ray
   if (obj->frontAndBack && dot(&ray->d, &n) > 0)
   {
      //printf("flip\n");
      n.x *= -1;
      n.y *= -1;
      n.z *= -1;
   }

   memcpy(&ray->p0, &p, sizeof(struct point));

   // ******************** Importance sampling ************************
   dice = drand48();
   if (dice <= diffPct)
   { //diffuse
#ifdef __USE_IS
      cosWeightedSample(&n, &d);
      ray->d = d;
#else
      //******************** Random sampling ************************
      hemiSphereCoordinates(&n, &(d.px), &(d.py), &(d.pz));

      double d_dot_n = dot(&d, &n);
      ray->d.px = d.px * d_dot_n;
      ray->d.py = d.py * d_dot_n;
      ray->d.pz = d.pz * d_dot_n;
#endif

#ifdef __USE_ES // explicit light sampling
      if (ray->isLightRay == 0)
      {
         // ray from intersection point to light source
         struct ray3D pToLight;
         pToLight.p0.x = p.x;
         pToLight.p0.y = p.y;
         pToLight.p0.z = p.z;
         pToLight.p0.w = 1;

         double prob = 0;
         dice = drand48();
         for (int i = 0; i < num_lights; i++)
         {
            prob += light_list[i]->LSweight;
            if (dice < prob)
            {
               curr_light = i;
               break;
            }
         }

         double x, y, z;
         light_list[curr_light]->randomPoint(light_list[curr_light], &x, &y, &z);
         pToLight.d.x = x - p.x; //0 - p.px;
         pToLight.d.y = y - p.y; //9.95 - p.py;
         pToLight.d.z = z - p.z; //5 - p.pz;
         pToLight.d.w = 1;

         double light_lambda;
         struct object3D *obstruction = NULL;
         struct point lightp;
         struct point nls;
         double La, Lb;

         findFirstHit(&pToLight, &light_lambda, obj, &obstruction, &lightp, &nls,
                      &La, &Lb);

         explt = light_list[curr_light];

#ifdef __DEBUG
         printf("source %s\n", obj->label);
         if (obstruction == NULL)
         {
            printf("obst == NULL\n");
         }
         else
         {
            printf("obst %s\n", obstruction->label);
         }
#endif

         if (obstruction == light_list[curr_light])
         {
            double A = pct * light_list[curr_light]->LSweight;
            double dxd = pToLight.d.x * pToLight.d.x + pToLight.d.y * pToLight.d.y + pToLight.d.pzz pToLight.d.pzz
            normalize(&pToLight.d);
            double n_dot_l = fabs(dot(&n, &pToLight.d));
            double nls_dot_l = fabs(dot(&nls, &pToLight.d));
            double w = MIN(1, (A * n_dot_l * nls_dot_l) / (dxd));

            double LR, LG, LB;
            if (light_list[curr_light]->texImg == NULL)
            { // Not textured, use object colour
               LR = light_list[curr_light]->col.R;
               LG = light_list[curr_light]->col.G;
               LB = light_list[curr_light]->col.B;
            }
            else
            {
               light_list[curr_light]->textureMap(light_list[curr_light]->texImg, La, Lb, &LR, &LG, &LB);
            }
            ray->Ir += w * LR * ray->R;
            ray->Ig += w * LG * ray->G;
            ray->Ib += w * LB * ray->B;

#ifdef __DEBUG
            if (1)
            {
               printf("ray R: %f, G: %f, B: %f\n", ray->R, ray->G, ray->B);
               printf("light R: %f, G: %f, B: %f\n", light_list[curr_light]->col.R, light_list[curr_light]->col.G, light_list[curr_light]->col.B);
               printf("intersection point: (%f,%f,%f)\n", p.px, p.py, p.pz);
               printf("light point: (%f,%f,%f)\n", lightp.px, lightp.py, lightp.pz);
               printf("w: %f a: %f n_dot_l %f nls_dot_l: %f dxd: %f\n", w, A, n_dot_l, nls_dot_l, dxd);
               printf("R: %f, G: %f, R: %f\n", ray->Ir, ray->Ig, ray->Ib);
            }
#endif
         }
#ifdef __USE_BD
         else
         {
            struct colourRGB EScol;
            biDirectionalSample(light_list[curr_light], ray, &EScol);
            ray->Ir += EScol.R * ray->R;
            ray->Ig += EScol.G * ray->G;
            ray->Ib += EScol.B * ray->B;
         }
#endif
      }
#endif // end of explicit light sampling
   }
   else if (dice <= diffPct + reflPct)
   { //reflective
      explt = NULL;
      struct ray3D ray_reflected;
      rayReflect(ray, &p, &n, &ray_reflected);
      //burnished reflection
      ray->d.x = randn(ray_reflected.d.x, obj->refl_sig);
      ray->d.y = randn(ray_reflected.d.y, obj->refl_sig);
      ray->d.z = randn(ray_reflected.d.z, obj->refl_sig);
      normalize(&ray->d);
   }
   else
   { //refractive
      double r_index = obj->r_index;
#ifdef __USE_DISP
      RGB2Hue(&ray->H, ray->R, ray->G, ray->B);
      r_index = hueToRef(ray->H, r_index);
#endif
      explt = NULL;
      double r, n1 = 1, n2 = 1, theta = -1;
      double c = dot(&n, &(ray->d));
      // moving from outside object to inside the object
      if (c > 0)
      {
         n.x *= -1;
         n.y *= -1;
         n.z *= -1;
         n1 = r_index;
         theta = c;
      }
      else
      {
         n2 = r_index;
         c *= -1;
         theta = c;
      }
      r = n1 / n2;
      //theta = c;

      struct ray3D refractRay;
      memcpy(&refractRay, ray, sizeof(struct ray3D));
      refractRay.p0.x = p.x - n.x * THR;
      refractRay.p0.y = p.y - n.y * THR;
      refractRay.p0.z = p.z - n.z * THR;
      refractRay.p0.w = 1;
      double s = 1 - (r * r) * (1 - (c * c));
      if (s > 0)
      {
         refractRay.d.x = ray->d.x * r + (r * c - sqrt(s)) * n.x;
         refractRay.d.y = ray->d.y * r + (r * c - sqrt(s)) * n.y;
         refractRay.d.z = ray->d.z * r + (r * c - sqrt(s)) * n.z;
         refractRay.d.w = 1;
         normalize(&refractRay.d);
         // Use Shlick's to figure out amount of reflected and refracted light
         double R0 = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2));
         // R(theta)=R0+((1-R0)*(1-cos(theta))^5)
         R_Shlick = MIN(1, R0 + (1 - R0) * pow(1 - theta, 5));
         dice = drand48();
         // randomly choose if refract or reflect on each iteration
         if (dice < (1 - R_Shlick))
         {
            dice = drand48();

            max_col = MAX(MAX(refractRay.R, refractRay.G), MAX(refractRay.R, refractRay.B));

            if (sqrt(dice) < max_col)
            {
               return PathTrace(&refractRay, depth + 1, col, obj, explt);
            }
         }
      }
      // reflect - total internal reflection or reflecting based on dice

      struct ray3D ray_reflected;
      rayReflect(ray, &p, &n, &ray_reflected);
      ray->d = ray_reflected.d;
   }

   dice = drand48();
   max_col = MAX(MAX(ray->R, ray->G), MAX(ray->R, ray->B));
   if (sqrt(dice) < max_col)
   {
      return PathTrace(ray, depth + 1, col, obj, explt);
   }
   else
   {
#ifdef __DEBUG
      printf("some russian stuff!\n");
#endif
      col->R = ray->Ir;
      col->G = ray->Ig;
      col->B = ray->Ib;
      return;
   }
}

int main(int argc, char *argv[])
{
   // Main function for the path tracer. Parses input parameters,
   // sets up the initial blank image, and calls the functions
   // that set up the scene and do the raytracing.
   struct image *im;       // Will hold the final image
   struct view *cam;       // Camera and view for this scene
   int sx;                 // Size of the  image
   int num_samples;        // Number of samples to use per pixel
   char output_name[1024]; // Name of the output file for the .ppm image file
   struct point e;       // Camera view parameters 'e', 'g', and 'up'
   struct point g;
   struct point up;
   double du, dv;        // Increase along u and v directions for pixel coordinates
   struct point pc, d; // Point structures to keep the coordinates of a pixel and
                         // the direction or a ray
   struct ray3D ray;     // Structure to keep the ray from e to a pixel
   struct color col; // Return colour for pixels
   int i, j, k;          // Counters for pixel coordinates and samples
   double *rgbIm;        // Image is now double precision floating point since we
                         // will be accumulating brightness differences with a
                         // wide dynamic range
   struct object3D *obj; // Will need this to process lightsource weights
   double *wght;         // Holds weights for each pixel - to provide log response
   double wt;

   time_t t1, t2;
   FILE *f;

   if (argc < 5)
   {
      fprintf(stderr, "PathTracer: Can not parse input parameters\n");
      fprintf(stderr, "USAGE: PathTracer size rec_depth num_samples output_name\n");
      fprintf(stderr, "   size = Image size (both along x and y)\n");
      fprintf(stderr, "   rec_depth = Recursion depth\n");
      fprintf(stderr, "   num_samples = Number of samples per pixel\n");
      fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
      exit(0);
   }
   sx = atoi(argv[1]);
   MAX_DEPTH = atoi(argv[2]);
   num_samples = atoi(argv[3]);
   strcpy(&output_name[0], argv[4]);

   fprintf(stderr, "Rendering image at %d x %d\n", sx, sx);
   fprintf(stderr, "Recursion depth = %d\n", MAX_DEPTH);
   fprintf(stderr, "Number of samples = %d\n", num_samples);
   fprintf(stderr, "Output file name: %s\n", output_name);

   object_list = NULL;
   texture_list = NULL;

   // Allocate memory for the new image
   im = newImage(sx, sx);
   wght = (double *)calloc(sx * sx, sizeof(double));
   if (!im || !wght)
   {
      fprintf(stderr, "Unable to allocate memory for image\n");
      exit(0);
   }
   else
      rgbIm = (double *)im->rgbdata;
   for (i = 0; i < sx * sx; i++)
      *(wght + i) = 1.0;

   buildScene(); // Create a scene.

   //count number of lights
   num_lights = 0;
   struct object3D *curr_obj = object_list;
   while (curr_obj != NULL)
   {
      if (curr_obj->isLightSource)
      {
         num_lights++;
      }
      curr_obj = curr_obj->next;
   }

   //create array of light pointers
   light_list = (struct object3D **)malloc(num_lights * sizeof(struct object3D *));
   num_lights = 0;
   curr_obj = object_list;
   while (curr_obj != NULL)
   {
      if (curr_obj->isLightSource)
      {
         light_list[num_lights] = curr_obj;
         num_lights++;
      }
      curr_obj = curr_obj->next;
   }
   curr_light = 0;

   // printf("Num Lights: %d\n", num_lights);
   // for(int i = 0; i<num_lights; i++){
   //   printf("light %s\n", light_list[i]->label);
   // }

   // Mind the homogeneous coordinate w of all vectors below. DO NOT
   // forget to set it to 1, or you'll get junk out of the
   // geometric transformations later on.

   // Camera center
   e.x = 0;
   e.y = 0;
   e.z = -15;
   e.w = 1;

   // To define the gaze vector, we choose a point 'pc' in the scene that
   // the camera is looking at, and do the vector subtraction pc-e.
   // Here we set up the camera to be looking at the origin.
   g.x = 0 - e.x;
   g.y = 0 - e.y;
   g.z = 0 - e.z;
   g.w = 1;
   // In this case, the camera is looking along the world Z axis, so
   // vector w should end up being [0, 0, -1]

   // Define the 'up' vector to be the Y axis
   up.x = 0;
   up.y = 1;
   up.z = 0;
   up.w = 1;

   // Set up view with given the above vectors, a 4x4 window,
   // and a focal length of -1 (why? where is the image plane?)
   // Note that the top-left corner of the window is at (-2, 2)
   // in camera coordinates.
   cam = setupView(&e, &g, &up, -3, -2, 2, 4);

   if (cam == NULL)
   {
      fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
      cleanup(object_list, texture_list);
      deleteImage(im);
      exit(0);
   }

   du = cam->wsize / (sx - 1);  // du and dv. In the notes in terms of wl and wr, wt and wb,
   dv = -cam->wsize / (sx - 1); // here we use wl, wt, and wsize. du=dv since the image is
                                // and dv is negative since y increases downward in pixel
                                // coordinates and upward in camera coordinates.

   //fprintf(stderr, "View parameters:\n");
   //fprintf(stderr, "Left=%f, Top=%f, Width=%f, f=%f\n", cam->wl, cam->wt, cam->wsize, cam->f);
   //fprintf(stderr, "Camera to world conversion matrix (make sure it makes sense!):\n");
   //printmatrix(cam->C2W);
   //fprintf(stderr, "World to camera conversion matrix:\n");
   //printmatrix(cam->W2C);
   //fprintf(stderr, "\n");

   // Update light source weights - will give you weights for each light source that add up to 1
   obj = object_list;
   pct = 0;
   while (obj != NULL)
   {
      if (obj->isLightSource)
         pct += obj->LSweight;
      obj = obj->next;
   }
   obj = object_list;
   while (obj != NULL)
   {
      if (obj->isLightSource)
      {
         obj->LSweight /= pct;
      }
      obj = obj->next;
   }
   fprintf(stderr, "\n");

   NUM_RAYS = 0;

   t1 = time(NULL);

   fprintf(stderr, "Rendering pass... ");
   for (k = 0; k < num_samples; k++)
   {
      fprintf(stderr, "%d/%d, ", k, num_samples);
#ifdef __DEBUG
      for (j = 0; j < sx; j++)
      { // For each of the pixels in the image
         for (i = 76; i < 77; i++)
         {
#else
#pragma omp parallel for schedule(dynamic, 1) private(i, j, pc, wt, ray, col, d)
      for (j = 0; j < sx; j++)
      { // For each of the pixels in the image
         for (i = 0; i < sx; i++)
         {
#endif

#ifdef __DEBUG
            printf("\n----------------------| PIXEL: (%d %d) |----------------------\n", i, j);
#endif
            col.R = 0;
            col.G = 0;
            col.B = 0;
            // Random sample within the pixel's area
            pc.x = (cam->wl + ((i + (drand48() - .5)) * du));
            pc.y = (cam->wt + ((j + (drand48() - .5)) * dv));
            pc.z = cam->f;
            pc.w = 1;

            // Convert image plane sample coordinates to world coordinates
            matVecMult(cam->C2W, &pc);

            // Now compute the ray direction
            memcpy(&d, &pc, sizeof(struct point));
            subVectors(&cam->e, &d); // Direction is d=pc-e
            normalize(&d);

            // Create a ray and do the raytracing for this pixel.
            initRay(&ray, &pc, &d);
            ray.isLightRay = 0;
#ifdef __USE_DISP
            ray.H = drand48();
            hue2RGB(ray.H, &(ray.R), &(ray.G), &(ray.B));
#endif

            wt = *(wght + i + (j * sx));
            PathTrace(&ray, 1, &col, NULL, NULL);

#ifdef __DEBUG
            printf("col returned: R: %f, G: %f, B: %f\n", col.R, col.G, col.B);
#endif

            (*(rgbIm + ((i + (j * sx)) * 3) + 0)) += col.R * pow(2, -log(wt));
            (*(rgbIm + ((i + (j * sx)) * 3) + 1)) += col.G * pow(2, -log(wt));
            (*(rgbIm + ((i + (j * sx)) * 3) + 2)) += col.B * pow(2, -log(wt));
            wt += col.R;
            wt += col.G;
            wt += col.B;
            *(wght + i + (j * sx)) = wt;
         } // end for i
      }    // end for j
      if (k % 25 == 0)
         dataOutput(rgbIm, sx, &output_name[0]); // Update output image every 25 passes
   }                                             // End for k
   t2 = time(NULL);

   fprintf(stderr, "\nDone!\n");

   dataOutput(rgbIm, sx, &output_name[0]);

   fprintf(stderr, "Total number of rays created: %ld\n", NUM_RAYS);
   fprintf(stderr, "Rays per second: %f\n", (double)NUM_RAYS / (double)difftime(t2, t1));

   // Exit section. Clean up and return.
   cleanup(object_list, texture_list); // Object and texture lists
   deleteImage(im);                    // Rendered image
   free(cam);                          // camera view
   free(wght);
   free(light_list);
   exit(0);
}
