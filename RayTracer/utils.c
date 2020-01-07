/*
   utils.c
   Basecode by F. J. Estrada
  
   Uses Tom F. El-Maraghi's code for computing inverse
   matrices.

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

#include "utils.h"
#include "affineTransforms.h"
// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4] = {{1.0, 0.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0, 0.0},
                       {0.0, 0.0, 1.0, 0.0},
                       {0.0, 0.0, 0.0, 1.0}};

#include "structures.c"

/////////////////////////////////////////////
// Primitive data structure section
/////////////////////////////////////////////
struct point *newPoint(double px, double py, double pz) {
   // Allocate a new point structure, initialize it to
   // the specified coordinates, and return a pointer
   // to it.

   struct point *pt = (struct point *)calloc(1, sizeof(struct point));
   if (!pt)
      fprintf(stderr, "Out of memory allocating point structure!\n");
   else {
      pt->x = px;
      pt->y = py;
      pt->z = pz;
      pt->w = 1.0;
   }
   return (pt);
}

struct pointLS *newPLS(struct point *p0, double r, double g, double b) {
   // Allocate a new point light sourse structure. Initialize the light
   // source to the specified RGB colour
   // Note that this is a point light source in that it is a single point
   // in space, if you also want a uniform direction for light over the
   // scene (a so-called directional light) you need to place the
   // light source really far away.

   struct pointLS *ls = (struct pointLS *)calloc(1, sizeof(struct pointLS));
   if (!ls)
      fprintf(stderr, "Out of memory allocating light source!\n");
   else {
      memcpy(&ls->p0, p0, sizeof(struct point));  // Copy light source location

      ls->col.R = r;  // Store light source colour and
      ls->col.G = g;  // intensity
      ls->col.B = b;
   }
   return (ls);
}

/////////////////////////////////////////////
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj) {
   // Transforms a ray using the inverse transform for the specified object. This is so that we can
   // use the intersection test for the canonical object. Note that this has to be done carefully!

   //copy original ray to new ray
   memcpy(ray_transformed, ray_orig, sizeof(struct ray3D));

   //inverse tranform ray origin
   ray_transformed->p0 = obj->Tinv * ray_transformed->p0;

   //inverse tranform ray direction without translation
   ray_transformed->d.w = 0;
   ray_transformed->d = obj->Tinv * ray_transformed->d;
   ray_transformed->d.w = 1;
}

inline void normalTransform(struct point *n_orig, struct point *n_transformed, struct object3D *obj) {
   // Computes the normal at an affinely transformed point given the original normal and the
   // object's inverse transformation. From the notes:
   // n_transformed=A^-T*n normalized.

   //take the transpose of Tinv
   struct transform inversetranspose;
   for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
         inversetranspose.T[i][j] = obj->Tinv.T[j][i];
      }
   }

   *n_transformed = inversetranspose * n_orig;
   n_transformed->w = 1;
   normalize(n_transformed);
}

void rayReflect(struct ray3D *ray_orig, struct point *p, struct point *n, struct ray3D *ray_reflected) {
   //this function assumes n is unit length!

   //reflection starts at point of intersection
   memcpy(&(ray_reflected->p0), p, sizeof(struct point));

   //r=d−2(d⋅n)n
   double ddotn = dot(&(ray_orig->d), n);
   ray_reflected->d.x = ray_orig->d.x - 2 * ddotn * n->x;
   ray_reflected->d.y = ray_orig->d.y - 2 * ddotn * n->y;
   ray_reflected->d.z = ray_orig->d.z - 2 * ddotn * n->z;
   ray_reflected->d.w = 1;
   normalize(&ray_reflected->d);
}

/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////
void insertObject(struct object3D *o, struct object3D **list) {
   if (o == NULL)
      return;
   // Inserts an object into the object list.
   if (*(list) == NULL) {
      *(list) = o;
      (*(list))->next = NULL;
   } else {
      o->next = (*(list))->next;
      (*(list))->next = o;
   }
}

struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
   // Intialize a new plane with the specified parameters:
   // ra, rd, rs, rg - Albedos for the components of the Phong model
   // r, g, b, - Colour for this plane
   // alpha - Transparency, must be set to 1 unless you are doing refraction
   // r_index - Refraction index if you are doing refraction.
   // shiny - Exponent for the specular component of the Phong model
   //
   // The plane is defined by the following vertices (CCW)
   // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
   // With normal vector (0,0,1) (i.e. parallel to the XY plane)

   struct object3D *plane = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!plane)
      fprintf(stderr, "Unable to allocate new plane, out of memory!\n");
   else {
      plane->alb.ra = ra;
      plane->alb.rd = rd;
      plane->alb.rs = rs;
      plane->alb.rg = rg;
      plane->col.R = r;
      plane->col.G = g;
      plane->col.B = b;
      plane->alpha = alpha;
      plane->r_index = r_index;
      plane->shinyness = shiny;
      plane->intersect = &planeIntersect;
      plane->surfaceCoords = &planeCoordinates;
      plane->randomPoint = &planeSample;
      plane->texImg = NULL;
      plane->photonMap = NULL;
      plane->normalMap = NULL;
      memcpy(&plane->T.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&plane->Tinv.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      plane->textureMap = &texMap;
      plane->frontAndBack = 1;
      plane->photonMapped = 0;
      plane->normalMapped = 0;
      plane->isCSG = 0;
      plane->isLightSource = 0;
      plane->CSGnext = NULL;
      plane->next = NULL;
   }
   return (plane);
}

struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
   // Intialize a new sphere with the specified parameters:
   // ra, rd, rs, rg - Albedos for the components of the Phong model
   // r, g, b, - Colour for this plane
   // alpha - Transparency, must be set to 1 unless you are doing refraction
   // r_index - Refraction index if you are doing refraction.
   // shiny -Exponent for the specular component of the Phong model
   //
   // This is assumed to represent a unit sphere centered at the origin.
   //

   struct object3D *sphere = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!sphere)
      fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
   else {
      sphere->alb.ra = ra;
      sphere->alb.rd = rd;
      sphere->alb.rs = rs;
      sphere->alb.rg = rg;
      sphere->col.R = r;
      sphere->col.G = g;
      sphere->col.B = b;
      sphere->alpha = alpha;
      sphere->r_index = r_index;
      sphere->shinyness = shiny;
      sphere->intersect = &sphereIntersect;
      sphere->surfaceCoords = &sphereCoordinates;
      sphere->randomPoint = &sphereSample;
      sphere->texImg = NULL;
      sphere->photonMap = NULL;
      sphere->normalMap = NULL;
      memcpy(&sphere->T.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&sphere->Tinv.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      sphere->textureMap = &texMap;
      sphere->frontAndBack = 0;
      sphere->photonMapped = 0;
      sphere->normalMapped = 0;
      sphere->isCSG = 0;
      sphere->isLightSource = 0;
      sphere->CSGnext = NULL;
      sphere->next = NULL;
   }
   return (sphere);
}

struct object3D *newBox(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
   // Intialize a new sphere with the specified parameters:
   // ra, rd, rs, rg - Albedos for the components of the Phong model
   // r, g, b, - Colour for this plane
   // alpha - Transparency, must be set to 1 unless you are doing refraction
   // r_index - Refraction index if you are doing refraction.
   // shiny -Exponent for the specular component of the Phong model
   //
   // This is assumed to represent a unit sphere centered at the origin.
   // Canonical object is a cube with side length = 1 it is defined by 6 planes
   // the cube is centered around the origin

   struct object3D *box = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!box)
      fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
   else {
      box->alb.ra = ra;
      box->alb.rd = rd;
      box->alb.rs = rs;
      box->alb.rg = rg;
      box->col.R = r;
      box->col.G = g;
      box->col.B = b;
      box->alpha = alpha;
      box->r_index = r_index;
      box->shinyness = shiny;
      box->intersect = &boxIntersect;
      box->texImg = NULL;
      box->photonMap = NULL;
      box->normalMap = NULL;
      memcpy(&box->T.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&box->Tinv.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      box->textureMap = &texMap;
      box->frontAndBack = 1;
      box->photonMapped = 0;
      box->normalMapped = 0;
      box->isCSG = 0;
      box->isLightSource = 0;
      box->CSGnext = NULL;
      box->next = NULL;
   }
   return (box);
}

struct object3D *newCyl(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
   ///////////////////////////////////////////////////////////////////////////////////////
   // TO DO:
   //	Complete the code to create and initialize a new cylinder object.
   // The canonical cylinder will be defined as x^2 + y^2 = 1 for |z| <= 1
   ///////////////////////////////////////////////////////////////////////////////////////
   struct object3D *cylinder = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!cylinder)
      fprintf(stderr, "Unable to allocate new cylinder, out of memory!\n");
   else {
      cylinder->alb.ra = ra;
      cylinder->alb.rd = rd;
      cylinder->alb.rs = rs;
      cylinder->alb.rg = rg;
      cylinder->col.R = r;
      cylinder->col.G = g;
      cylinder->col.B = b;
      cylinder->alpha = alpha;
      cylinder->r_index = r_index;
      cylinder->shinyness = shiny;
      cylinder->intersect = &cylIntersect;
      cylinder->surfaceCoords = &cylCoordinates;
      cylinder->randomPoint = &cylSample;
      cylinder->texImg = NULL;
      cylinder->photonMap = NULL;
      cylinder->normalMap = NULL;
      memcpy(&cylinder->T.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&cylinder->Tinv.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      cylinder->textureMap = &texMap;
      cylinder->frontAndBack = 0;
      cylinder->photonMapped = 0;
      cylinder->normalMapped = 0;
      cylinder->isCSG = 0;
      cylinder->isLightSource = 0;
      cylinder->CSGnext = NULL;
      cylinder->next = NULL;
   }
   return (cylinder);
}

struct object3D *newCone(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny) {
   ///////////////////////////////////////////////////////////////////////////////////////
   // TO DO:
   //	Complete the code to create and initialize a new cylinder object.
   // The canonical cylinder will be defined as x^2 + y^2 = 1 for |z| <= 1
   ///////////////////////////////////////////////////////////////////////////////////////
   struct object3D *cone = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!cone)
      fprintf(stderr, "Unable to allocate new cone, out of memory!\n");
   else {
      cone->alb.ra = ra;
      cone->alb.rd = rd;
      cone->alb.rs = rs;
      cone->alb.rg = rg;
      cone->col.R = r;
      cone->col.G = g;
      cone->col.B = b;
      cone->alpha = alpha;
      cone->r_index = r_index;
      cone->shinyness = shiny;
      cone->intersect = &coneIntersect;
      cone->texImg = NULL;
      cone->photonMap = NULL;
      cone->normalMap = NULL;
      memcpy(&cone->T.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&cone->Tinv.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      cone->textureMap = &texMap;
      cone->frontAndBack = 0;
      cone->photonMapped = 0;
      cone->normalMapped = 0;
      cone->isCSG = 0;
      cone->isLightSource = 0;
      cone->CSGnext = NULL;
      cone->next = NULL;
   }
   return (cone);
}

void importMesh(struct triangle **tg_list, struct point **v) {
   FILE *file = fopen("teapot.obj", "r");
   char t;
   double x, y, z;
   struct point v2_v1, v3_v1;
   int num_vertices = 0, num_lines = 0;
   ;
   while (fscanf(file, "%c", &t) && t == 'v' && fscanf(file, " %lf %lf %lf\n", &x, &y, &z) != EOF) {
      printf("x %f, y %f, z %f\n", x, y, z);
      num_vertices++;
   }
   struct point *vertices = (struct point *)calloc(num_vertices, sizeof(struct point));
   rewind(file);
   num_vertices = 0;
   while (fscanf(file, "%c", &t) && t == 'v' && fscanf(file, " %lf %lf %lf\n", &x, &y, &z) != EOF) {
      (vertices + num_vertices)->x = x;
      (vertices + num_vertices)->y = y;
      (vertices + num_vertices)->z = z;
      num_vertices++;
      num_lines++;
   }
   int v1, v2, v3;
   struct triangle *tg = NULL;
   struct triangle *head = NULL;
   fseek(file, -1, SEEK_CUR);
   while (fscanf(file, "%c %d %d %d\n", &t, &v1, &v2, &v3) != EOF && t == 'f') {
      printf("v1 %d, v2 %d, v3 %d\n", v1, v2, v3);
      tg = (struct triangle *)calloc(1, sizeof(struct triangle));
      if (tg == NULL) {
         fprintf(stderr, "Unable to allocate new triangle, out of memory!\n");
      }
      tg->v1 = v1;
      tg->v2 = v2;
      tg->v3 = v3;
      v2_v1.x = (vertices + v2 - 1)->x - (vertices + v1 - 1)->x;
      v2_v1.y = (vertices + v2 - 1)->y - (vertices + v1 - 1)->y;
      v2_v1.z = (vertices + v2 - 1)->z - (vertices + v1 - 1)->z;

      v3_v1.x = (vertices + v3 - 1)->x - (vertices + v1 - 1)->x;
      v3_v1.y = (vertices + v3 - 1)->y - (vertices + v1 - 1)->y;
      v3_v1.z = (vertices + v3 - 1)->z - (vertices + v1 - 1)->z;

      tg->n.x = ((v2_v1.y * v3_v1.z) - (v3_v1.y * v2_v1.z));
      tg->n.y = ((v3_v1.x * v2_v1.z) - (v2_v1.x * v3_v1.z));
      tg->n.z = ((v2_v1.x * v3_v1.y) - (v3_v1.x * v2_v1.y));
      tg->n.w = 1;
      //printf("before px: %f py: %f pz: %f\n", tg->n.px, tg->n.py, tg->n.pz);
      normalize(&(tg->n));
      //printf("after px: %f py: %f pz: %f\n", tg->n.px, tg->n.py, tg->n.pz);

      if (head == NULL) {
         head = tg;
      } else {
         tg->next = head;
         head = tg;
      }
      num_lines++;
   }
   *v = vertices;
   *tg_list = head;
   printf("Number of lines in file: %d\n", num_lines);
   fclose(file);
}

struct object3D *newMesh(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny, struct triangle *triangles, struct point *vertices) {
   ///////////////////////////////////////////////////////////////////////////////////////
   // TO DO:
   //	Complete the code to create and initialize a new cylinder object.
   // The canonical cylinder will be defined as x^2 + y^2 = 1 for |z| <= 1
   ///////////////////////////////////////////////////////////////////////////////////////
   struct object3D *mesh = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!mesh)
      fprintf(stderr, "Unable to allocate new mesh, out of memory!\n");
   else {
      mesh->alb.ra = ra;
      mesh->alb.rd = rd;
      mesh->alb.rs = rs;
      mesh->alb.rg = rg;
      mesh->col.R = r;
      mesh->col.G = g;
      mesh->col.B = b;
      mesh->alpha = alpha;
      mesh->r_index = r_index;
      mesh->shinyness = shiny;
      mesh->intersect = &meshIntersect;
      mesh->texImg = NULL;
      mesh->photonMap = NULL;
      mesh->normalMap = NULL;
      memcpy(&mesh->T.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&mesh->Tinv.T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      mesh->textureMap = &texMap;
      mesh->frontAndBack = 0;
      mesh->photonMapped = 0;
      mesh->normalMapped = 0;
      mesh->isCSG = 0;
      mesh->isLightSource = 0;
      mesh->CSGnext = NULL;
      mesh->next = NULL;
      mesh->triangles = triangles;
      mesh->vertices = vertices;
   }
   return (mesh);
}

void meshIntersect(struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
   struct triangle *t = mesh->triangles;
   *lambda = -1;
   double l = -1;
   struct point p1, n1;
   while (t != NULL) {
      //fprintf(stderr, "v1 %d v2 %d v3 %d\n", t->v1, t->v2, t->v3);
      triangleIntersect(t, mesh, ray, &l, &p1, &n1, a, b);
      if (l > THR && (l < *lambda || *lambda == -1)) {
         *lambda = l;
         memcpy(n, &n1, sizeof(struct point));
         memcpy(p, &p1, sizeof(struct point));

         printf("lambda: %f px: %f py: %f pz: %f\n", *lambda, p->x, p->y, p->z);
      }
      t = t->next;
   }
}

void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical plane.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, plane);
   *lambda = -1;
   struct point norm;
   // normal of canonical plane
   norm.x = 0;
   norm.y = 0;
   norm.z = 1;
   struct point p1;

   double l;
   double d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0) {
      p1.x = -ray_transformed.p0.x;
      p1.y = -ray_transformed.p0.y;
      p1.z = -ray_transformed.p0.z;

      l = dot(&(p1), &norm) / d_dot_n;
      // Check if the intersection point is inside the plane
      rayPosition(&ray_transformed, l, p);
      double x = p->x;
      double y = p->y;
      if (fabs(p->z) < THR && fabs(p->x) <= 1 + THR && fabs(p->y) <= 1 + THR) {
         *lambda = l;
         rayPosition(ray, l, p);
         //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
         normalTransform(&norm, n, plane);
         //printf("nt: (%f, %f, %f)\n", n->px, n->py, n->pz);
      }
      //TODO: compute a and b
      if (plane->texImg != NULL) {
         *a = (x + 1) / 2;
         *b = (-y + 1) / 2;
      }
   }
}

void triangleIntersect(struct triangle *triangle, struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical plane.
   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, mesh);

   struct point v1, v2, v3;
   double M, betta, gamma;
   struct point norm;

   double A, B, C, D, E, F, G, H, I, J, K, L;
   v1 = *(mesh->vertices + triangle->v1 - 1);
   v2 = *(mesh->vertices + triangle->v2 - 1);
   v3 = *(mesh->vertices + triangle->v3 - 1);

   A = v1.x - v2.x;
   B = v1.y - v2.y;
   C = v1.z - v2.z;
   D = v1.x - v3.x;
   E = v1.y - v3.y;
   F = v1.z - v3.z;
   G = ray_transformed.d.x;
   H = ray_transformed.d.y;
   I = ray_transformed.d.z;
   J = v1.x - ray_transformed.p0.x;
   K = v1.y - ray_transformed.p0.y;
   L = v1.z - ray_transformed.p0.z;
   //fprintf(stderr, "A: %f B: %f C: %f D: %f E: %f F: %f \n", A, B, C, D, E, F);

   *lambda = -1;
   M = A * (E * I - H * F) + B * (G * F - D * I) + C * (D * H - E * G);

   if (abs(M) > THR) {
      gamma = I * (A * K - J * B) + H * (J * C - A * L) + G * (B * L - K * C);
      gamma /= M;

      if (gamma < 0 || gamma > 1) {
         return;
      }
      betta = J * (E * I - H * F) + K * (G * F - D * I) + L * (D * H - E * G);
      betta /= M;

      if (betta < 0 || (betta > 1 - gamma)) {
         return;
      }
      *lambda = -(F * (A * K - J * B) + E * (J * C - A * L) + D * (B * L - K * C));
      *lambda /= M;
      normalTransform(&(triangle->n), n, mesh);

      rayPosition(ray, *lambda, p);
   }
}

void solveQuadratic(struct ray3D *ray, double *l1, double *l2) {
   struct point a_sub_c;
   double A, B, C, D;
   a_sub_c = ray->p0;
   A = dot(&(ray->d), &(ray->d));
   B = dot(&a_sub_c, &(ray->d));
   C = dot(&a_sub_c, &a_sub_c) - 1;
   D = B * B - A * C;

   if (D < 0) {
      return;
   }

   if (D > 0) {
      *l1 = -B / A - sqrt(D) / A;
      *l2 = -B / A + sqrt(D) / A;
   }
}

void boxIntersect(struct object3D *box, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical sphere.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
   double l = -1, temp_a, temp_b;
   struct point temp_pt, temp_n;
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, box);
   *lambda = -1;
   //TODO: compute a and b

   struct object3D *plane;
   plane = newPlane(box->alb.ra, box->alb.rd, box->alb.rs, box->alb.rg, box->col.R, box->col.G, box->col.B, box->alpha, box->r_index, box->shinyness);
   // get the intersection between 6 planes and see which one is closest
   // xy plane at z = 0.5
   plane->T *= Tr(0, 0, 1);
   invert(&plane->T.T[0][0], &plane->Tinv.T[0][0]);
   planeIntersect(plane, &ray_transformed, &l, &temp_pt, &temp_n, &temp_a, &temp_b);
   if (l != -1 && l > THR) {
      *lambda = l;
      rayPosition(ray, l, p);
      //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
      normalTransform(&temp_n, n, box);
   }
   l = -1;
   plane->T *= RotX(PI / 2);
   invert(&plane->T.T[0][0], &plane->Tinv.T[0][0]);
   planeIntersect(plane, &ray_transformed, &l, &temp_pt, &temp_n, &temp_a, &temp_b);

   if (l != -1 && l > THR && (*lambda == -1 || l < *lambda)) {
      *lambda = l;
      rayPosition(ray, l, p);
      //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
      normalTransform(&temp_n, n, box);
   }
   l = -1;

   plane->T *= RotX(PI / 2);
   invert(&plane->T.T[0][0], &plane->Tinv.T[0][0]);
   planeIntersect(plane, &ray_transformed, &l, &temp_pt, &temp_n, &temp_a, &temp_b);
   if (l != -1 && l > THR && (*lambda == -1 || l < *lambda)) {
      *lambda = l;
      rayPosition(ray, l, p);
      //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
      normalTransform(&temp_n, n, box);
   }
   l = -1;

   plane->T *= RotX(PI / 2);
   invert(&plane->T.T[0][0], &plane->Tinv.T[0][0]);
   planeIntersect(plane, &ray_transformed, &l, &temp_pt, &temp_n, &temp_a, &temp_b);
   if (l != -1 && l > THR && (*lambda == -1 || l < *lambda)) {
      *lambda = l;
      rayPosition(ray, l, p);
      //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
      normalTransform(&temp_n, n, box);
   }
   l = -1;

   plane->T *= RotZ(PI / 2);
   invert(&plane->T.T[0][0], &plane->Tinv.T[0][0]);
   planeIntersect(plane, &ray_transformed, &l, &temp_pt, &temp_n, &temp_a, &temp_b);
   if (l != -1 && l > THR && (*lambda == -1 || l < *lambda)) {
      *lambda = l;
      rayPosition(ray, l, p);
      //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
      normalTransform(&temp_n, n, box);
   }
   l = -1;

   plane->T *= RotZ(PI);
   invert(&plane->T.T[0][0], &plane->Tinv.T[0][0]);
   planeIntersect(plane, &ray_transformed, &l, &temp_pt, &temp_n, &temp_a, &temp_b);
   if (l != -1 && l > THR && (*lambda == -1 || l < *lambda)) {
      *lambda = l;
      rayPosition(ray, l, p);
      //printf("n: (%f, %f, %f)\n", n->px, n->py, n->pz);
      normalTransform(&temp_n, n, box);
   }
   free(plane);
}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical sphere.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
   struct ray3D ray_transformed;
   double l1 = -1, l2 = -1;
   rayTransform(ray, &ray_transformed, sphere);
   solveQuadratic(&ray_transformed, &l1, &l2);

   *lambda = -1;
   if (MAX(l1, l2) >= THR) {
      if (MIN(l1, l2) <= THR) {
         *lambda = MAX(l1, l2);
      } else {
         *lambda = MIN(l1, l2);
      }
   }
#ifdef DEBUG
   if (0) {
      printf("l1: %f, l2: %f\n", l1, l2);

      printf("sphere tray d: %f %f %f\n", ray_transformed.d.px, ray_transformed.d.py, ray_transformed.d.pz);
      printf("sphere tray p0: %f %f %f\n", ray_transformed.p0.px, ray_transformed.p0.py, ray_transformed.p0.pz);

      printf("sphere ray d: %f %f %f\n", ray->d.px, ray->d.py, ray->d.pz);
      printf("sphere ray p0: %f %f %f\n", ray->p0.px, ray->p0.py, ray->p0.pz);
   }
#endif  // DEBUG

   double x, y, z;
   if (*lambda != -1) {
      rayPosition(&ray_transformed, *lambda, p);
      n->x = p->x;
      n->y = p->y;
      n->z = p->z;
      x = p->x;
      y = p->y;
      z = p->z;
      n->w = 1;

      normalize(n);

      normalTransform(n, n, sphere);

      rayPosition(ray, *lambda, p);
   }
   //TODO: compute a and b
   *a = 0.5 + (atan2(z, x)) / (2 * PI);
   *b = 0.5 - (asin(y)) / (PI);
}

void cylIntersect(struct object3D *cylinder, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical cylinder.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
   struct ray3D ray_transformed;
   struct ray3D ray_copy;
   double l1 = -1, l2 = -1;
   rayTransform(ray, &ray_transformed, cylinder);
   *lambda = -1;
   rayTransform(ray, &ray_copy, cylinder);
   ray_copy.d.z = 0;
   ray_copy.p0.z = 0;
   solveQuadratic(&ray_copy, &l1, &l2);
   // Check if the z component of the intersection point falls within the range
   if (l1 > THR && fabs(l1 * ray_transformed.d.z + ray_transformed.p0.z) <= 1) {
      *lambda = l1;
      rayPosition(&ray_transformed, *lambda, p);
      n->x = 2 * p->x;
      n->y = 2 * p->y;
      n->z = 0;
   }

   // Check if the z component of the intersection point falls within the range
   if (l2 > THR && fabs(l2 * ray_transformed.d.z + ray_transformed.p0.z) <= 1 && (l2 < l1 || *lambda == -1)) {
      *lambda = l2;
      rayPosition(&ray_transformed, *lambda, p);
      n->x = 2 * p->x;
      n->y = 2 * p->y;
      n->z = 0;
   }

   struct point p1, cap_n, cap_p;
   cap_n.x = 0;
   cap_n.y = 0;
   cap_n.z = -1;

   double l;
   double d_dot_n;
   d_dot_n = dot(&(ray_transformed.d), &cap_n);
   if (d_dot_n != 0) {
      p1.x = -ray_transformed.p0.x;
      p1.y = -ray_transformed.p0.y;
      p1.z = -1 - ray_transformed.p0.z;

      l = dot(&(p1), &cap_n) / d_dot_n;

      rayPosition(&ray_transformed, l, &cap_p);
      // check if the cap intersection is within the cylinder boundaries
#ifdef DEBUG
      printf("l: %f, lambda: %f, x^2 + y^2: %f\n", l, *lambda, cap_p.px * cap_p.px + cap_p.py * cap_p.py);
#endif

      if (l > THR && cap_p.x * cap_p.x + cap_p.y * cap_p.y <= 1) {
         if ((*lambda != -1 && l < *lambda) || *lambda == -1) {
            *lambda = l;
            memcpy(n, &cap_n, sizeof(struct point));
         }
      }
   }

   cap_n.x = 0;
   cap_n.y = 0;
   cap_n.z = 1;

   d_dot_n = dot(&(ray_transformed.d), &cap_n);
   if (d_dot_n != 0) {
      p1.x = -ray_transformed.p0.x;
      p1.y = -ray_transformed.p0.y;
      p1.z = 1 - ray_transformed.p0.z;

      l = dot(&(p1), &cap_n) / d_dot_n;
      rayPosition(&ray_transformed, l, &cap_p);
#ifdef DEBUG
      printf("l: %f, lambda: %f, x^2 + y^2: %f\n", l, *lambda, cap_p.px * cap_p.px + cap_p.py * cap_p.py);
#endif

      if (l > THR && cap_p.x * cap_p.x + cap_p.y * cap_p.y <= 1) {
         if ((*lambda != -1 && l < *lambda) || *lambda == -1) {
            *lambda = l;
            memcpy(n, &cap_n, sizeof(struct point));
         }
      }
   }

#ifdef DEBUG
   printf("lambda: %f\n", *lambda);
#endif

   if (*lambda != -1) {
      normalTransform(n, n, cylinder);
      rayPosition(ray, *lambda, p);
   }
}

void coneIntersect(struct object3D *cone, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b) {
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical cylinder.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
   struct ray3D ray_transformed;
   struct ray3D ray_copy;
   double l1 = -1, l2 = -1;
   rayTransform(ray, &ray_transformed, cone);
   *lambda = -1;
   rayTransform(ray, &ray_copy, cone);
   ray_copy.d.z = 0;
   ray_copy.p0.z = 0;
   double A, B, C, D;
   A = ray_transformed.d.x * ray_transformed.d.x + ray_transformed.d.y * ray_transformed.d.y - ray_transformed.d.z * ray_transformed.d.z;
   B = ray_transformed.d.x * ray_transformed.p0.x + ray_transformed.d.y * ray_transformed.p0.y - ray_transformed.d.z * ray_transformed.p0.z;
   C = ray_transformed.p0.x * ray_transformed.p0.x + ray_transformed.p0.y * ray_transformed.p0.y - ray_transformed.p0.z * ray_transformed.p0.z;
   D = B * B - A * C;

   if (D < 0) {
      return;
   }

   if (D > 0) {
      l1 = -B / A - sqrt(D) / A;
      l2 = -B / A + sqrt(D) / A;
   }

   // Check if the z component of the intersection point falls within the range
   if (l1 > THR && l1 * ray_transformed.d.z + ray_transformed.p0.z <= 1 && l1 * ray_transformed.d.z + ray_transformed.p0.z >= 0) {
      *lambda = l1;
      rayPosition(&ray_transformed, *lambda, p);
      n->x = p->x;
      n->y = p->y;
      n->z = -p->z;
   }

   // Check if the z component of the intersection point falls within the range
   if (l2 > THR && (l2 * ray_transformed.d.z + ray_transformed.p0.z) <= 1 && (l2 * ray_transformed.d.z + ray_transformed.p0.z) >= 0 && (l2 < l1 || *lambda == -1)) {
      *lambda = l2;
      rayPosition(&ray_transformed, *lambda, p);
      n->x = p->x;
      n->y = p->y;
      n->z = -p->z;
   }

   struct point p1, cap_n, cap_p;
   cap_n.x = 0;
   cap_n.y = 0;
   cap_n.z = -1;

   double l;
   double d_dot_n;

   cap_n.x = 0;
   cap_n.y = 0;
   cap_n.z = 1;

   d_dot_n = dot(&(ray_transformed.d), &cap_n);
   if (d_dot_n != 0) {
      p1.x = -ray_transformed.p0.x;
      p1.y = -ray_transformed.p0.y;
      p1.z = 1 - ray_transformed.p0.z;

      l = dot(&(p1), &cap_n) / d_dot_n;
      rayPosition(&ray_transformed, l, &cap_p);
#ifdef DEBUG
      printf("l: %f, lambda: %f, x^2 + y^2: %f\n", l, *lambda, cap_p.px * cap_p.px + cap_p.py * cap_p.py);
#endif

      if (l > THR && cap_p.x * cap_p.x + cap_p.y * cap_p.y <= cap_p.z * cap_p.z) {
         if ((*lambda != -1 && l < *lambda) || *lambda == -1) {
            *lambda = l;
            memcpy(n, &cap_n, sizeof(struct point));
         }
      }
   }

#ifdef DEBUG
   printf("lambda: %f\n", *lambda);
#endif

   if (*lambda != -1) {
      normalTransform(n, n, cone);
      rayPosition(ray, *lambda, p);
   }
}

/////////////////////////////////////////////////////////////////
// Surface coordinates & random sampling on object surfaces
/////////////////////////////////////////////////////////////////
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z) {
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
   // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
   // bottom of the plane.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
}

void sphereCoordinates(struct object3D *sphere, double a, double b, double *x, double *y, double *z) {
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
   // 'a' in [0, 2*PI] corresponds to the spherical coordinate theta
   // 'b' in [-PI/2, PI/2] corresponds to the spherical coordinate phi

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
}

void cylCoordinates(struct object3D *cyl, double a, double b, double *x, double *y, double *z) {
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
   // 'a' in [0, 2*PI] corresponds to angle theta around the cylinder
   // 'b' in [0, 1] corresponds to height from the bottom

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
}

void planeSample(struct object3D *plane, double *x, double *y, double *z) {
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
   // Sapling should be uniform, meaning there should be an equal change of gedtting
   // any spot on the plane

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
}

void sphereSample(struct object3D *sphere, double *x, double *y, double *z) {
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the sphere
   // Sampling should be uniform - note that this is tricky for a sphere, do some
   // research and document in your report what method is used to do this, along
   // with a reference to your source.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
}

void cylSample(struct object3D *cyl, double *x, double *y, double *z) {
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the cylinder
   // Sampling should be uniform over the cylinder.

   /////////////////////////////////
   // TO DO: Complete this function.
   /////////////////////////////////
}

/////////////////////////////////
// Texture mapping functions
/////////////////////////////////
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list) {
   // Load a texture or normal map image from file and assign it to the
   // specified object.
   // type:   1  ->  Texture map  (RGB, .ppm)
   //         2  ->  Normal map   (RGB, .ppm)
   //         3  ->  Alpha map    (grayscale, .pgm)
   // Stores loaded images in a linked list to avoid replication
   struct image *im;
   struct textureNode *p;

   if (o != NULL) {
      // Check current linked list
      p = *(t_list);
      while (p != NULL) {
         if (strcmp(&p->name[0], filename) == 0) {
            // Found image already on the list
            if (type == 1)
               o->texImg = p->im;
            else if (type == 2)
               o->normalMap = p->im;
            else
               o->alphaMap = p->im;
            return;
         }
         p = p->next;
      }

      // Load this texture image
      if (type == 1 || type == 2)
         im = readPPMimage(filename);
      else if (type == 3)
         im = readPGMimage(filename);

      // Insert it into the texture list
      if (im != NULL) {
         p = (struct textureNode *)calloc(1, sizeof(struct textureNode));
         strcpy(&p->name[0], filename);
         p->type = type;
         p->im = im;
         p->next = NULL;
         // Insert into linked list
         if ((*(t_list)) == NULL)
            *(t_list) = p;
         else {
            p->next = (*(t_list))->next;
            (*(t_list))->next = p;
         }
         // Assign to object
         if (type == 1)
            o->texImg = im;
         else if (type == 2)
            o->normalMap = im;
         else
            o->alphaMap = im;
      }

   }  // end if (o != NULL)
}

void texMap(struct image *img, double a, double b, double *R, double *G, double *B) {
   /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
   a given object.

  The colour is returned in R, G, B. Uses bi-linear interpolation
  to determine texture colour.
 */

   a = MIN(1, MAX(0, a));
   b = MIN(1, MAX(0, b));
   a = a * (img->sx - 1);
   b = b * (img->sy - 1);

   int x1 = (int)floor(a);
   int y1 = (int)floor(b);
   int x2 = MIN(img->sy - 1, (int)ceil(a));
   int y2 = MIN(img->sy - 1, (int)ceil(b));
   //printf("a: %f b: %f\n", a, b);
   //printf("x: %d y: %d\n", x, y);
   double *rgbIm = (double *)img->rgbdata;
   *(R) = ((double)rgbIm[3 * (y1 * img->sx + x1) + 0]);
   *(G) = ((double)rgbIm[3 * (y1 * img->sx + x1) + 1]);
   *(B) = ((double)rgbIm[3 * (y1 * img->sx + x1) + 2]);
   //printf("r: %f g: %f b: %f\n", *R, *G, *B);
   //here
   return;
}

void texMapN(struct image *img, double a, double b, double *R, double *G, double *B) {
   /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
  a given object.

  The colour is returned in R, G, B. Uses nearest neighbour
  to determine texture colour.
 */

   a = MIN(1, MAX(0, a));
   b = MIN(1, MAX(0, b));
   //printf("a: %f b: %f\n", a, b);
   int x = (int)floor(a * (img->sx - 1));
   int y = (int)floor(b * (img->sy - 1));
   //printf("x: %d y: %d\n", x, y);
   double *rgbIm = (double *)img->rgbdata;
   *(R) = ((double)rgbIm[3 * (y * img->sx + x)]);
   *(G) = ((double)rgbIm[3 * (y * img->sx + x) + 1]);
   *(B) = ((double)rgbIm[3 * (y * img->sx + x) + 2]);
   //printf("r: %f g: %f b: %f\n", *R, *G, *B);
   //here
   return;
}

void alphaMap(struct image *img, double a, double b, double *alpha) {
   // Just like texture map but returns the alpha value at a,b,
   // notice that alpha maps are single layer grayscale images, hence
   // the separate function.

   *(alpha) = 1;  // Returns 1 which means fully opaque. Replace
   return;        // with your code if implementing alpha maps.
}

/////////////////////////////
// Light sources
/////////////////////////////
void insertPLS(struct pointLS *l, struct pointLS **list) {
   if (l == NULL)
      return;
   // Inserts a light source into the list of light sources
   if (*(list) == NULL) {
      *(list) = l;
      (*(list))->next = NULL;
   } else {
      l->next = (*(list))->next;
      (*(list))->next = l;
   }
}

void addAreaLight(double sx, double sy, double nx, double ny, double nz,
                  double tx, double ty, double tz, int N,
                  double r, double g, double b, struct object3D **o_list, struct pointLS **l_list) {
   /*
   This function sets up and inserts a rectangular area light source
   with size (sx, sy)
   orientation given by the normal vector (nx, ny, nz)
   centered at (tx, ty, tz)
   consisting of (N) point light sources (uniformly sampled)
   and with colour/intensity (r,g,b)

   Note that the light source must be visible as a uniformly colored rectangle which
   casts no shadows. If you require a lightsource to shade another, you must
   make it into a proper solid box with a back and sides of non-light-emitting
   material
 */

   // NOTE: The best way to implement area light sources is to random sample from the
   //       light source's object surface within rtShade(). This is a bit more tricky
   //       but reduces artifacts significantly. If you do that, then there is no need
   //       to insert a series of point lightsources in this function.
}

///////////////////////////////////
// Geometric transformation section
///////////////////////////////////

void invert(double *T, double *Tinv) {
   // Computes the inverse of transformation matrix T.
   // the result is returned in Tinv.

   double *U, *s, *V, *rv1;
   int singFlag, i;

   // Invert the affine transform
   U = NULL;
   s = NULL;
   V = NULL;
   rv1 = NULL;
   singFlag = 0;

   SVD(T, 4, 4, &U, &s, &V, &rv1);
   if (U == NULL || s == NULL || V == NULL) {
      fprintf(stderr, "Error: Matrix not invertible for this object, returning identity\n");
      memcpy(Tinv, eye4x4, 16 * sizeof(double));
      return;
   }

   // Check for singular matrices...
   for (i = 0; i < 4; i++)
      if (*(s + i) < 1e-9)
         singFlag = 1;
   if (singFlag) {
      fprintf(stderr, "Error: Transformation matrix is singular, returning identity\n");
      memcpy(Tinv, eye4x4, 16 * sizeof(double));
      return;
   }

   // Compute and store inverse matrix
   InvertMatrix(U, s, V, 4, Tinv);

   free(U);
   free(s);
   free(V);
}

void printmatrix(double mat[4][4]) {
   fprintf(stderr, "Matrix contains:\n");
   fprintf(stderr, "%f %f %f %f\n", mat[0][0], mat[0][1], mat[0][2], mat[0][3]);
   fprintf(stderr, "%f %f %f %f\n", mat[1][0], mat[1][1], mat[1][2], mat[1][3]);
   fprintf(stderr, "%f %f %f %f\n", mat[2][0], mat[2][1], mat[2][2], mat[2][3]);
   fprintf(stderr, "%f %f %f %f\n", mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
}

/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point *e, struct point *g, struct point *up, double f, double wl, double wt, double wsize) {
   /*
   This function sets up the camera axes and viewing direction as discussed in the
   lecture notes.
   e - Camera center
   g - Gaze direction
   up - Up vector
   fov - Field of view in degrees
   f - focal length
 */
   struct view *c;
   struct point *u, *v;

   u = v = NULL;

   // Allocate space for the camera structure
   c = (struct view *)calloc(1, sizeof(struct view));
   if (c == NULL) {
      fprintf(stderr, "Out of memory setting up camera model!\n");
      return (NULL);
   }

   // Set up camera center and axes
   c->e.x = e->x;  // Copy camera center location, note we must make sure
   c->e.y = e->y;  // the camera center provided to this function has pw=1
   c->e.z = e->z;
   c->e.w = 1;

   // Set up w vector (camera's Z axis). w=-g/||g||
   c->w.x = -g->x;
   c->w.y = -g->y;
   c->w.z = -g->z;
   c->w.w = 1;
   normalize(&c->w);

   // Set up the horizontal direction, which must be perpenticular to w and up
   u = cross(&c->w, up);
   normalize(u);
   c->u.x = u->x;
   c->u.y = u->y;
   c->u.z = u->z;
   c->u.w = 1;

   // Set up the remaining direction, v=(u x w)  - Mind the signs
   v = cross(&c->u, &c->w);
   normalize(v);
   c->v.x = v->x;
   c->v.y = v->y;
   c->v.z = v->z;
   c->v.w = 1;

   // Copy focal length and window size parameters
   c->f = f;
   c->wl = wl;
   c->wt = wt;
   c->wsize = wsize;

   // Set up coordinate conversion matrices
   // Camera2World matrix (M_cw in the notes)
   // Mind the indexing convention [row][col]
   c->C2W.T[0][0] = c->u.x;
   c->C2W.T[1][0] = c->u.y;
   c->C2W.T[2][0] = c->u.z;
   c->C2W.T[3][0] = 0;

   c->C2W.T[0][1] = c->v.x;
   c->C2W.T[1][1] = c->v.y;
   c->C2W.T[2][1] = c->v.z;
   c->C2W.T[3][1] = 0;

   c->C2W.T[0][2] = c->w.x;
   c->C2W.T[1][2] = c->w.y;
   c->C2W.T[2][2] = c->w.z;
   c->C2W.T[3][2] = 0;

   c->C2W.T[0][3] = c->e.x;
   c->C2W.T[1][3] = c->e.y;
   c->C2W.T[2][3] = c->e.z;
   c->C2W.T[3][3] = 1;

   // World2Camera matrix (M_wc in the notes)
   // Mind the indexing convention [row][col]
   c->W2C.T[0][0] = c->u.x;
   c->W2C.T[1][0] = c->v.x;
   c->W2C.T[2][0] = c->w.x;
   c->W2C.T[3][0] = 0;

   c->W2C.T[0][1] = c->u.y;
   c->W2C.T[1][1] = c->v.y;
   c->W2C.T[2][1] = c->w.y;
   c->W2C.T[3][1] = 0;

   c->W2C.T[0][2] = c->u.z;
   c->W2C.T[1][2] = c->v.z;
   c->W2C.T[2][2] = c->w.z;
   c->W2C.T[3][2] = 0;

   c->W2C.T[0][3] = -dot(&c->u, &c->e);
   c->W2C.T[1][3] = -dot(&c->v, &c->e);
   c->W2C.T[2][3] = -dot(&c->w, &c->e);
   c->W2C.T[3][3] = 1;

   free(u);
   free(v);
   return (c);
}

/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename) {
   // Reads an image from a .ppm file. A .ppm file is a very simple image representation
   // format with a text header followed by the binary RGB data at 24bits per pixel.
   // The header has the following form:
   //
   // P6
   // # One or more comment lines preceded by '#'
   // 340 200
   // 255
   //
   // The first line 'P6' is the .ppm format identifier, this is followed by one or more
   // lines with comments, typically used to inidicate which program generated the
   // .ppm file.
   // After the comments, a line with two integer values specifies the image resolution
   // as number of pixels in x and number of pixels in y.
   // The final line of the header stores the maximum value for pixels in the image,
   // usually 255.
   // After this last header line, binary data stores the RGB values for each pixel
   // in row-major order. Each pixel requires 3 bytes ordered R, G, and B.
   //
   // NOTE: Windows file handling is rather crotchetty. You may have to change the
   //       way this file is accessed if the images are being corrupted on read
   //       on Windows.
   //
   // readPPMdata converts the image colour information to floating point. This is so that
   // the texture mapping function doesn't have to do the conversion every time
   // it is asked to return the colour at a specific location.
   //

   FILE *f;
   struct image *im;
   char line[1024];
   int sizx, sizy;
   int i;
   unsigned char *tmp;
   double *fRGB;
   int tmpi;
   char *tmpc;

   im = (struct image *)calloc(1, sizeof(struct image));
   if (im != NULL) {
      im->rgbdata = NULL;
      f = fopen(filename, "rb+");
      if (f == NULL) {
         fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
         free(im);
         return (NULL);
      }
      tmpc = fgets(&line[0], 1000, f);
      if (strcmp(&line[0], "P6\n") != 0) {
         fprintf(stderr, "Wrong file format, not a .ppm file or header end-of-line characters missing\n");
         free(im);
         fclose(f);
         return (NULL);
      }
      // fprintf(stderr,"%s\n",line);
      // Skip over comments
      tmpc = fgets(&line[0], 511, f);
      while (line[0] == '#') {
         //fprintf(stderr,"%s",line);
         tmpc = fgets(&line[0], 511, f);
      }
      sscanf(&line[0], "%d %d\n", &sizx, &sizy);  // Read file size
      //fprintf(stderr,"nx=%d, ny=%d\n\n",sizx,sizy);
      im->sx = sizx;
      im->sy = sizy;

      tmpc = fgets(&line[0], 9, f);  // Read the remaining header line
      //fprintf(stderr,"%s\n",line);
      tmp = (unsigned char *)calloc(sizx * sizy * 3, sizeof(unsigned char));
      fRGB = (double *)calloc(sizx * sizy * 3, sizeof(double));
      if (tmp == NULL || fRGB == NULL) {
         fprintf(stderr, "Out of memory allocating space for image\n");
         free(im);
         fclose(f);
         return (NULL);
      }

      tmpi = fread(tmp, sizx * sizy * 3 * sizeof(unsigned char), 1, f);
      fclose(f);

      // Conversion to floating point
      for (i = 0; i < sizx * sizy * 3; i++)
         *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
      free(tmp);
      im->rgbdata = (void *)fRGB;

      return (im);
   }

   fprintf(stderr, "Unable to allocate memory for image structure\n");
   return (NULL);
}

struct image *readPGMimage(const char *filename) {
   // Just like readPPMimage() except it is used to load grayscale alpha maps. In
   // alpha maps, a value of 255 corresponds to alpha=1 (fully opaque) and 0
   // correspondst to alpha=0 (fully transparent).
   // A .pgm header of the following form is expected:
   //
   // P5
   // # One or more comment lines preceded by '#'
   // 340 200
   // 255
   //
   // readPGMdata converts the image grayscale data to double floating point in [0,1].

   FILE *f;
   struct image *im;
   char line[1024];
   int sizx, sizy;
   int i;
   unsigned char *tmp;
   double *fRGB;
   int tmpi;
   char *tmpc;

   im = (struct image *)calloc(1, sizeof(struct image));
   if (im != NULL) {
      im->rgbdata = NULL;
      f = fopen(filename, "rb+");
      if (f == NULL) {
         fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
         free(im);
         return (NULL);
      }
      tmpc = fgets(&line[0], 1000, f);
      if (strcmp(&line[0], "P5\n") != 0) {
         fprintf(stderr, "Wrong file format, not a .pgm file or header end-of-line characters missing\n");
         free(im);
         fclose(f);
         return (NULL);
      }
      // Skip over comments
      tmpc = fgets(&line[0], 511, f);
      while (line[0] == '#')
         tmpc = fgets(&line[0], 511, f);
      sscanf(&line[0], "%d %d\n", &sizx, &sizy);  // Read file size
      im->sx = sizx;
      im->sy = sizy;

      tmpc = fgets(&line[0], 9, f);  // Read the remaining header line
      tmp = (unsigned char *)calloc(sizx * sizy, sizeof(unsigned char));
      fRGB = (double *)calloc(sizx * sizy, sizeof(double));
      if (tmp == NULL || fRGB == NULL) {
         fprintf(stderr, "Out of memory allocating space for image\n");
         free(im);
         fclose(f);
         return (NULL);
      }

      tmpi = fread(tmp, sizx * sizy * sizeof(unsigned char), 1, f);
      fclose(f);

      // Conversion to double floating point
      for (i = 0; i < sizx * sizy; i++)
         *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
      free(tmp);
      im->rgbdata = (void *)fRGB;

      return (im);
   }

   fprintf(stderr, "Unable to allocate memory for image structure\n");
   return (NULL);
}

struct image *newImage(int size_x, int size_y) {
   // Allocates and returns a new image with all zeros. Assumes 24 bit per pixel,
   // unsigned char array.
   struct image *im;

   im = (struct image *)calloc(1, sizeof(struct image));
   if (im != NULL) {
      im->rgbdata = NULL;
      im->sx = size_x;
      im->sy = size_y;
      im->rgbdata = (void *)calloc(size_x * size_y * 3, sizeof(unsigned char));
      if (im->rgbdata != NULL)
         return (im);
   }
   fprintf(stderr, "Unable to allocate memory for new image\n");
   return (NULL);
}

void imageOutput(struct image *im, const char *filename) {
   // Writes out a .ppm file from the image data contained in 'im'.
   // Note that Windows typically doesn't know how to open .ppm
   // images. Use Gimp or any other seious image processing
   // software to display .ppm images.
   // Also, note that because of Windows file format management,
   // you may have to modify this file to get image output on
   // Windows machines to work properly.
   //
   // Assumes a 24 bit per pixel image stored as unsigned chars
   //

   FILE *f;

   if (im != NULL)
      if (im->rgbdata != NULL) {
         f = fopen(filename, "wb+");
         if (f == NULL) {
            fprintf(stderr, "Unable to open file %s for output! No image written\n", filename);
            return;
         }
         fprintf(f, "P6\n");
         fprintf(f, "# Output from RayTracer.c\n");
         fprintf(f, "%d %d\n", im->sx, im->sy);
         fprintf(f, "255\n");
         fwrite((unsigned char *)im->rgbdata, im->sx * im->sy * 3 * sizeof(unsigned char), 1, f);
         fclose(f);
         return;
      }
   fprintf(stderr, "imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im) {
   // De-allocates memory reserved for the image stored in 'im'
   if (im != NULL) {
      if (im->rgbdata != NULL)
         free(im->rgbdata);
      free(im);
   }
}

void cleanup(struct object3D *o_list, struct pointLS *l_list, struct textureNode *t_list) {
   // De-allocates memory reserved for the object list and the point light source
   // list. Note that *YOU* must de-allocate any memory reserved for images
   // rendered by the raytracer.
   struct object3D *p, *q;
   struct pointLS *r, *s;
   struct textureNode *t, *u;

   p = o_list;  // De-allocate all memory from objects in the list
   while (p != NULL) {
      q = p->next;
      if (p->photonMap != NULL)  // If object is photon mapped, free photon map memory
      {
         if (p->photonMap->rgbdata != NULL)
            free(p->photonMap->rgbdata);
         free(p->photonMap);
      }
      free(p);
      p = q;
   }

   r = l_list;  // Delete light source list
   while (r != NULL) {
      s = r->next;
      free(r);
      r = s;
   }

   t = t_list;  // Delete texture Images
   while (t != NULL) {
      u = t->next;
      if (t->im->rgbdata != NULL)
         free(t->im->rgbdata);
      free(t->im);
      free(t);
      t = u;
   }
}

struct point path(double lambda) {
   struct point p;
   double a = 2, b = 10, c = 2;
   p.x = (a * lambda) * cos(10 * PI * lambda);
   p.y = (b * lambda);
   p.z = (c * lambda) * sin(10 * PI * lambda);
   return p;
}