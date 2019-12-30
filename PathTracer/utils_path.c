/*
   PathTracer utilities.
   Basecode by F. J. Estrada
  
   Uses Tom F. El-Maraghi's code for computing inverse
   matrices.

   Written in part by Lioudmila Tishkina

   Silviu Draghici
*/

#include "utils_path.h"

// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4] = {{1.0, 0.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0, 0.0},
                       {0.0, 0.0, 1.0, 0.0},
                       {0.0, 0.0, 0.0, 1.0}};

/////////////////////////////////////////////
// Primitive data structure section
/////////////////////////////////////////////
struct point *newPoint(double px, double py, double pz)
{
   // Allocate a new point structure, initialize it to
   // the specified coordinates, and return a pointer
   // to it.

   struct point *pt = (struct point *)calloc(1, sizeof(struct point));
   if (!pt)
      fprintf(stderr, "Out of memory allocating point structure!\n");
   else
   {
      pt->x = px;
      pt->y = py;
      pt->z = pz;
      pt->w = 1.0;
   }
   return (pt);
}

void rayReflect(struct ray3D *ray_orig, struct point *p, struct point *n, struct ray3D *ray_reflected)
{
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
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj)
{
   // Transforms a ray using the inverse transform for the specified object. This is so that we can
   // use the intersection test for the canonical object. Note that this has to be done carefully!

   //copy original ray to new ray
   memcpy(ray_transformed, ray_orig, sizeof(struct ray3D));
   //apply negative translation

   //inverse tranform ray origin
   matVecMult(obj->Tinv, &(ray_transformed->p0));

   //inverse tranform ray direction without translation
   ray_transformed->d.w = 0;
   matVecMult(obj->Tinv, &(ray_transformed->d));
   ray_transformed->d.w = 1;
}

inline void normalTransform(struct point *n_orig, struct point *n_transformed, struct object3D *obj)
{
   // Computes the normal at an affinely transformed point given the original normal and the
   // object's inverse transformation. From the notes:
   // n_transformed=A^-T*n normalized.

   memcpy(n_transformed, n_orig, sizeof(struct point));
   n_transformed->w = 0;

   //take the transpose of Tinv without the transformation
   double A[4][4];
   for (int i = 0; i < 4; i++)
   {
      for (int j = 0; j < 4; j++)
      {
         A[i][j] = obj->Tinv[j][i];
      }
   }

   matVecMult(A, n_transformed);
   n_transformed->w = 1;
   normalize(n_transformed);
}

/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////
void insertObject(struct object3D *o, struct object3D **list)
{
   if (o == NULL)
      return;
   // Inserts an object into the object list.
   if (*(list) == NULL)
   {
      *(list) = o;
      (*(list))->next = NULL;
   }
   else
   {
      o->next = (*(list))->next;
      (*(list))->next = o;
   }
}

struct object3D *newPlane(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
   // Intialize a new plane with the specified parameters:
   // diffPct, reflPct, tranPct - specify the amount of diffuse, reflective, and
   //   refracting properties of the material. They *must* sum to 1.0
   // r, g, b, - Colour for this plane
   // refl_sig - Determines the amount of spread for reflection directions. If zero
   //   rays are reflected only along the perfect reflection direction, if non-zero,
   //   the perfect reflection direction is bent a bit (the amount is drawn from a
   //   zero-mean Gaussian distribution with sigma refl_sig). This makes the reflection
   //   component less sharp, and makes the material look more 'matte'
   // r_index - Refraction index for the refraction component.
   //
   // The plane is defined by the following vertices (CCW)
   // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
   // With normal vector (0,0,1) (i.e. parallel to the XY plane)

   struct object3D *plane = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!plane)
      fprintf(stderr, "Unable to allocate new plane, out of memory!\n");
   else
   {
      plane->diffPct = diffPct;
      plane->reflPct = reflPct;
      plane->tranPct = tranPct;
      plane->col.R = r;
      plane->col.G = g;
      plane->col.B = b;
      plane->refl_sig = refl_sig;
      plane->r_index = r_index;
      plane->intersect = &planeIntersect;
      plane->surfaceCoords = &planeCoordinates;
      plane->randomPoint = &planeSample;
      plane->texImg = NULL;
      plane->normalMap = NULL;
      plane->alphaMap = NULL;
      plane->intersectMap = NULL;
      memcpy(&plane->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&plane->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
      plane->textureMap = &texMapN;
      plane->frontAndBack = 1;
      plane->isCSG = 0;
      plane->isLightSource = 0;
      plane->LSweight = 1.0;
      plane->CSGnext = NULL;
      plane->next = NULL;
   }
   return (plane);
}

struct object3D *newSphere(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
   // Intialize a new sphere with the specified parameters. The parameters have the same meaning
   // as for planes, so have a look above at the newPlane() function comments.
   //
   // This is assumed to represent a unit sphere centered at the origin.

   struct object3D *sphere = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!sphere)
      fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
   else
   {
      sphere->diffPct = diffPct;
      sphere->reflPct = reflPct;
      sphere->tranPct = tranPct;
      sphere->col.R = r;
      sphere->col.G = g;
      sphere->col.B = b;
      sphere->refl_sig = refl_sig;
      sphere->r_index = r_index;
      sphere->intersect = &sphereIntersect;
      sphere->surfaceCoords = &sphereCoordinates;
      sphere->randomPoint = &sphereSample;
      sphere->texImg = NULL;
      sphere->normalMap = NULL;
      sphere->intersectMap = NULL;
      memcpy(&sphere->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&sphere->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
      sphere->textureMap = &texMap;
      sphere->frontAndBack = 0;
      sphere->isCSG = 0;
      sphere->isLightSource = 0;
      sphere->LSweight = 1.0;
      sphere->CSGnext = NULL;
      sphere->next = NULL;
   }
   return (sphere);
}

struct object3D *newCyl(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
   struct object3D *cylinder = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!cylinder)
      fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
   else
   {
      cylinder->diffPct = diffPct;
      cylinder->reflPct = reflPct;
      cylinder->tranPct = tranPct;
      cylinder->col.R = r;
      cylinder->col.G = g;
      cylinder->col.B = b;
      cylinder->refl_sig = refl_sig;
      cylinder->r_index = r_index;
      cylinder->intersect = &cylIntersect;
      cylinder->surfaceCoords = &cylCoordinates;
      cylinder->randomPoint = &cylSample;
      cylinder->texImg = NULL;
      cylinder->normalMap = NULL;
      cylinder->intersectMap = NULL;
      cylinder->alphaMap = NULL;
      memcpy(&cylinder->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&cylinder->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
      cylinder->textureMap = &texMap;
      cylinder->frontAndBack = 0;
      cylinder->isCSG = 0;
      cylinder->isLightSource = 0;
      cylinder->LSweight = 1.0;
      cylinder->CSGnext = NULL;
      cylinder->next = NULL;
   }
   return (cylinder);
}

struct object3D *newBox(double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
   // Intialize a new sphere with the specified parameters:
   // ra, rd, rs, rg - Albedos for the components of the Phong model
   // r, g, b, - Colour for this plane
   // alpha - Transparency, must be set to 1 unless you are doing refraction
   // r_index - Refraction index if you are doing refraction.
   // shiny -Exponent for the specular component of the Phong model
   //
   // Canonical box is a cube with sides of size 1
   // The cube defined by 2 points that define the extent of the box
   // max_pt (0.5, 0.5, 0.5) and min_pt (-0.5, -0.5, -0.5)

   struct object3D *box = (struct object3D *)calloc(1, sizeof(struct object3D));

   if (!box)
      fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
   else
   {
      box->diffPct = diffPct;
      box->reflPct = reflPct;
      box->tranPct = tranPct;
      box->col.R = r;
      box->col.G = g;
      box->col.B = b;
      box->refl_sig = refl_sig;
      box->r_index = r_index;
      box->intersect = &boxIntersect;
      box->texImg = NULL;
      box->normalMap = NULL;
      memcpy(&box->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&box->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
      box->textureMap = &texMap;
      box->frontAndBack = 0;
      box->isCSG = 0;
      box->isLightSource = 0;
      box->CSGnext = NULL;
      box->next = NULL;
   }
   return (box);
}

struct object3D *newMesh(char *file, double diffPct, double reflPct, double tranPct, double r, double g, double b, double refl_sig, double r_index)
{
   struct object3D *mesh = (struct object3D *)calloc(1, sizeof(struct object3D));
   struct triangle *triangles;
   struct point *vertices;
   struct point *normals;
   struct point *min_pt = (struct point *)calloc(1, sizeof(struct point));
   struct point *max_pt = (struct point *)calloc(1, sizeof(struct point));
   importMesh(file, &triangles, &vertices, &normals, max_pt, min_pt);

   if (!mesh)
   {
      fprintf(stderr, "Unable to allocate new mesh, out of memory!\n");
   }
   else
   {
      mesh->diffPct = diffPct;
      mesh->reflPct = reflPct;
      mesh->tranPct = tranPct;
      mesh->col.R = r;
      mesh->col.G = g;
      mesh->col.B = b;
      mesh->refl_sig = refl_sig;
      mesh->r_index = r_index;
      mesh->intersect = &meshIntersect;
      mesh->texImg = NULL;
      mesh->normalMap = NULL;
      mesh->intersectMap = NULL;
      mesh->alphaMap = NULL;
      memcpy(&mesh->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
      memcpy(&mesh->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
      mesh->textureMap = &texMap;
      mesh->frontAndBack = 0;
      mesh->isCSG = 0;
      mesh->isLightSource = 0;
      mesh->CSGnext = NULL;
      mesh->next = NULL;
      mesh->LSweight = 1.0;
      mesh->triangles = triangles;
      mesh->vertices = vertices;
      mesh->normals = normals;
      mesh->min = min_pt;
      mesh->max = max_pt;
   }
   return (mesh);
}

void importMesh(char *filename, struct triangle **tg_list, struct point **v, struct point **n, struct point *max_pt, struct point *min_pt)
{
   FILE *f = fopen(filename, "r");

   if (f == NULL)
   {
      fprintf(stderr, "Unable to open obj file %s\n", filename);
      return;
   }
   char t1, t2;
   double x, y, z;
   struct point v2_v1, v3_v1;
   int num_vertices = 0, num_lines = 0;
   // will store the extent of the box
   min_pt->w = 1;
   max_pt->w = 1;
   while (fscanf(f, "%c", &t1) && t1 == 'v' && fscanf(f, " %lf %lf %lf\n", &x, &y, &z) != EOF)
   {
      min_pt->x = x;
      min_pt->y = y;
      min_pt->z = z;
      max_pt->x = x;
      max_pt->y = y;
      max_pt->z = z;

      //printf("x %f, y %f, z %f\n", x, y, z);
      num_vertices++;
   }
   printf("min x %f y %f z %f\n", min_pt->x, min_pt->y, min_pt->pz;
   printf("max x %f y %f z %f\n", max_pt->x, max_pt->y, max_pt->pz;

   struct point *vertices = (struct point *)calloc(num_vertices, sizeof(struct point));
   struct point *normals = (struct point *)calloc(num_vertices, sizeof(struct point));
   rewind(f);
   num_vertices = 0;
   while (fscanf(f, "%c", &t1) && t1 == 'v' && fscanf(f, " %lf %lf %lf\n", &x, &y, &z) != EOF)
   {
      (vertices + num_vertices)->x = x;
      (vertices + num_vertices)->y = y;
      (vertices + num_vertices)->z = z;
      min_pt->x = MIN(x, min_pt->x);
      min_pt->y = MIN(y, min_pt->y);
      min_pt->z = MIN(z, min_pt->z);

      max_pt->x = MAX(x, max_pt->x);
      max_pt->y = MAX(y, max_pt->y);
      max_pt->z = MAX(z, max_pt->z);

      num_vertices++;
      num_lines++;
   }

   fseek(f, -4, SEEK_CUR);
   num_vertices = 0;
   while (fscanf(f, "%c%c", &t1, &t2) && t1 == 'v' && t2 == 'n' && fscanf(f, " %lf %lf %lf\n", &x, &y, &z) != EOF)
   {
      (normals + num_vertices)->x = x;
      (normals + num_vertices)->y = y;
      (normals + num_vertices)->z = z;
      num_vertices++;

      num_lines++;
   }

   // printf("min x %f y %f z %f\n", min_pt->px, min_pt->py, min_pt->pz);
   // printf("max x %f y %f z %f\n", max_pt->px, max_pt->py, max_pt->pz);
   int v1, v2, v3, vt1, vt2, vt3, n1, n2, n3;
   struct triangle *tg = NULL;
   struct triangle *head = NULL;
   fseek(f, -2, SEEK_CUR);
   while (fscanf(f, "%c %d %d %d\n", &t1, &v1, &v2, &v3) != EOF && t1 == 'f')
   {
      tg = (struct triangle *)calloc(1, sizeof(struct triangle));
      if (tg == NULL)
      {
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

      tg->n.x = ((v2_v1.y * v3_v1.pz - (v3_v1.y * v2_v1.pzz;
      tg->n.y = ((v3_v1.x * v2_v1.pz - (v2_v1.x * v3_v1.pz);
      tg->n.z = ((v2_v1.x * v3_v1.y) - (v3_v1.x * v2_v1.y));
      tg->n.w = 1;
      //printf("before px: %f py: %f pz: %f\n", tg->n.px, tg->n.py, tg->n.pz);
      normalize(&(tg->n));
      //printf("after px: %f py: %f pz: %f\n", tg->n.px, tg->n.py, tg->n.pz);

      if (head == NULL)
      {
         head = tg;
      }
      else
      {
         tg->next = head;
         head = tg;
      }
      num_lines++;
   }
   *v = vertices;
   *tg_list = head;
   *n = normals;
   printf("Number of lines in file: %d\n", num_lines);
   fclose(f);
}

void meshIntersect(struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b)
{
   // check if the bounding box is intersected before checking the mesh
   // struct ray3D ray_transformed;
   // rayTransform(ray, &ray_transformed, mesh);
   double temp_l, temp_a, temp_b;
   struct point temp_n, temp_p;
   *lambda = -1;

   boundingBoxIntersect(mesh, ray, &temp_l);
   if (temp_l < 0)
   {
      return;
   }
   struct triangle *t = mesh->triangles;
   double l = -1;
   struct point p1, n1;
   while (t != NULL)
   {
      //fprintf(stderr, "v1 %d v2 %d v3 %d\n", t->v1, t->v2, t->v3);
      triangleIntersect(t, mesh, ray, &l, &p1, &n1, a, b);
      if (l > THR && (l < *lambda || *lambda == -1))
      {
         *lambda = l;
         memcpy(n, &n1, sizeof(struct point));
         memcpy(p, &p1, sizeof(struct point));
      }
      t = t->next;
   }
}

void triangleIntersect(struct triangle *triangle, struct object3D *mesh, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b)
{
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical plane.
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

   if (abs(M) > THR)
   {
      gamma = I * (A * K - J * B) + H * (J * C - A * L) + G * (B * L - K * C);
      gamma /= M;

      if (gamma < 0 || gamma > 1)
      {
         return;
      }
      betta = J * (E * I - H * F) + K * (G * F - D * I) + L * (D * H - E * G);
      betta /= M;

      if (betta < 0 || (betta > 1 - gamma))
      {
         return;
      }
      *lambda = -(F * (A * K - J * B) + E * (J * C - A * L) + D * (B * L - K * C));
      *lambda /= M;

      rayPosition(ray, *lambda, p);
      struct point temp_p;
      rayPosition(&ray_transformed, *lambda, &temp_p);

      struct point n1, n2, n3;
      n1 = *(mesh->normals + triangle->v1 - 1);
      n2 = *(mesh->normals + triangle->v2 - 1);
      n3 = *(mesh->normals + triangle->v3 - 1);

      //Barycentric
      struct point v1_p, v2_p, v3_p, v1_v2, v1_v3;
      v1_p.x = v1.x - temp_p.x;
      v1_p.y = v1.y - temp_p.y;
      v1_p.z = v1.z - temp_p.z;
      v2_p.x = v2.x - temp_p.x;
      v2_p.y = v2.y - temp_p.y;
      v2_p.z = v2.z - temp_p.z;
      v3_p.x = v3.x - temp_p.x;
      v3_p.y = v3.y - temp_p.y;
      v3_p.z = v3.z - temp_p.z;
      v1_v2.x = v1.x - v2.x;
      v1_v2.y = v1.y - v2.y;
      v1_v2.z = v1.z - v2.z;
      v1_v3.x = v1.x - v3.x;
      v1_v3.y = v1.y - v3.y;
      v1_v3.z = v1.z - v3.z;

      double a1, a2, a3, abig;
      abig = fabs((v1_v2.y * v1_v3.pz- v1_v2.pz* v1_v3.y) - (v1_v2.x * v1_v3.pzz v1_v2.pzz v1_v3.x) + (v1_v2.x * v1_v3.y - v1_v2.y * v1_v3.x));

      a3 = fabs((v1_p.y * v2_p.pz- v1_p.pz* v2_p.y) - (v1_p.x * v2_p.pzz v1_p.pzz v2_p.x) + (v1_p.x * v2_p.y - v1_p.y * v2_p.x)) / abig;
      a2 = fabs((v3_p.y * v1_p.pz- v3_p.pz* v1_p.y) - (v3_p.x * v1_p.pzz v3_p.pzz v1_p.x) + (v3_p.x * v1_p.y - v3_p.y * v1_p.x)) / abig;
      a1 = fabs((v2_p.y * v3_p.pz- v2_p.pz* v3_p.y) - (v2_p.x * v3_p.pzz v2_p.pzz v3_p.x) + (v2_p.x * v3_p.y - v2_p.y * v3_p.x)) / abig;
      struct point temp_n;
      temp_n.x = a1 * (n1.x) + a2 * (n2.x) + a3 * (n3.x);
      temp_n.y = a1 * (n1.y) + a2 * (n2.y) + a3 * (n3.y);
      temp_n.z = a1 * (n1.z) + a2 * (n2.z) + a3 * (n3.z);
      normalize(&temp_n);
      normalTransform(&temp_n, n, mesh);
   }
}

void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b)
{
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical plane.

   struct ray3D ray_transformed;
   double n_delta_x = 0, n_delta_y = 0, n_delta_z = 0;

   rayTransform(ray, &ray_transformed, plane);
   *lambda = -1;
   struct point norm;
   // normal of canonical plane
   norm.x = 0;
   norm.y = 0;
   norm.z = 1;
   norm.w = 1;
   struct point p1;

   double l;
   double d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -ray_transformed.p0.x;
      p1.y = -ray_transformed.p0.y;
      p1.z = -ray_transformed.p0.z;
      p1.w = 1;

      l = dot(&(p1), &norm) / d_dot_n;
      // Check if the intersection point is inside the plane
      rayPosition(&ray_transformed, l, p);
      double x = p->x;
      double y = p->y;
      if (fabs(p->z) < THR && fabs(p->x) <= 1 + THR && fabs(p->y) <= 1 + THR)
      {
         *lambda = l;
         rayPosition(ray, l, p);
         *a = (x + 1) / 2;
         *b = (-y + 1) / 2;
         normalTransform(&norm, n, plane);
      }
   }
}

// checks intersection with bounding box using extent of the bounding box
// assumes a canonical position of box and object
void boundingBoxIntersect(struct object3D *mesh, struct ray3D *ray, double *lambda)
{
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, mesh);
   *lambda = -1;
   struct point norm;
   struct point temp_p;
   struct point p1;
   struct point e;
   double x;
   double y;
   // printf("ray px %f py %f pz %f dx %f dy %f dz %f \n", ray->d.px, ray->d.py, ray->d.pz, ray->p0.px, ray->p0.py, ray->p0.pz);

   double l = -1;
   double d_dot_n;
   e.x = ray_transformed.p0.x;
   e.y = ray_transformed.p0.y;
   e.z = ray_transformed.p0.z;
   double min_x = mesh->min->x;
   double min_y = mesh->min->y;
   double min_z = mesh->min->z;

   double max_x = mesh->max->x;
   double max_y = mesh->max->y;
   double max_z = mesh->max->z;

   norm.x = 1;
   norm.y = 0;
   norm.z = 0;
   norm.w = 1;
   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = max_x - e.x;
      p1.y = -e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.x - max_x) < THR && temp_p.z < max_z + THR && temp_p.z > min_z + THR && temp_p.y < max_y + THR && temp_p.y > min_y + THR && (*lambda == -1 || l < *lambda))
         {
            // printf("l1 ray px %f py %f pz %f dx %f dy %f dz %f \n", ray->d.px, ray->d.py, ray->d.pz, ray->p0.px, ray->p0.py, ray->p0.pz);
            *lambda = l;
            // printf("l1 %f px %f py %f pz %f\n", l, p->px, p->py, p->pz);
         }
      }
   }

   l = -1;
   norm.x = -1;
   norm.y = 0;
   norm.z = 0;
   norm.w = 1;
   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = min_x - e.x;
      p1.y = -e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l >= 0)
      {

         // Check if the intersection point is inside the plane

         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.x - min_x) < THR && temp_p.z < max_z + THR && temp_p.z > min_z + THR && temp_p.y < max_y + THR && temp_p.y > min_y + THR && (*lambda == -1 || l < *lambda))
         {
            // printf("l2 ray px %f py %f pz %f dx %f dy %f dz %f \n", ray->d.px, ray->d.py, ray->d.pz, ray->p0.px, ray->p0.py, ray->p0.pz);
            *lambda = l;
            // printf("l2 %f px %f py %f pz %f\n", l, p->px, p->py, p->pz);
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = 1;
   norm.z = 0;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = max_y - e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (temp_p.x < max_x + THR && temp_p.x > min_x + THR && temp_p.z < max_z + THR && temp_p.z > min_z + THR && fabs(temp_p.y - max_y) < THR && (*lambda == -1 || l < *lambda))
         {

            *lambda = l;
            //printf("l3 %f px %f py %f pz %f\n", l, p->px, p->py, p->pz);
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = -1;
   norm.z = 0;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = min_y - e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (temp_p.x < max_x + THR && temp_p.x > min_x + THR && temp_p.z < max_z + THR && temp_p.z > min_z + THR && fabs(temp_p.y - min_y) < THR && (*lambda == -1 || l < *lambda))
         {
            *lambda = l;
            //printf("l4 %f px %f py %f pz %f\n", l, p->px, p->py, p->pz);
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = 0;
   norm.z = 1;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = -e.y;
      p1.z = max_z - e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (temp_p.x < max_x + THR && temp_p.x > min_x + THR && temp_p.y < max_y + THR && temp_p.y > min_y + THR && fabs(temp_p.pzz max_z) < THR && (*lambda == -1 || l < *lambda))
         {
            *lambda = l;
            //printf("l5 %f px %f py %f pz %f\n", l, p->px, p->py, p->pz);
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = 0;
   norm.z = -1;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = -e.y;
      p1.z = min_z - e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (temp_p.x < max_x + THR && temp_p.x > min_x + THR && temp_p.y < max_y + THR && temp_p.y > min_y + THR && fabs(temp_p.pzz min_z) < THR && (*lambda == -1 || l < *lambda))
         {
            *lambda = l;
            //printf("l6 %f px %f py %f pz %f\n", l, p->px, p->py, p->pz);
         }
      }
   }

   // if (*lambda > -1)
   // {
   //   printf("lambda: %f\n", *lambda);
   //   printf("%f px %f py %f pz %f\n", *lambda, p->px, p->py, p->pz);
   //   printf("%f nx %f ny %f nz %f\n", *lambda, n->px, n->py, n->pz);

   //   printf("-----\n");
   // }
}

void boxIntersect(struct object3D *box, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b)
{
   struct ray3D ray_transformed;
   rayTransform(ray, &ray_transformed, box);
   *lambda = -1;
   struct point norm;
   struct point temp_p;
   struct point p1;
   struct point e;
   double x;
   double y;

   double l = -1;
   double d_dot_n;
   e.x = ray_transformed.p0.x;
   e.y = ray_transformed.p0.y;
   e.z = ray_transformed.p0.z;

   norm.x = 1;
   norm.y = 0;
   norm.z = 0;
   norm.w = 1;
   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = 0.5 - e.x;
      p1.y = -e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l > THR)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.x - 0.5) < THR && fabs(temp_p.z) < 0.5 + THR && fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
         {
            rayPosition(ray, l, p);
            *lambda = l;
            normalTransform(&norm, n, box);
         }
      }
   }

   l = -1;
   norm.x = -1;
   norm.y = 0;
   norm.z = 0;
   norm.w = 1;
   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -0.5 - e.x;
      p1.y = -e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l > THR)
      {
         // Check if the intersection point is inside the plane

         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.x + 0.5) < THR && fabs(temp_p.z) < 0.5 + THR && fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
         {
            rayPosition(ray, l, p);
            *lambda = l;
            normalTransform(&norm, n, box);
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = 1;
   norm.z = 0;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = 0.5 - e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;

      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.z) < 0.5 + THR && fabs(temp_p.x) < 0.5 + THR && fabs(temp_p.y - 0.5) < THR && (*lambda == -1 || l < *lambda))
         {
            rayPosition(ray, l, p);

            *lambda = l;
            normalTransform(&norm, n, box);
            x = p->x;
            y = p->y;
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = -1;
   norm.z = 0;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = -0.5 - e.y;
      p1.z = -e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.z) < 0.5 + THR && fabs(temp_p.x) < 0.5 + THR && fabs(temp_p.y + 0.5) < THR && (*lambda == -1 || l < *lambda))
         {
            rayPosition(ray, l, p);
            *lambda = l;
            normalTransform(&norm, n, box);
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = 0;
   norm.z = 1;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = -e.y;
      p1.z = 0.5 - e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.z - 0.5) < THR && fabs(temp_p.x) < 0.5 + THR && fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
         {
            rayPosition(ray, l, p);
            *lambda = l;
            normalTransform(&norm, n, box);
         }
      }
   }

   l = -1;
   norm.x = 0;
   norm.y = 0;
   norm.z = -1;
   norm.w = 1;

   d_dot_n = dot(&(ray_transformed.d), &norm);
   if (d_dot_n != 0)
   {
      p1.x = -e.x;
      p1.y = -e.y;
      p1.z = -0.5 - e.z;

      l = dot(&(p1), &norm) / d_dot_n;
      if (l >= 0)
      {
         // Check if the intersection point is inside the plane
         rayPosition(&ray_transformed, l, &temp_p);

         if (fabs(temp_p.z + 0.5) < THR && fabs(temp_p.x) < 0.5 + THR && fabs(temp_p.y) < 0.5 + THR && (*lambda == -1 || l < *lambda))
         {
            rayPosition(ray, l, p);
            *lambda = l;

            normalTransform(&norm, n, box);
         }
      }
   }

   if (*lambda > -1)
   {
      x = p->x;
      y = p->y;
      if (box->texImg != NULL)
      {
         *a = (x + 1) / 2;
         *b = (-y + 1) / 2;
      }
   }
}

void solveQuadratic(struct ray3D *ray, double *l1, double *l2)
{
   struct point a_sub_c;
   double A, B, C, D;
   a_sub_c = ray->p0;
   A = dot(&(ray->d), &(ray->d));
   B = dot(&a_sub_c, &(ray->d));
   C = dot(&a_sub_c, &a_sub_c) - 1;
   D = B * B - A * C;

   if (D < 0)
   {
      return;
   }

   if (D > 0)
   {
      *l1 = -B / A - sqrt(D) / A;
      *l2 = -B / A + sqrt(D) / A;
   }
}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point *p, struct point *n, double *a, double *b)
{
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical sphere.

   struct ray3D ray_transformed;
   double l1 = -1, l2 = -1;
   rayTransform(ray, &ray_transformed, sphere);
   solveQuadratic(&ray_transformed, &l1, &l2);

   *lambda = -1;
   if (MAX(l1, l2) >= THR)
   {
      if (MIN(l1, l2) <= THR)
      {
         *lambda = MAX(l1, l2);
      }
      else
      {
         *lambda = MIN(l1, l2);
      }
   }
#ifdef DEBUG
   if (0)
   {

      printf("l1: %f, l2: %f\n", l1, l2);

      printf("sphere tray d: %f %f %f\n", ray_transformed.d.px, ray_transformed.d.py, ray_transformed.d.pz);
      printf("sphere tray p0: %f %f %f\n", ray_transformed.p0.px, ray_transformed.p0.py, ray_transformed.p0.pz);

      printf("sphere ray d: %f %f %f\n", ray->d.px, ray->d.py, ray->d.pz);
      printf("sphere ray p0: %f %f %f\n", ray->p0.px, ray->p0.py, ray->p0.pz);
   }
#endif // DEBUG

   double x = 0, y = 0, z = 0;
   if (*lambda != -1)
   {

      rayPosition(&ray_transformed, *lambda, p);
      n->x = p->x;
      n->y = p->y;
      n->z = p->z;
      x = p->x;
      y = p->y;
      z = p->z;
      n->w = 1;

      //Source: https://en.wikipedia.org/wiki/UV_mapping
      *a = 0.5 + (atan2(z, x)) / (2 * PI);
      *b = 0.5 - (asin(y)) / (PI);

      normalize(n);

      normalTransform(n, n, sphere);

      rayPosition(ray, *lambda, p);
   }
}

void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda, struct point *p, struct point *n, double *a, double *b)
{
   // Computes and returns the value of 'lambda' at the intersection
   // between the specified ray and the specified canonical cylinder.

   struct ray3D ray_transformed;
   struct ray3D ray_copy;
   double l1 = -1, l2 = -1;
   rayTransform(r, &ray_transformed, cylinder);
   *lambda = -1;
   rayTransform(r, &ray_copy, cylinder);
   ray_copy.d.z = 0;
   ray_copy.p0.z = 0;
   (*a) = 0;
   (*b) = 0;
   solveQuadratic(&ray_copy, &l1, &l2);
   // Check if the z component of the intersection point falls within the range
   if (l1 > THR && fabs(l1 * ray_transformed.d.z + ray_transformed.p0.z) <= 1)
   {
      *lambda = l1;
      rayPosition(&ray_transformed, *lambda, p);
      *a = atan2(p->x, p->y) / (2 * PI);
      *b = (1 + p->z) / 2;
      n->x = 2 * p->x;
      n->y = 2 * p->y;
      n->z = 0;
   }

   // Check if the z component of the intersection point falls within the range
   if (l2 > THR && fabs(l2 * ray_transformed.d.z + ray_transformed.p0.z) <= 1 && (l2 < l1 || *lambda == -1))
   {
      *lambda = l2;
      rayPosition(&ray_transformed, *lambda, p);
      n->x = p->x;
      n->y = p->y;
      n->z = 0;
   }

   struct point p1, cap_n, cap_p;
   cap_n.x = 0;
   cap_n.y = 0;
   cap_n.z = -1;

   double l;
   double d_dot_n;
   d_dot_n = dot(&(ray_transformed.d), &cap_n);
   if (d_dot_n != 0)
   {
      p1.x = -ray_transformed.p0.x;
      p1.y = -ray_transformed.p0.y;
      p1.z = -1 - ray_transformed.p0.z;

      l = dot(&(p1), &cap_n) / d_dot_n;

      rayPosition(&ray_transformed, l, &cap_p);
      // check if the cap intersection is within the cylinder boundaries
#ifdef DEBUG
      printf("l: %f, lambda: %f, x^2 + y^2: %f\n", l, *lambda, cap_p.px * cap_p.px + cap_p.py * cap_p.py);
#endif

      if (l > THR && cap_p.x * cap_p.x + cap_p.y * cap_p.y <= 1)
      {
         if ((*lambda != -1 && l < *lambda) || *lambda == -1)
         {
            *lambda = l;
            memcpy(n, &cap_n, sizeof(struct point));
         }
      }
   }

   cap_n.x = 0;
   cap_n.y = 0;
   cap_n.z = 1;

   d_dot_n = dot(&(ray_transformed.d), &cap_n);
   if (d_dot_n != 0)
   {
      p1.x = -ray_transformed.p0.x;
      p1.y = -ray_transformed.p0.y;
      p1.z = 1 - ray_transformed.p0.z;

      l = dot(&(p1), &cap_n) / d_dot_n;
      rayPosition(&ray_transformed, l, &cap_p);
#ifdef DEBUG
      printf("l: %f, lambda: %f, x^2 + y^2: %f\n", l, *lambda, cap_p.px * cap_p.px + cap_p.py * cap_p.py);
#endif

      if (l > THR && cap_p.x * cap_p.x + cap_p.y * cap_p.y <= 1)
      {
         if ((*lambda != -1 && l < *lambda) || *lambda == -1)
         {
            *lambda = l;
            memcpy(n, &cap_n, sizeof(struct point));
         }
      }
   }

#ifdef DEBUG
   printf("lambda: %f\n", *lambda);
#endif

   if (*lambda != -1)
   {
      normalTransform(n, n, cylinder);
      rayPosition(r, *lambda, p);
   }
}

/////////////////////////////////////////////////////////////////
// Surface coordinates & random sampling on object surfaces
/////////////////////////////////////////////////////////////////
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z)
{
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
   // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
   // bottom of the plane.

   struct point p;
   p.x = -2 * a + 1;
   p.y = -2 * b + 1;
   p.z = 0;
   p.w = 1;

   matVecMult(plane->T, &p);

   *x = p.x;
   *y = p.y;
   *z = p.z;
}

void hemiSphereCoordinates(struct point *n, double *x, double *y, double *z)
{
   // Gives a random (x, y, z) coordinates of a point on a hemisphere for
   // the path tracing algorithm of diffuse surfaces
   double a = drand48() * 2 * PI;
   double b = drand48() * 2 - 1;
   b = acos(b);

   struct point p;
   p.x = cos(a) * sin(b);
   p.y = sin(a) * sin(b);
   p.z = cos(b);
   p.w = 1;
   while (dot(n, &p) < 0)
   {
      a = drand48() * 2 * PI;
      b = drand48() * 2 - 1;
      b = acos(b);
      p.x = cos(a) * sin(b);
      p.y = sin(a) * sin(b);
      p.z = cos(b);
      p.w = 1;
   }
   *x = p.x;
   *y = p.y;
   *z = p.z;
}

void sphereCoordinates(struct object3D *sphere, double a, double b, double *x, double *y, double *z)
{
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
   // 'a' in [0, 2*PI] corresponds to the spherical coordinate theta
   // 'b' in [-PI/2, PI/2] corresponds to the spherical coordinate phi

   struct point p;
   p.x = cos(a) * sin(b);
   p.y = sin(a) * sin(b);
   p.z = cos(b);
   p.w = 1;
   matVecMult(sphere->T, &p);

   *x = p.x;
   *y = p.y;
   *z = p.z;
}

void cylCoordinates(struct object3D *cyl, double a, double b, double *x, double *y, double *z)
{
   // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
   // 'a' in [0, 2*PI] corresponds to angle theta around the cylinder
   // 'b' in [0, 1] corresponds to height from the bottom

   struct point p;
   p.x = sin(a);
   p.y = cos(a);
   p.z = b * 2;
   p.w = 1;
   matVecMult(cyl->T, &p);
   *x = p.x;
   *y = p.y;
   *z = p.z;
}

void planeSample(struct object3D *plane, double *x, double *y, double *z)
{
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
   // Sapling should be uniform, meaning there should be an equal change of gedtting
   // any spot on the plane

   double a = drand48();
   double b = drand48();
   planeCoordinates(plane, a, b, x, y, z);
}

void sphereSample(struct object3D *sphere, double *x, double *y, double *z)
{
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the sphere

   double a = drand48() * 2 * PI;
   double b = drand48() * 2 - 1;
   b = acos(b);
   sphereCoordinates(sphere, a, b, x, y, z);
}

void cylSample(struct object3D *cyl, double *x, double *y, double *z)
{
   // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the cylinder
   // Sampling should be uniform over the cylinder.

   double a = drand48() * PI * 2;
   double b = drand48();
   cylCoordinates(cyl, a, b, x, y, z);
}

//////////////////////////////////
// Importance sampling for BRDF
//////////////////////////////////
void cosWeightedSample(struct point *n, struct point *d)
{
   // This function returns a randomly sampled direction over
   // a hemisphere whose pole is the normal direction n. The
   // sampled direction comes from a distribution weighted
   // by the cosine of the angle between n and d.
   // Use this for importance sampling for diffuse surfaces.

   double u1, r, theta, phi;
   double x, y, z, c;
   double v[4][4], R[4][4];
   struct point nz, *cr;
   char line[1024];

   // Random sample on hemisphere with cosine-weighted distribution
   u1 = drand48();
   r = sqrt(u1);
   theta = 2 * PI * drand48();
   x = r * cos(theta);
   y = r * sin(theta);
   z = sqrt(1.0 - (x * x) - (y * y));

   // Need a rotation matrix - start with identity
   memset(&R[0][0], 0, 4 * 4 * sizeof(double));
   R[0][0] = 1.0;
   R[1][1] = 1.0;
   R[2][2] = 1.0;
   R[3][3] = 1.0;

   // Rotation based on cylindrical coordinate conversion
   theta = atan2(n->y, n->x);
   phi = acos(n->z);
   RotateYMat(R, phi);
   RotateZMat(R, theta);

   // Rotate d to align with normal
   d->x = x;
   d->y = y;
   d->z = z;
   d->w = 1.0;
   matVecMult(R, d);

   return;
}

/////////////////////////////////
// Texture mapping functions
/////////////////////////////////
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list)
{
   // Load a texture or normal map image from file and assign it to the
   // specified object.
   // type:   1  ->  Texture map   (RGB, .ppm)
   //         2  ->  Normal map    (RGB, .ppm)
   //         3  ->  Alpha map     (grayscale, .pgm)
   //         4  ->  Intersect map (grayscale, .pgm)
   // Stores loaded images in a linked list to avoid replication
   struct image *im;
   struct textureNode *p;

   if (o != NULL)
   {
      // Check current linked list
      p = *(t_list);
      while (p != NULL)
      {
         if (strcmp(&p->name[0], filename) == 0)
         {
            // Found image already on the list
            if (type == 1)
               o->texImg = p->im;
            else if (type == 2)
               o->normalMap = p->im;
            else if (type == 3)
               o->alphaMap = p->im;
            else
               o->intersectMap = p->im;
            return;
         }
         p = p->next;
      }

      // Load this texture image
      if (type == 1 || type == 2)
         im = readPPMimage(filename);
      else if (type == 3 || type == 4)
         im = readPGMimage(filename);

      // Insert it into the texture list
      if (im != NULL)
      {
         p = (struct textureNode *)calloc(1, sizeof(struct textureNode));
         strcpy(&p->name[0], filename);
         p->type = type;
         p->im = im;
         p->next = NULL;
         // Insert into linked list
         if ((*(t_list)) == NULL)
            *(t_list) = p;
         else
         {
            p->next = (*(t_list))->next;
            (*(t_list))->next = p;
         }
         // Assign to object
         if (type == 1)
            o->texImg = im;
         else if (type == 2)
            o->normalMap = im;
         else if (type == 3)
            o->alphaMap = im;
         else
            o->intersectMap = im;
      }

   } // end if (o != NULL)
}

void texMap(struct image *img, double a, double b, double *R, double *G, double *B)
{
   /*
  Function to determine the colour of a textured object at
  the normalized texture coordinates (a,b).

  a and b are texture coordinates in [0 1].
  img is a pointer to the image structure holding the texture for
   a given object.

  The colour is returned in R, G, B. Uses bi-linear interpolation
  to determine texture colour.
 */
   struct color xc1, xc2;

   a = MIN(1, MAX(0, a));
   b = MIN(1, MAX(0, b));
   a = a * (img->sx - 1);
   b = b * (img->sy - 1);

   int x1 = MAX(0, (int)floor(a));
   int y1 = MAX(0, (int)floor(b));
   int x2 = MIN(img->sx - 1, (int)ceil(a));
   int y2 = MIN(img->sy - 1, (int)ceil(b));

   double *rgbIm = (double *)img->rgbdata;
   double ax1, ax2, ay1, ay2;
   ax1 = (x2 - a) / (x2 - x1);
   ax2 = (a - x1) / (x2 - x1);
   xc1.R = ax1 * ((double)rgbIm[3 * (y1 * img->sx + x1) + 0]) + ax2 * ((double)rgbIm[3 * (y1 * img->sx + x2) + 0]);
   //printf("R x1: %f x2: %f xc1: %f\n", ((double)rgbIm[3*(y1*img->sx + x1) + 0]), ((double)rgbIm[3*(y1*img->sx + x2) + 0]), xc1.R);
   xc1.G = ax1 * ((double)rgbIm[3 * (y1 * img->sx + x1) + 1]) + ax2 * ((double)rgbIm[3 * (y1 * img->sx + x2) + 1]);
   xc1.B = ax1 * ((double)rgbIm[3 * (y1 * img->sx + x1) + 2]) + ax2 * ((double)rgbIm[3 * (y1 * img->sx + x2) + 2]);

   xc2.R = ax1 * ((double)rgbIm[3 * (y2 * img->sx + x1) + 0]) + ax2 * ((double)rgbIm[3 * (y2 * img->sx + x2) + 0]);
   //printf("R x1: %f x2: %f xc1: %f\n", ((double)rgbIm[3*(y1*img->sx + x1) + 0]), ((double)rgbIm[3*(y1*img->sx + x2) + 0]), xc1.R);
   xc2.G = ax1 * ((double)rgbIm[3 * (y2 * img->sx + x1) + 1]) + ax2 * ((double)rgbIm[3 * (y2 * img->sx + x2) + 1]);
   xc2.B = ax1 * ((double)rgbIm[3 * (y2 * img->sx + x1) + 2]) + ax2 * ((double)rgbIm[3 * (y2 * img->sx + x2) + 2]);

   //printf("a: %f b: %f\n", a, b);
   //printf("x: %d y: %d\n", x, y);
   ay1 = (y2 - b) / (y2 - y1);
   ay2 = (b - y1) / (y2 - y1);
   *(R) = MIN(1, MAX(0, ay1 * xc1.R + ay2 * xc2.R));
   *(G) = MIN(1, MAX(0, ay1 * xc1.G + ay2 * xc2.G));
   *(B) = MIN(1, MAX(0, ay1 * xc1.B + ay2 * xc2.B));
   //printf("r: %f g: %f b: %f\n", *R, *G, *B);
   return;
}

void texMapN(struct image *img, double a, double b, double *R, double *G, double *B)
{
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

void alphaMap(struct image *img, double a, double b, double *alpha)
{
   // Just like texture map but returns the alpha value at a,b,
   // notice that alpha maps are single layer grayscale images, hence
   // the separate function.
   double xc1, xc2;

   a = MIN(1, MAX(0, a));
   b = MIN(1, MAX(0, b));
   a = a * (img->sx - 1);
   b = b * (img->sy - 1);

   int x1 = (int)floor(a);
   int y1 = (int)floor(b);
   int x2 = MIN(img->sx - 1, (int)ceil(a));
   int y2 = MIN(img->sy - 1, (int)ceil(b));

   double *rgbIm = (double *)img->rgbdata;
   double ax1, ax2, ay1, ay2;
   ax1 = (x2 - a) / (x2 - x1);
   ax2 = (a - x1) / (x2 - x1);
   xc1 = ax1 * ((double)rgbIm[(y1 * img->sx + x1)]) + ax2 * ((double)rgbIm[(y1 * img->sx + x2)]);
   xc2 = ax1 * ((double)rgbIm[(y2 * img->sx + x1)]) + ax2 * ((double)rgbIm[(y2 * img->sx + x2)]);

   //printf("a: %f b: %f\n", a, b);
   //printf("x: %d y: %d\n", x, y);
   ay1 = (y2 - b) / (y2 - y1);
   ay2 = (b - y1) / (y2 - y1);
   //printf("r: %f g: %f b: %f\n", *R, *G, *B);
   //here

   *(alpha) = ay1 * xc1 + ay2 * xc2;
   //printf("alpha: %f", *alpha);
   return; // with your code if implementing alpha maps.
}

///////////////////////////////////
// Geometric transformation section
///////////////////////////////////

void invert(double *T, double *Tinv)
{
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
   if (U == NULL || s == NULL || V == NULL)
   {
      fprintf(stderr, "Error: Matrix not invertible for this object, returning identity\n");
      memcpy(Tinv, eye4x4, 16 * sizeof(double));
      return;
   }

   // Check for singular matrices...
   for (i = 0; i < 4; i++)
      if (*(s + i) < 1e-9)
         singFlag = 1;
   if (singFlag)
   {
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

void RotateXMat(double T[4][4], double theta)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // X axis.

   double R[4][4];
   memset(&R[0][0], 0, 16 * sizeof(double));

   R[0][0] = 1.0;
   R[1][1] = cos(theta);
   R[1][2] = -sin(theta);
   R[2][1] = sin(theta);
   R[2][2] = cos(theta);
   R[3][3] = 1.0;

   matMult(R, T);
}

void RotateX(struct object3D *o, double theta)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // X axis.

   double R[4][4];
   memset(&R[0][0], 0, 16 * sizeof(double));

   R[0][0] = 1.0;
   R[1][1] = cos(theta);
   R[1][2] = -sin(theta);
   R[2][1] = sin(theta);
   R[2][2] = cos(theta);
   R[3][3] = 1.0;

   matMult(R, o->T);
}

void RotateYMat(double T[4][4], double theta)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Y axis.

   double R[4][4];
   memset(&R[0][0], 0, 16 * sizeof(double));

   R[0][0] = cos(theta);
   R[0][2] = sin(theta);
   R[1][1] = 1.0;
   R[2][0] = -sin(theta);
   R[2][2] = cos(theta);
   R[3][3] = 1.0;

   matMult(R, T);
}

void RotateY(struct object3D *o, double theta)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Y axis.

   double R[4][4];
   memset(&R[0][0], 0, 16 * sizeof(double));

   R[0][0] = cos(theta);
   R[0][2] = sin(theta);
   R[1][1] = 1.0;
   R[2][0] = -sin(theta);
   R[2][2] = cos(theta);
   R[3][3] = 1.0;

   matMult(R, o->T);
}

void RotateZMat(double T[4][4], double theta)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Z axis.

   double R[4][4];
   memset(&R[0][0], 0, 16 * sizeof(double));

   R[0][0] = cos(theta);
   R[0][1] = -sin(theta);
   R[1][0] = sin(theta);
   R[1][1] = cos(theta);
   R[2][2] = 1.0;
   R[3][3] = 1.0;

   matMult(R, T);
}

void RotateZ(struct object3D *o, double theta)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that rotates the object theta *RADIANS* around the
   // Z axis.

   double R[4][4];
   memset(&R[0][0], 0, 16 * sizeof(double));

   R[0][0] = cos(theta);
   R[0][1] = -sin(theta);
   R[1][0] = sin(theta);
   R[1][1] = cos(theta);
   R[2][2] = 1.0;
   R[3][3] = 1.0;

   matMult(R, o->T);
}

void TranslateMat(double T[4][4], double tx, double ty, double tz)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that translates the object by the specified amounts.

   double tr[4][4];
   memset(&tr[0][0], 0, 16 * sizeof(double));

   tr[0][0] = 1.0;
   tr[1][1] = 1.0;
   tr[2][2] = 1.0;
   tr[0][3] = tx;
   tr[1][3] = ty;
   tr[2][3] = tz;
   tr[3][3] = 1.0;

   matMult(tr, T);
}

void Translate(struct object3D *o, double tx, double ty, double tz)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that translates the object by the specified amounts.

   double tr[4][4];
   memset(&tr[0][0], 0, 16 * sizeof(double));

   tr[0][0] = 1.0;
   tr[1][1] = 1.0;
   tr[2][2] = 1.0;
   tr[0][3] = tx;
   tr[1][3] = ty;
   tr[2][3] = tz;
   tr[3][3] = 1.0;

   matMult(tr, o->T);
}

void ScaleMat(double T[4][4], double sx, double sy, double sz)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that scales the object as indicated.

   double S[4][4];
   memset(&S[0][0], 0, 16 * sizeof(double));

   S[0][0] = sx;
   S[1][1] = sy;
   S[2][2] = sz;
   S[3][3] = 1.0;

   matMult(S, T);
}

void Scale(struct object3D *o, double sx, double sy, double sz)
{
   // Multiply the current object transformation matrix T in object o
   // by a matrix that scales the object as indicated.

   double S[4][4];
   memset(&S[0][0], 0, 16 * sizeof(double));

   S[0][0] = sx;
   S[1][1] = sy;
   S[2][2] = sz;
   S[3][3] = 1.0;

   matMult(S, o->T);
   o->LSweight *= (sx * sy * sz); // Update object volume! careful
                                  // won't work for hierarchical
                                  // objects!
}

void printmatrix(double mat[4][4])
{
   fprintf(stderr, "Matrix contains:\n");
   fprintf(stderr, "%f %f %f %f\n", mat[0][0], mat[0][1], mat[0][2], mat[0][3]);
   fprintf(stderr, "%f %f %f %f\n", mat[1][0], mat[1][1], mat[1][2], mat[1][3]);
   fprintf(stderr, "%f %f %f %f\n", mat[2][0], mat[2][1], mat[2][2], mat[2][3]);
   fprintf(stderr, "%f %f %f %f\n", mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
}

/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point *e, struct point *g, struct point *up, double f, double wl, double wt, double wsize)
{
   /*
   This function sets up the camera axes and viewing direction as discussed in the
   lecture notes.
   e - Camera center
   g - Gaze direction
   up - Up vector
   fov - Fild of view in degrees
   f - focal length
 */
   struct view *c;
   struct point *u, *v;

   u = v = NULL;

   // Allocate space for the camera structure
   c = (struct view *)calloc(1, sizeof(struct view));
   if (c == NULL)
   {
      fprintf(stderr, "Out of memory setting up camera model!\n");
      return (NULL);
   }

   // Set up camera center and axes
   c->e.x = e->x; // Copy camera center location, note we must make sure
   c->e.y = e->y; // the camera center provided to this function has pw=1
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
   c->C2W[0][0] = c->u.x;
   c->C2W[1][0] = c->u.y;
   c->C2W[2][0] = c->u.z;
   c->C2W[3][0] = 0;

   c->C2W[0][1] = c->v.x;
   c->C2W[1][1] = c->v.y;
   c->C2W[2][1] = c->v.z;
   c->C2W[3][1] = 0;

   c->C2W[0][2] = c->w.x;
   c->C2W[1][2] = c->w.y;
   c->C2W[2][2] = c->w.z;
   c->C2W[3][2] = 0;

   c->C2W[0][3] = c->e.x;
   c->C2W[1][3] = c->e.y;
   c->C2W[2][3] = c->e.z;
   c->C2W[3][3] = 1;

   // World2Camera matrix (M_wc in the notes)
   // Mind the indexing convention [row][col]
   c->W2C[0][0] = c->u.x;
   c->W2C[1][0] = c->v.x;
   c->W2C[2][0] = c->w.x;
   c->W2C[3][0] = 0;

   c->W2C[0][1] = c->u.y;
   c->W2C[1][1] = c->v.y;
   c->W2C[2][1] = c->w.y;
   c->W2C[3][1] = 0;

   c->W2C[0][2] = c->u.z;
   c->W2C[1][2] = c->v.z;
   c->W2C[2][2] = c->w.z;
   c->W2C[3][2] = 0;

   c->W2C[0][3] = -dot(&c->u, &c->e);
   c->W2C[1][3] = -dot(&c->v, &c->e);
   c->W2C[2][3] = -dot(&c->w, &c->e);
   c->W2C[3][3] = 1;

   free(u);
   free(v);
   return (c);
}

/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename)
{
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
   if (im != NULL)
   {
      im->rgbdata = NULL;
      f = fopen(filename, "rb+");
      if (f == NULL)
      {
         fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
         free(im);
         return (NULL);
      }
      tmpc = fgets(&line[0], 1000, f);
      if (strcmp(&line[0], "P6\n") != 0)
      {
         fprintf(stderr, "Wrong file format, not a .ppm file or header end-of-line characters missing\n");
         free(im);
         fclose(f);
         return (NULL);
      }
      fprintf(stderr, "%s\n", line);
      // Skip over comments
      tmpc = fgets(&line[0], 511, f);
      while (line[0] == '#')
      {
         fprintf(stderr, "%s", line);
         tmpc = fgets(&line[0], 511, f);
      }
      sscanf(&line[0], "%d %d\n", &sizx, &sizy); // Read file size
      fprintf(stderr, "nx=%d, ny=%d\n\n", sizx, sizy);
      im->sx = sizx;
      im->sy = sizy;

      tmpc = fgets(&line[0], 9, f); // Read the remaining header line
      fprintf(stderr, "%s\n", line);
      tmp = (unsigned char *)calloc(sizx * sizy * 3, sizeof(unsigned char));
      fRGB = (double *)calloc(sizx * sizy * 3, sizeof(double));
      if (tmp == NULL || fRGB == NULL)
      {
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

struct image *readPGMimage(const char *filename)
{
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
   if (im != NULL)
   {
      im->rgbdata = NULL;
      f = fopen(filename, "rb+");
      if (f == NULL)
      {
         fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
         free(im);
         return (NULL);
      }
      tmpc = fgets(&line[0], 1000, f);
      if (strcmp(&line[0], "P5\n") != 0)
      {
         fprintf(stderr, "Wrong file format, not a .pgm file or header end-of-line characters missing\n");
         free(im);
         fclose(f);
         return (NULL);
      }
      // Skip over comments
      tmpc = fgets(&line[0], 511, f);
      while (line[0] == '#')
         tmpc = fgets(&line[0], 511, f);
      sscanf(&line[0], "%d %d\n", &sizx, &sizy); // Read file size
      im->sx = sizx;
      im->sy = sizy;

      tmpc = fgets(&line[0], 9, f); // Read the remaining header line
      tmp = (unsigned char *)calloc(sizx * sizy, sizeof(unsigned char));
      fRGB = (double *)calloc(sizx * sizy, sizeof(double));
      if (tmp == NULL || fRGB == NULL)
      {
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

struct image *newImage(int size_x, int size_y)
{
   // Allocates and returns a new image with all zeros. This allocates a double
   // precision floating point image! MIND the difference with the raytracer
   // code that uses 24bpp images.
   struct image *im;

   im = (struct image *)calloc(1, sizeof(struct image));
   if (im != NULL)
   {
      im->rgbdata = NULL;
      im->sx = size_x;
      im->sy = size_y;
      im->rgbdata = (void *)calloc(size_x * size_y * 3, sizeof(double));
      if (im->rgbdata != NULL)
         return (im);
   }
   fprintf(stderr, "Unable to allocate memory for new image\n");
   return (NULL);
}

void imageOutput(struct image *im, const char *filename)
{
   // Writes out a .ppm file from the image data contained in 'im'.
   // Note that Windows typically doesn't know how to open .ppm
   // images. Use Gimp or any other seious image processing
   // software to display .ppm images.
   // Also, note that because of Windows file format management,
   // you may have to modify this file to get image output on
   // Windows machines to work properly.
   //

   FILE *f;
   unsigned char *bits24;
   double *rgbIm;

   if (im != NULL)
      if (im->rgbdata != NULL)
      {
         rgbIm = (double *)im->rgbdata;
         bits24 = (unsigned char *)calloc(im->sx * im->sy * 3, sizeof(unsigned char));
         for (int i = 0; i < im->sx * im->sy * 3; i++)
            *(bits24 + i) = (unsigned char)(255.0 * (*(rgbIm + i)));
         f = fopen(filename, "wb+");
         if (f == NULL)
         {
            fprintf(stderr, "Unable to open file %s for output! No image written\n", filename);
            return;
         }
         fprintf(f, "P6\n");
         fprintf(f, "# Output from PathTracer.c\n");
         fprintf(f, "%d %d\n", im->sx, im->sy);
         fprintf(f, "255\n");
         fwrite(bits24, im->sx * im->sy * 3 * sizeof(unsigned char), 1, f);
         fclose(f);
         return;

         free(bits24);
      }
   fprintf(stderr, "imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im)
{
   // De-allocates memory reserved for the image stored in 'im'
   if (im != NULL)
   {
      if (im->rgbdata != NULL)
         free(im->rgbdata);
      free(im);
   }
}

void dataOutput(double *im, int sx, char *name)
{
   FILE *f;
   double *imT;
   double HDRhist[1000];
   int i, j;
   double mx, mi, biw, pct;
   unsigned char *bits24;
   char pfmname[1024];

   imT = (double *)calloc(sx * sx * 3, sizeof(double));
   memcpy(imT, im, sx * sx * 3 * sizeof(double));
   strcpy(&pfmname[0], name);
   strcat(&pfmname[0], ".pfm");

   // Output the floating point data so we can post-process externally
   f = fopen(pfmname, "w");
   fprintf(f, "PF\n");
   fprintf(f, "%d %d\n", sx, sx);
   fprintf(f, "%1.1f\n", -1.0);
   fwrite(imT, sx * sx * 3 * sizeof(double), 1, f);
   fclose(f);

   // Post processing HDR map - find reasonable cutoffs for normalization
   for (j = 0; j < 1000; j++)
      HDRhist[j] = 0;

   mi = 10e6;
   mx = -10e6;
   for (i = 0; i < sx * sx * 3; i++)
   {
      if (*(imT + i) < mi)
         mi = *(imT + i);
      if (*(imT + i) > mx)
         mx = *(imT + i);
   }

   for (i = 0; i < sx * sx * 3; i++)
   {
      *(imT + i) = *(imT + i) - mi;
      *(imT + i) = *(imT + i) / (mx - mi);
   }
   fprintf(stderr, "Image stats: Minimum=%f, maximum=%f\n", mi, mx);
   biw = 1.000001 / 1000.0;
   // Histogram
   for (i = 0; i < sx * sx * 3; i++)
   {
      for (j = 0; j < 1000; j++)
         if (*(imT + i) >= (biw * j) && *(imT + i) < (biw * (j + 1)))
         {
            HDRhist[j]++;
            break;
         }
   }

   pct = .005 * (sx * sx * 3);
   mx = 0;
   for (j = 5; j < 990; j++)
   {
      mx += HDRhist[j];
      if (HDRhist[j + 5] - HDRhist[j - 5] > pct)
         break;
      if (mx > pct)
         break;
   }
   mi = (biw * (.90 * j));

   for (j = 990; j > 5; j--)
   {
      if (HDRhist[j - 5] - HDRhist[j + 5] > pct)
         break;
   }
   mx = (biw * (j + (.25 * (999 - j))));

   fprintf(stderr, "Limit values chosen at min=%f, max=%f... normalizing image\n", mi, mx);

   for (i = 0; i < sx * sx * 3; i++)
   {
      *(imT + i) = *(imT + i) - mi;
      *(imT + i) = *(imT + i) / (mx - mi);
      if (*(imT + i) < 0.0)
         *(imT + i) = 0.0;
      if (*(imT + i) > 1.0)
         *(imT + i) = 1.0;
      *(imT + i) = pow(*(imT + i), .75);
   }

   bits24 = (unsigned char *)calloc(sx * sx * 3, sizeof(unsigned char));
   for (int i = 0; i < sx * sx * 3; i++)
      *(bits24 + i) = (unsigned char)(255.0 * (*(imT + i)));
   f = fopen(name, "wb+");
   if (f == NULL)
   {
      fprintf(stderr, "Unable to open file %s for output! No image written\n", name);
      return;
   }
   fprintf(f, "P6\n");
   fprintf(f, "# Output from PathTracer.c\n");
   fprintf(f, "%d %d\n", sx, sx);
   fprintf(f, "255\n");
   fwrite(bits24, sx * sx * 3 * sizeof(unsigned char), 1, f);
   fclose(f);
   return;

   free(bits24);
   free(imT);
}

void cleanup(struct object3D *o_list, struct textureNode *t_list)
{
   // De-allocates memory reserved for the object list and for any loaded textures
   // Note that *YOU* must de-allocate any memory reserved for images
   // rendered by the raytracer.
   struct object3D *p, *q;
   struct textureNode *t, *u;

   p = o_list; // De-allocate all memory from objects in the list
   while (p != NULL)
   {
      q = p->next;
      free(p);
      p = q;
   }

   t = t_list; // Delete texture Images
   while (t != NULL)
   {
      u = t->next;
      if (t->im->rgbdata != NULL)
         free(t->im->rgbdata);
      free(t->im);
      free(t);
      t = u;
   }
}

double randn(double mu, double sigma)
{
   // https://en.wikipedia.org/wiki/Marsaglia_polar_method
   double U1, U2, W, mult;
   static double X1, X2;
   static int call = 0;

   if (call == 1)
   {
      call = !call;
      return (mu + sigma * X2);
   }

   do
   {
      U1 = -1 + drand48() * 2;
      U2 = -1 + drand48() * 2;
      W = pow(U1, 2) + pow(U2, 2);
   } while (W >= 1 || W == 0);

   mult = sqrt((-2 * log(W)) / W);
   X1 = U1 * mult;
   X2 = U2 * mult;

   call = !call;

   return (mu + sigma * X1);
}

void RGB2Hue(double *H, double R, double G, double B)
{
   //https://stackoverflow.com/questions/23090019/fastest-formula-to-get-hue-from-rgb

   double col_max; //0 = red is max, 1 = green is max, 2 = blue is max
   double min = 10, max = -10;
   //find max color
   if (R > max)
   {
      max = R;
      col_max = 0;
   }
   if (G > max)
   {
      max = G;
      col_max = 1;
   }
   if (B > max)
   {
      max = B;
      col_max = 2;
   }

   //find min color
   if (R < min)
   {
      min = R;
   }
   if (G < min)
   {
      min = G;
   }
   if (B < min)
   {
      min = B;
   }

   //calculate Hue
   if (col_max == 0)
   { //red is max
      *H = (G - B) / (max - min);
   }
   else if (col_max == 1)
   { //green is max
      *H = 2.0 + (B - R) / (max - min);
   }
   else
   { //blue is max
      *H = 4.0 + (R - G) / (max - min);
   }

   if (*H < 0)
   {
      *H += 6;
   }

   *H /= 6;

   // if(!(0 <= *H && *H <= 1)){
   //   printf("hue: %f\n", *H);
   // }
}

void hue2RGB(double H, double *R, double *G, double *B)
{
   /* Given a HUE value in [0 1], with 0 for deep red and
  * 1 for purple, obtains the corresponding RGB values
  * that give the equivalent colour
  */

   double C, X;

   C = 1.0;
   X = C * (1.0 - fabs(fmod(6.0 * H, 2.0) - 1.0));

   if (H < 1.0 / 6.0)
   {
      *R = 1.0;
      *G = X;
      *B = 0;
   }
   else if (H < 2.0 / 6.0)
   {
      *R = X;
      *G = C;
      *B = 0;
   }
   else if (H < 3.0 / 6.0)
   {
      *R = 0;
      *G = C;
      *B = X;
   }
   else if (H < 4.0 / 6.0)
   {
      *R = 0;
      *G = X;
      *B = C;
   }
   else if (H < 5.0 / 6.0)
   {
      *R = X;
      *G = 0;
      *B = C;
   }
   else
   {
      *R = C;
      *G = 0;
      *B = X;
   }
}

double hueToRef(double H, double r_idx)
{
   double color_spread = 0.6;
   return r_idx + color_spread * (H - 0.5);
}

void biDirectionalSample(struct object3D *light, struct ray3D *cam_ray, struct color *retcol)
{
   int ldepth = 5;

   struct ray3D light_ray;
   light_ray.rayPos = &rayPosition;
   double x, y, z;
   light->randomPoint(light, &x, &y, &z);
   light_ray.p0.x = x;
   light_ray.p0.y = y;
   light_ray.p0.z = z;
   light_ray.p0.w = 1;

   light_ray.d.x = drand48();
   light_ray.d.y = drand48();
   light_ray.d.z = drand48();
   light_ray.d.w = 1;
   normalize(&light_ray.d);

   light_ray.R = light->col.R;
   light_ray.G = light->col.G;
   light_ray.B = light->col.B;
   light_ray.Ir = 0;
   light_ray.Ig = 0;
   light_ray.Ib = 0;
   light_ray.isLightRay = 1;

   struct color col; // colour for the light ray
   PathTrace(&light_ray, 1, &col, NULL, NULL);
   struct ray3D Cam_to_Light;
   Cam_to_Light.p0 = cam_ray->p0;
   Cam_to_Light.d.x = light_ray.p0.x - cam_ray->p0.x;
   Cam_to_Light.d.y = light_ray.p0.y - cam_ray->p0.y;
   Cam_to_Light.d.z = light_ray.p0.z - cam_ray->p0.z;
   Cam_to_Light.d.w = 1;
   double lambda;
   struct object3D *obstruction = NULL;
   struct point lightp;
   struct point nls;
   double La, Lb;

   findFirstHit(&Cam_to_Light, &lambda, NULL, &obstruction, &lightp, &nls, &La, &Lb);
   if (THR < lambda && lambda < 1)
   { //the path is obstructed.
      retcol->R = 0;
      retcol->G = 0;
      retcol->B = 0;
      return;
   }

   double dxd = Cam_to_Light.d.x * Cam_to_Light.d.x + Cam_to_Light.d.y * Cam_to_Light.d.y + Cam_to_Light.d.pzz Cam_to_Light.d.pzz
   normalize(&Cam_to_Light.d);
   double nls_dot_l = fabs(dot(&nls, &Cam_to_Light.d));
   double w = MIN(1, (nls_dot_l) / (dxd));
   retcol->R = w * col.R;
   retcol->G = w * col.G;
   retcol->B = w * col.B;
}