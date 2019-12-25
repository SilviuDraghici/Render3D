/*
   L-system code
   
   Silviu Draghici
*/

#define b_to_f 0.666
#define stemR 0.2
#define stemG 0.9
#define stemB 0.25

void newTree(struct object3D **object_list, double hierarchyMat[4][4], struct color *col, double distFromC, double numBranches, double maxRotation, double depth)
{
   struct object3D *o;

   //this is what the default shape of the cylinder should be
   //not part of the hierarchy
   double yscale = 2;
   double cylTransforms[4][4];
   memcpy(cylTransforms, eye4x4, 16 * sizeof(double));
   RotateXMat(cylTransforms, -PI / 2);
   ScaleMat(cylTransforms, 0.5, yscale, 0.5);

   o = newCyl(.05, .95, 0, 0, stemR, stemG, stemB, 1, 1, 6);
   memcpy(o->T, cylTransforms, 16 * sizeof(double));

   matMult(hierarchyMat, o->T);

   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, object_list);
   if (1 >= depth)
   {
      return;
   }

   double newTransforms[4][4];
   double angle = drand48() * 2 * maxRotation;
   for (int i = 0; i < numBranches; i++)
   {
      memcpy(newTransforms, eye4x4, 16 * sizeof(double));
      ScaleMat(newTransforms, 0.8, 0.8, 0.8);
      TranslateMat(newTransforms, 0, yscale / 2 + 0.5, 0);

      RotateXMat(newTransforms, PI / 4);
      RotateYMat(newTransforms, angle + i * 2 * PI / numBranches);
      TranslateMat(newTransforms, 0, yscale / 2 + 0.5, 0);
      matMult(hierarchyMat, newTransforms);
      newBranch(object_list, newTransforms, col, numBranches, maxRotation, 1, depth);
   }
}

void newBranch(struct object3D **object_list, double hierarchyMat[4][4], struct color *col, double numBranches, double maxRotation, double curr_depth, double depth)
{

   //this is what the default shape of the cylinder should be
   //not part of the hierarchy
   double yscale = 2;
   double cylTransforms[4][4];
   memcpy(cylTransforms, eye4x4, 16 * sizeof(double));
   RotateXMat(cylTransforms, -PI / 2);
   ScaleMat(cylTransforms, 0.5, yscale, 0.5);

   struct object3D *o;
   o = newCyl(.05, .95, 0, 0, stemR, stemG, stemB, 1, 1, 6);
   memcpy(o->T, cylTransforms, 16 * sizeof(double));

   matMult(hierarchyMat, o->T);

   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, object_list);

   if (curr_depth + 1 >= depth)
   {
      return;
   }

   double newTransforms[4][4];
   double angle = drand48() * 2 * maxRotation;
   double dice;
   struct color f;
   for (int i = 0; i < numBranches; i++)
   {
      memcpy(newTransforms, eye4x4, 16 * sizeof(double));
      ScaleMat(newTransforms, 0.8, 0.8, 0.8);
      TranslateMat(newTransforms, 0, yscale / 2 + 0.5, 0);

      RotateXMat(newTransforms, PI / 4);
      RotateYMat(newTransforms, angle + i * 2 * PI / numBranches);
      TranslateMat(newTransforms, 0, yscale / 2 + 0.5, 0);
      matMult(hierarchyMat, newTransforms);

      dice = drand48();
      if (dice <= b_to_f)
      { //new flower
         newFlower(object_list, col, newTransforms);
      }
      else
      {
         newBranch(object_list, newTransforms, col, numBranches, maxRotation, curr_depth + 1, depth);
      }
   }
}

void newFlower(struct object3D **object_list, struct color *petalCol, double hierarchyMat[4][4])
{
   //printf("new flower!\n");
   struct color f;
   if (petalCol == NULL)
   {
      //printf("making color\n");
      f.R = drand48();
      f.G = drand48();
      f.B = drand48();
      petalCol = &f;
   }
   //printf("flower.R = %f, flower.G = %f, flower.B = %f;\n", petalCol->R,petalCol->G, petalCol->B);
   struct object3D *o;
   //stem
   double fscale = 3; // the initial flower was too small
   double stemscale = 0.5;
   o = newCyl(.05, .95, 0, 0, stemR, stemG, stemB, 1, 1, 6);
   RotateX(o, -PI / 2);
   Scale(o, fscale * 0.1, fscale * stemscale, fscale * 0.1);

   matMult(hierarchyMat, o->T);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, object_list);

   //newSphere(ra,  rd, rs, rg,  r, g, b, alpha, r_index, shiny)
   o = newSphere(.05, .95, .95, .75, 1 - petalCol->R, 1 - petalCol->G, 1 - petalCol->B, 1, 1, 6);
   Scale(o, fscale * 0.3, fscale * 0.1, fscale * 0.3);
   Translate(o, 0, fscale * stemscale, 0);
   matMult(hierarchyMat, o->T);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, object_list);

   for (int i = 0; i < 5; i++)
   {
      o = newSphere(.05, .95, 0, 0, petalCol->R, petalCol->G, petalCol->B, 1, 1, 6);
      Translate(o, 0.5, 0, 0);
      Scale(o, fscale * 0.5, fscale * 0.005, fscale * 0.3);
      RotateZ(o, 0.25);
      Translate(o, 0.1, fscale * stemscale - 0.04, 0);
      RotateY(o, i * 2 * PI / 5);
      matMult(hierarchyMat, o->T);
      invert(&o->T[0][0], &o->Tinv[0][0]);
      insertObject(o, object_list);
   }
}
