/*
   L-system code
   
   Silviu Draghici
*/

#define b_to_f 0.666
#define stemR 0.2
#define stemG 0.9
#define stemB 0.25

void newTree(struct object3D **object_list, struct transform hierarchy, struct color *col, double distFromC, double numBranches, double maxRotation, double depth) {
   struct object3D *o;

   //this is what the default shape of the cylinder should be
   //not part of the hierarchy
   double yscale = 2;
   struct transform cylTransforms;
   cylTransforms = Sc(0.5, yscale, 0.5) * RotX(-PI / 2);

   o = newCyl(.05, .95, 0, 0, stemR, stemG, stemB, 1, 1, 6);
   o->T = hierarchy * cylTransforms;
   invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
   insertObject(o, object_list);

   if (1 >= depth) {
      return;
   }

   transform branch = RotX(PI / 4) * Tr(0, (yscale / 2) + 0.5, 0) * Sc(0.8);
   transform newTransforms;
   double angle = drand48() * 2 * maxRotation;
   for (int i = 0; i < numBranches; i++) {
      newTransforms = branch;
      newTransforms *= Tr(0, (yscale / 2) + 0.5, 0) * RotY(angle + i * 2 * PI / numBranches);
      newTransforms *= hierarchy;
      newBranch(object_list, newTransforms, col, numBranches, maxRotation, 1, depth);
   }
}

void newBranch(struct object3D **object_list, struct transform hierarchy, struct color *col, double numBranches, double maxRotation, double curr_depth, double depth) {
   struct object3D *o;

   //this is what the default shape of the cylinder should be
   //not part of the hierarchy
   double yscale = 2;
   struct transform cylTransforms;
   cylTransforms = Sc(0.5, yscale, 0.5) * RotX(-PI / 2);

   o = newCyl(.05, .95, 0, 0, stemR, stemG, stemB, 1, 1, 6);
   o->T = hierarchy * cylTransforms;
   invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
   insertObject(o, object_list);

   if (curr_depth + 1 >= depth) {
      return;
   }

   transform branch = RotX(PI / 4) * Tr(0, (yscale / 2) + 0.5, 0) * Sc(0.8);
   transform newTransforms;
   double angle = drand48() * 2 * maxRotation;
   double dice;
   struct color f;
   for (int i = 0; i < numBranches; i++) {
      newTransforms = branch;
      newTransforms *= Tr(0, (yscale / 2) + 0.5, 0) * RotY(angle + i * 2 * PI / numBranches);
      newTransforms *= hierarchy;

      dice = drand48();
      if (dice <= b_to_f) {  //new flower
         newFlower(object_list, col, newTransforms);
      } else {
         newBranch(object_list, newTransforms, col, numBranches, maxRotation, curr_depth + 1, depth);
      }
   }
}

void newFlower(struct object3D **object_list, struct color *petalCol, struct transform hierarchy) {
   //printf("new flower!\n");
   struct color f;
   if (petalCol == NULL) {
      //printf("making color\n");
      f.R = drand48();
      f.G = drand48();
      f.B = drand48();
      petalCol = &f;
   }
   //printf("flower.R = %f, flower.G = %f, flower.B = %f;\n", petalCol->R,petalCol->G, petalCol->B);
   struct object3D *o;

   //stem
   double fscale = 3;  // the initial flower was too small
   double stemscale = 0.5;
   o = newCyl(.05, .95, 0, 0, stemR, stemG, stemB, 1, 1, 6);
   o->T = Sc(fscale * 0.1, fscale * stemscale, fscale * 0.1) * RotX(-PI / 2);
   o->T *= hierarchy;
   invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
   insertObject(o, object_list);

   //newSphere(ra,  rd, rs, rg,  r, g, b, alpha, r_index, shiny)
   o = newSphere(.05, .95, .95, .75, 1 - petalCol->R, 1 - petalCol->G, 1 - petalCol->B, 1, 1, 6);
   o->T = Tr(0, fscale * stemscale, 0) * Sc(fscale * 0.3, fscale * 0.1, fscale * 0.3);
   o->T *= hierarchy;
   invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
   insertObject(o, object_list);

   transform petal = Sc(fscale * 0.5, fscale * 0.005, fscale * 0.3) * Tr(0.5, 0, 0);
   petal *= Tr(0.1, fscale * stemscale - 0.04, 0) * RotZ(0.25);
   for (int i = 0; i < 5; i++) {
      o = newSphere(.05, .95, 0, 0, petalCol->R, petalCol->G, petalCol->B, 1, 1, 6);
      o->T = RotY(i * 2 * PI / 5) * petal;
      o->T *= hierarchy;
      invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
      insertObject(o, object_list);
   }
}
