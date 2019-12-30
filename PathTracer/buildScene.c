void buildScene(void)
{
   // Sets up all objets in the scene. This involves creating each object,
   // defining the transformations needed to shape and position it as
   // desired, specifying the reflectance properties (albedos and colours)
   // and setting up textures where needed.
   //
   // NOTE: Light sources are now EXCLUSIVELY area light sources. They are
   //       defined by regular objects whose 'isLightSource' flag is set
   //       to 1. Therefore, you can create light sources with any shape
   //       and colour using the same object primitives and transforms
   //       you're using to set up the scene.
   //
   // To create hierarchical objects:
   //    You must keep track of transformations carried out by parent objects
   //    as you move through the hierarchy. Declare and manipulate your own
   //    transformation matrices (use the provided functions in utils.c to
   //    compound transformations on these matrices). When declaring a new
   //    object within the hierarchy
   //    - Initialize the object
   //    - Apply any object-level transforms to shape/rotate/resize/move
   //      the object using regular object transformation functions
   //    - Apply the transformations passed on from the parent object
   //      by pre-multiplying the matrix containing the parent's transforms
   //      with the object's own transformation matrix.
   //    - Compute and store the object's inverse transform as usual.
   //
   // NOTE: After setting up the transformations for each object, don't
   //       forget to set up the inverse transform matrix!

   struct object3D *o;
   struct point p;

   // Cornell box
   o = newSphere(1.0, 0.0, 0.0, .75, .25, .25, .05, 1.4); // Left
   o->label[0] = 'L';
   o->label[1] = 'E';
   o->label[2] = 'F';
   o->label[3] = 'T';
   o->label[4] = '\0';
   Scale(o, 500, 500, 500);
   Translate(o, -510, 0, 5);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, &object_list);

   o = newSphere(1.0, 0.0, 0.0, .25, .25, .75, .05, 1.4); // Right
   o->label[0] = 'R';
   o->label[1] = 'I';
   o->label[2] = 'G';
   o->label[3] = 'H';
   o->label[4] = 'T';
   o->label[5] = '\0';
   Scale(o, 500, 500, 500);
   Translate(o, 510, 0, 5);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, &object_list);

   o = newSphere(1.0, 0.0, 0.0, .75, .75, .75, .05, 1.4); // Back
   o->label[0] = 'B';
   o->label[1] = 'A';
   o->label[2] = 'C';
   o->label[3] = 'K';
   o->label[4] = '\0';
   Scale(o, 500, 500, 500);
   Translate(o, 0, 0, 515);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, &object_list);

   o = newSphere(1.0, 0.0, 0.0, .75, .75, .75, .02, 1.4); // Bottom
   o->label[0] = 'B';
   o->label[1] = 'O';
   o->label[2] = 'T';
   o->label[3] = 'T';
   o->label[4] = 'O';
   o->label[5] = 'M';
   o->label[6] = '\0';
   Scale(o, 500, 500, 500);
   Translate(o, 0, -510, 5);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, &object_list);

   o = newSphere(1.0, 0.0, 0.0, .75, .75, .75, .05, 1.4); // Top
   o->label[0] = 'T';
   o->label[1] = 'O';
   o->label[2] = 'P';
   o->label[3] = '\0';
   Scale(o, 500, 500, 500);
   Translate(o, 0, 510, 5);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, &object_list);

   // Two spheres scene
   o = newSphere(0.0, 0.0, 1.0, .99, .99, .99, .01, 1.54); // Refract
   o->label[0] = 'S';
   o->label[1] = 'P';
   o->label[2] = 'H';
   o->label[3] = '_';
   o->label[4] = 'L';
   o->label[5] = 'F';
   o->label[6] = 'T';
   o->label[7] = '\0';
   Scale(o, 3.75, 3.75, 3.75);
   Translate(o, -5, -4.0, 4.5);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, &object_list);

   o = newSphere(0.0, 1.0, 0.0, .99, .99, .99, .05, 2.47); // Reflect
   o->label[0] = 'S';
   o->label[1] = 'P';
   o->label[2] = 'H';
   o->label[3] = '_';
   o->label[4] = 'R';
   o->label[5] = 'H';
   o->label[6] = 'T';
   o->label[7] = '\0';
   Scale(o, 3.75, 3.75, 3.75);
   Translate(o, 4, -3.75, 6.5);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   insertObject(o, &object_list);

   /*
 // Ring of refracting spheres
 for (int i=0; i<5;i++)
 {
  o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,1.45+(.1*i));
  Scale(o,1.75,1.75,1.75);
  Translate(o,3.25*cos(2*PI*i/5),-2.45,3+3.25*sin(2*PI*i/5));
  invert(&o->T[0][0],&o->Tinv[0][0]);
  insertObject(o,&object_list);
 }

 for (int i=0; i<7;i++)
 {
  o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,2.00+(.05*i));
  Scale(o,1.75,1.75,1.75);
  Translate(o,4.60*cos(2*PI*i/7),-6.35,3+4.60*sin(2*PI*i/7));
  invert(&o->T[0][0],&o->Tinv[0][0]);
  insertObject(o,&object_list);
 }
*/

   // Planar light source at top
   o = newPlane(1.00, 0.00, 0.0, 1.0, 1.0, 1.0, 0.0, 1.54);
   o->label[0] = 'L';
   o->label[1] = '1';
   o->label[2] = '\0';
   Scale(o, .5, 2.5, 1);
   RotateX(o, PI / 2);
   Translate(o, 0, 9.995, 5);
   invert(&o->T[0][0], &o->Tinv[0][0]);
   o->isLightSource = 1;
   insertObject(o, &object_list);
}