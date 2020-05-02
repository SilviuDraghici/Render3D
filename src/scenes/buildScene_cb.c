
 Object *o;
 struct point p;
 struct pointLS *l;

 struct textureNode *t_list = NULL;
 
 // Cornell box

   cam_pos = point(0, 0, -15);
   //cam_gaze_point = point(0, 0, 0);
   //cam_gaze = cam_gaze_point - cam_pos;
   //cam_up = point(0, 1, 0);
   cam_focal = -3;

 // Left
 o = new Plane(.75,.25,.25);
 o->set_pathTrace_properties(1.0, 0.0, 0.0);
 o->r_index = 1.4;
 strcpy(o->label, "Left Wall");
 //loadTexture(o, nfile, 2, &t_list);
 o->T *= RotY(PI / 2);
 o->T *= Sc(25);
 o->T *= Tr(-10, 0, 5);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 insertObject(o,&object_list);

 // Right
 o = new Plane(.25,.25,.75);
 o->set_pathTrace_properties(1.0, 0.0, 0.0);
 o->r_index = 1.4;
 strcpy(o->label, "Right Wall");
 o->T *= RotY(PI / 2);
 o->T *= Sc(25);
 o->T *= Tr(10, 0, 5);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 insertObject(o,&object_list);

 // Back
 o = new Plane(.75,.75,.75);
 o->set_pathTrace_properties(1.0, 0.0, 0.0);
 o->r_index = 1.4;
 strcpy(o->label, "Back Wall");
 //loadTexture(o, file, 1, &t_list);
 //o->T *= RotateZ(o, PI/4);
 o->T *= Sc(10);
 o->T *= Tr(0,0,15);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 insertObject(o,&object_list);

 // Bottom
 o = new Plane(.75,.75,.75);
 o->set_pathTrace_properties(1.0, 0.0, 0.0);
 o->r_index = 1.4;
 strcpy(o->label, "Bottom Wall");
 o->T *= RotX(PI/2);
 o->T *= Sc(25);
 o->T *= Tr(0,-10,5);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 insertObject(o,&object_list);

 // Top
 o = new Plane(.75,.75,.75);
 o->set_pathTrace_properties(1.0, 0.0, 0.0);
 o->r_index = 1.4;
 strcpy(o->label, "Top Wall");
 o->T *= RotX(-PI/2);
 o->T *= Sc(25);
 o->T *= Tr(0,10,5);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 insertObject(o,&object_list);

 // Two spheres scene
 // Refract
 o = new Sphere(.99,.99,.99);
 strcpy(o->label, "Refract Shpere");
 o->set_pathTrace_properties(0.0, 0.0, 1.0);
 o->r_index = 1.54;

 o->T *= RotY(PI);
 o->T *= Sc(3.75);
 o->T *= Tr(-5,-4.0,4.5);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 insertObject(o,&object_list);

 // Reflect
 o = new Sphere(.99,.99,.99);
 strcpy(o->label, "Reflect Sphere");
 o->set_pathTrace_properties(0.0, 1.0, 0.0);
 o->refl_sig = 0.05;
 o->r_index = 2.47;
 strcpy(o->label, "Right Sphere");

 o->T *= Sc(3.75);
 o->T *= Tr(4,-3.75,6.5);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 insertObject(o,&object_list);

 // Planar light source at top
 o = new Plane(1.0,1.0,1.0);
 o->set_pathTrace_properties(1.0, 0.0, 0.0);
 o->refl_sig = 0.0;
 o->r_index = 1.54;
 strcpy(o->label, "Top Light");
 o->T *= Sc(0.5, 2.5, 1);
 o->T *= RotX(PI/2);
 o->T *= Tr(0,9.995,5);
 invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
 o->isLightSource=1;
 o->pt.LSweight *= 0.5 * 2.5 * 1; // <- scale weight by scale
 insertObject(o,&object_list);

 p.x = 0;
    p.y = 9.9;
    p.z = 5;
    p.w = 1;
    l = newPLS(&p, .95, .95, .95);
    insertPLS(l, &light_list);
