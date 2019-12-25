void buildScene(void)
{
 struct object3D *o;
 struct point3D p;

 struct textureNode *t_list = NULL;

 const char *file = "marble.ppm";
 const char *jupfile = "jupiter.ppm";
 const char *nfile = "stone_normal.ppm";
 const char *alphafile = "alphame.pgm";
 
 // Cornell box
 o=newPlane(1.0,0.0,0.0,0.75,0.25,0.25,.05,1.4);	// Left
 strcpy(o->label,"Left\0");
 RotateY(o,PI/2);
 Scale(o,25,25,25);
 Translate(o,-10,0,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.25,.25,.75,.05,1.4);		// Right
 strcpy(o->label,"Right\0");
 RotateY(o,PI/2);
 Scale(o,25,25,25);
 Translate(o,10,0,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.75,.75,.75,.05,1.4);		// Back
 strcpy(o->label,"Back\0");
 loadTexture(o, file, 1, &t_list);
 Scale(o,10,10,10);
 Translate(o,0,0,15);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.75,.75,.75,.02,1.4);	// Bottom
 strcpy(o->label,"Bottom\0");
 RotateX(o,PI/2);
 Scale(o,25,25,25);
 Translate(o,0,-10,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.75,.75,.75,.05,1.4);		// Top
 strcpy(o->label,"Top\0");
 RotateX(o,PI/2);
 Scale(o,500,500,500);
 Translate(o,0,10,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // Planar light source at top
 o=newPlane(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
 strcpy(o->label,"L1\0");
 Scale(o,.5,2.5,1);
 RotateX(o,PI/2);
 Translate(o,0,9.995,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 o->isLightSource=1;
 //insertObject(o,&object_list);


 //Laptop components
 const char *pupfile = "laptop/puppywallpaper.ppm";
 const char *keyboard = "laptop/keyboard.ppm";
 const char *keyboardn = "laptop/keyboardNormal.ppm";
 char laptopobj[] = "laptop/Lowpoly_Notebook_2.obj";

 double laptopMat[4][4] = {{1.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0}};
 ScaleMat(laptopMat, 4, 4, 4);
 RotateYMat(laptopMat, PI/3);
 //TranslateMat(laptopMat, 0, -8, 5);

//  struct object3D *m;
//  m = newMesh(laptopobj, 1.00, 0.00, 0.0, 1.0, 1.0, 1.0, 0.0, 1.54);
//  invert(&m->T[0][0],&m->Tinv[0][0]);
//  insertObject(m,&object_list);

//  double xscale = m->max->px - m->min->px;
//     double yscale = m->max->py - m->min->py;
//     double zscale = m->max->pz - m->min->pz;

//     double xtrans = (m->max->px - m->min->px)/2;
//     double ytrans = (m->max->py - m->min->py)/2;
//     double ztrans = (m->max->pz - m->min->pz)/2;

//     printf("max x %f y %f z %f\n", m->max->px, m->max->py, m->max->pz);
//     printf("min x %f y %f z %f\n", m->min->px, m->min->py, m->min->pz);

//     printf("scale x %f y %f z %f\n", xscale, yscale, zscale);
//     o=newBox(1.00, 0.00, 0.0, 1.0, 1.0, 1.0, 0.0, 1.54);
//     //Translate(o, xtrans, ytrans, ztrans);
//     memcpy(&o->T[0][0], &m->T[0][0], 16 * sizeof(double));

//     struct point3D translate;
//     translate.px = m->T[0][3];
//     translate.py = m->T[1][3];
//     translate.pz = m->T[2][3];

//     o->T[0][3] = 0;
//     o->T[1][3] = 0;
//     o->T[2][3] = 0;
//     Scale(o, xscale, yscale, zscale);

    // o->T[0][3] = translate.px;
    // o->T[1][3] = translate.py;
    // o->T[2][3] = translate.pz;
    // Translate(o, 0, ytrans, 0);
    // //Translate(o, translate.px, translate.py, translate.pz);

    // invert(&o->T[0][0], &o->Tinv[0][0]);
 
 o=newPlane(1.0,0.0,0.0,1,1,1,.05,1.4);	// Screen
 strcpy(o->label,"Screen\0");
 o->isLightSource=1;
 printf("label %s\n", o->label);
 loadTexture(o, pupfile, 1, &t_list);
 //RotateY(o, -PI/4);
 //RotateX(o, -PI/6);
 Scale(o,1.7,1,1);
 //Translate(o,0,0,5);
 matMult(laptopMat, o->T);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.05,.05,.05,.05,1.4);	// behind screen
 strcpy(o->label,"bScreen\0");
 printf("label %s\n", o->label);
 Scale(o,1.9,1.2,1.2);
 Translate(o,0,0,0.00001);
 matMult(laptopMat, o->T);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 //insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.75,.25,.25,.05,1.4);	// keys
 strcpy(o->label,"Keys\0");
 loadTexture(o, keyboard, 1, &t_list);
 loadTexture(o, keyboardn, 2, &t_list);
 //RotateY(o, -PI/2);
 RotateX(o, PI/2);
 //Scale(o,5,5,5);
 //Translate(o,0,0,-1);
 matMult(laptopMat, o->T);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 //insertObject(o,&object_list);
}