void buildScene(void)
{
 struct object3D *o;
 struct point p;

 struct textureNode *t_list = NULL;

 const char *file = "marble.ppm";
 const char *jupfile = "jupiter.ppm";
 const char *nfile = "stone_normal.ppm";
 const char *alphafile = "earthalpha.pgm";

 const char *fence = "fence.ppm";
 const char *fencenormal = "fencenormal.ppm";
 const char *fencealpha = "fencealpha.pgm";
 
 // Cornell box
 o=newPlane(1.0,0.0,0.0,.75,.25,.25,.05,1.4);	// Left
 o->label[0] = 'L';
 o->label[1] = 'E';
 o->label[2] = 'F';
 o->label[3] = 'T';
 o->label[4] = '\0';
 loadTexture(o, nfile, 2, &t_list);
 RotateY(o,PI/2);
 Scale(o,25,25,25);
 Translate(o,-10,0,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.25,.25,.75,.05,1.4);		// Right
 o->label[0] = 'R';
 o->label[1] = 'I';
 o->label[2] = 'G';
 o->label[3] = 'H';
 o->label[4] = 'T';
 o->label[5] = '\0';
 RotateY(o,PI/2);
 Scale(o,25,25,25);
 Translate(o,10,0,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.75,.75,.75,.05,1.4);		// Back
 o->label[0] = 'B';
 o->label[1] = 'A';
 o->label[2] = 'C';
 o->label[3] = 'K';
 o->label[4] = '\0';
 loadTexture(o, file, 1, &t_list);
 //RotateZ(o, PI/4);
 Scale(o,9,9,9);
 Translate(o,0,0,15);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.75,.75,.75,.02,1.4);	// Bottom
 o->label[0] = 'B';
 o->label[1] = 'O';
 o->label[2] = 'T';
 o->label[3] = 'T';
 o->label[4] = 'O';
 o->label[5] = 'M';
 o->label[6] = '\0';
 RotateX(o,PI/2);
 Scale(o,25,25,25);
 Translate(o,0,-10,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(1.0,0.0,0.0,.75,.75,.75,.05,1.4);		// Top
 o->label[0] = 'T';
 o->label[1] = 'O';
 o->label[2] = 'P';
 o->label[3] = '\0';
 RotateX(o,PI/2);
 Scale(o,500,500,500);
 Translate(o,0,10,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // Two spheres scene
 o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,1.54);		// Refract
 loadTexture(o, fence, 1, &t_list);//texture map
 loadTexture(o, fencenormal, 2, &t_list);//normal map
 loadTexture(o, fencealpha, 3, &t_list);//alpha map
 loadTexture(o, fencealpha, 4, &t_list);//intersect map
 RotateY(o, PI);
 Scale(o,3.75,3.75,3.75);
 Translate(o,-5,-4.0,4.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,1.54);		// Refract
 loadTexture(o, fence, 1, &t_list);//texture map
 //loadTexture(o, fencenormal, 2, &t_list);//normal map
 loadTexture(o, fencealpha, 3, &t_list);//alpha map
 //loadTexture(o, fencealpha, 4, &t_list);//intersect map
 RotateY(o, PI);
 Scale(o,3.75,3.75,3.75);
 Translate(o,-5,6.0,8.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(0.0,1.0,0.0,.99,.99,.99,.05,2.47);		// Reflect
 o->label[0] = 'S';
 o->label[1] = 'P';
 o->label[2] = 'H';
 o->label[3] = '_';
 o->label[4] = 'R';
 o->label[5] = 'H';
 o->label[6] = 'T';
 o->label[7] = '\0';
 loadTexture(o, nfile, 2, &t_list);
 Scale(o,3.75,3.75,3.75);
 Translate(o,4,-3.75,6.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);		// Jupiter
 o->label[0] = 'J';
 o->label[1] = 'U';
 o->label[2] = 'P';
 o->label[3] = 'I';
 o->label[4] = 'T';
 o->label[5] = 'E';
 o->label[6] = 'R';
 o->label[7] = '\0';
 loadTexture(o, jupfile, 1, &t_list);
 Scale(o,2,2,2);
 Translate(o,7,7,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // Planar light source at top
 o=newPlane(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
 o->label[0] = 'L';
 o->label[1] = '1';
 o->label[2] = '\0';
 Scale(o,.5,2.5,1);
 RotateX(o,PI/2);
 Translate(o,0,9.995,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 o->isLightSource=1;
 insertObject(o,&object_list);

 // Planar light source at bottom
 o=newPlane(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
 o->label[0] = 'L';
 o->label[1] = '2';
 o->label[2] = '\0';
 Scale(o,.5,2.5,1);
 RotateX(o,-PI/2);
 Translate(o,0,-9.995,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 o->isLightSource=1;
 //insertObject(o,&object_list);
}