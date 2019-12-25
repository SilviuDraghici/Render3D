
 struct object3D *o;
 struct pointLS *l;
 struct textureNode *t_list = NULL;
 struct point3D p, pl;
 p.pw = 1, pl.pw = 1;

 double f = 1 ;

 const char* file = "grass.ppm";

 double angle;
 if(iter_val != -1){
	 angle = iter_val*2*PI/num_frames/3;
 }else{
	 angle = 0;
 }

 double cameraMat[4][4] = {{1.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0}};
 RotateYMat(cameraMat, angle);
 TranslateMat(cameraMat, 0, -8, 5);

  struct color flower;
  //blue & copper      <-- nice color
  //flower.R = 0.019842, flower.G = 0.378377, flower.B = 0.678877;

  //darker red
  //flower.R = 0.364602, flower.G = 0.091331, flower.B = 0.092298;

  //red & cyan         <-- nice color
  flower.R = 0.767051, flower.G = 0.018915, flower.B = 0.252360;

  //light purple and green
  //flower.R = 0.875981, flower.G = 0.531557, flower.B = 0.920261;

  //purple and green   <-- nice color
  //flower.R = 0.810429, flower.G = 0.188420, flower.B = 0.886314;

  //yellow and blue
  //flower.R = 0.815274, flower.G = 0.984891, flower.B = 0.118352;

  newTree(&object_list, cameraMat, NULL, 0, 4, PI/2, 4);


  o = newSphere(.05,.95,.95,.75, 1, 1, 1,1,1,6);
  Scale(o, 3, 3, 3);
  matMult(cameraMat, o->T);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  //insertObject(o,&object_list);


  p = path(iter_val/num_frames);

//   //red firework
//   o = newSphere(1, 0, 0, 0, 1, 0, 0,1,1,6);
//   Scale(o, 0.5, 0.5, 0.5);
//   Translate(o, p.px, p.py, p.pz);
//   Translate(o, 5, 0, 0);

//   matMult(cameraMat, o->T);
//   invert(&o->T[0][0],&o->Tinv[0][0]);
//   insertObject(o,&object_list);

//   pl.px = p.px + 5, pl.py = p.py - 0.6, pl.pz = p.pz;
//   matVecMult(cameraMat, &pl);
//   l=newPLS(&pl, f*1, 0, 0);
//   //insertPLS(l,&light_list);

//   //green firework
//   o = newSphere(1, 0, 0, 0, 0, 1, 0,1,1,6);
//   Scale(o, 0.5, 0.5, 0.5);
//   Translate(o, p.px, p.py, p.pz);
//   Translate(o, -2.5, 0, 4.33);

//   matMult(cameraMat, o->T);
//   invert(&o->T[0][0],&o->Tinv[0][0]);
//   insertObject(o,&object_list);

//   pl.px = p.px -2.5, pl.py = p.py - 0.6, pl.pz = p.pz + 4.33;
//   matVecMult(cameraMat, &pl);
//   l=newPLS(&pl, 0, f*1, 0);
  //insertPLS(l,&light_list);

  double x, y, z, ff = 200;
  double b;
  double sp = 12;
  for (int t = 0; t < ff; t++){
    x = drand48()*sp - sp/2;
    y = drand48()*2 + 7;
    z = drand48()*sp - sp/2;
    b = abs(sin(2*PI*iter_val/40 + drand48()*PI));

    p.px = sin(iter_val*drand48()/40);
    p.py = sin(iter_val*drand48()/30);
    p.pz *= drand48();
    o = newSphere(1, 0, 0, 0, 1*b, 1*b, 1*b,1,1,6);
    Scale(o, 0.05, 0.05, 0.05);
    Translate(o, p.px, p.py, p.pz);
    Translate(o, x, y, z);
    matMult(cameraMat, o->T);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    pl.px = p.px + x, pl.py = p.py + y - 0.6, pl.pz = p.pz + z;
    matVecMult(cameraMat, &pl);
    l=newPLS(&pl, 1/ff, 1/ff, 1/ff);
    insertPLS(l,&light_list);
  }

  //blue firework
//   o = newSphere(1, 0, 0, 0, 0, 0, 1,1,1,6);
//   Scale(o, 0.5, 0.5, 0.5);
//   Translate(o, p.px, p.py, p.pz);
//   Translate(o, -2.5, 0, -4.33);

//   matMult(cameraMat, o->T);
//   invert(&o->T[0][0],&o->Tinv[0][0]);
//   insertObject(o,&object_list);

//   pl.px = p.px -2.5, pl.py = p.py - 0.6, pl.pz = p.pz - 4.33;
//   matVecMult(cameraMat, &pl);
//   l=newPLS(&pl, 0, 0, f*1);
//  insertPLS(l,&light_list);

//   newFlower(&object_list, &flower, cameraMat);

  o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
  loadTexture(o, file, 1, &t_list);
  o->label = 'P';
  double s = 2;
  Scale(o,s*11,s*11,s*11);
  RotateZ(o,PI/4);
  RotateX(o,PI/2);

  matMult(cameraMat, o->T);

  invert(&o->T[0][0],&o->Tinv[0][0]);
  insertObject(o,&object_list);

 // Insert a single point light source. We set up its position as a point structure, and specify its
 // colour in terms of RGB (in [0,1]).
//  p.px=0;
//  p.py=25.5;
//  p.pz=-3.5;
//  p.pw=1;
//  l=newPLS(&p,(1-f)*1,(1-f)*1,(1-f)*1);
//  insertPLS(l,&light_list);
