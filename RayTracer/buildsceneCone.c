
 struct object3D *o;
 struct pointLS *l;
 struct textureNode *t_list = NULL;
 struct point p, pl;
 p.pw = 1, pl.pw = 1;

 double f = 0.1 ;

 const char* file = "grass.ppm";

 double angle;
 if(iter_val != -1){
	 angle = iter_val*2*PI/num_frames;
 }else{
	 angle = 0;
 }

 double cameraMat[4][4] = {{1.0, 0.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0, 0.0},
                    {0.0, 0.0, 1.0, 0.0},
                    {0.0, 0.0, 0.0, 1.0}};


// RotateXMat(cameraMat, angle);
 RotateYMat(cameraMat, angle);

 TranslateMat(cameraMat, 0, 0, 3);
  o = newCone(.05,.95,.95,.75,.75,.95, 0.25,1,1,6);
  Scale(o, 2, 2, 4);


  matMult(cameraMat, o->T);
  invert(&o->T[0][0],&o->Tinv[0][0]);
  insertObject(o,&object_list);
 p.px=0;
 p.py=0;
 p.pz=0;
 p.pw=1;
 l=newPLS(&p,(1-f)*1,(1-f)*1,(1-f)*1);
 insertPLS(l,&light_list);

