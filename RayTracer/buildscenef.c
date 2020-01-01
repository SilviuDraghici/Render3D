struct object3D *o;
struct pointLS *l;
struct point p;

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
double angle = PI;
transform cameraMat = Tr(0,-4,5)*Sc(0.5) * RotY(angle);
newTree(&object_list, cameraMat, NULL, 0, 3, PI / 2, 3);

o = newSphere(.05, .95, .95, .75, .75, .95, .55, 1, 1, 6);
o->T = Tr(-2.2, 1.75, 1.35) * RotZ(-PI / 1.5) * Sc(.95, 1.65, .65);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .05, .55, .8, .75, 1, 1, 2);
//can still specify transforms line by line if I want
o->T = Sc(11, 11, 11);
o->T *= RotZ(PI / 4);
o->T *= RotX(PI / 2);
o->T *= Tr(0, -4, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &object_list);

// Insert a single point light source. We set up its position as a point structure, and specify its
// colour in terms of RGB (in [0,1]).
p.x = 0;
p.y = 25.5;
p.z = -3.5;
p.w = 1;
l = newPLS(&p, .95, .95, .95);
insertPLS(l, &light_list);
