
Object *o;
struct point p;
struct pointLS *l;

struct textureNode *t_list = NULL;

const char *file = "scenes/marble.ppm";
const char *jupfile = "scenes/jupiter.ppm";
const char *nfile = "scenes/stone_normal.ppm";
const char *alphafile = "scenes/earthalpha.pgm";

const char *fence = "fence.ppm";
const char *fencenormal = "fencenormal.ppm";
const char *fencealpha = "fencealpha.pgm";

// Cornell box

scene->cam_pos = point(0, 0, -15);
//scene->cam_gaze_point = point(0, 0, 0);
//scene->cam_gaze = cam_gaze_point - cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;

// Left
o = new Plane(.75, .25, .25);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Left Wall");
loadTexture(o, nfile, 2, &t_list);
o->T *= RotY(PI / 2);
o->T *= Sc(25);
o->T *= Tr(-10, 0, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Right
o = new Plane(.25, .25, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Right Wall");
o->T *= RotY(PI / 2);
o->T *= Sc(25);
o->T *= Tr(10, 0, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Back
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Back Wall");
loadTexture(o, file, 1, &t_list);
//o->T *= RotateZ(o, PI/4);
o->T *= Sc(10);
o->T *= Tr(0, 0, 15);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Bottom
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Bottom Wall");
o->T *= RotX(PI / 2);
o->T *= Sc(25);
o->T *= Tr(0, -10, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Top
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Top Wall");
o->T *= RotX(PI / 2);
o->T *= Sc(25);
o->T *= Tr(0, 10, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Two spheres scene
// Refract
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(0.0, 0.0, 1.0);
o->r_index = 1.54;
o->T *= RotY(PI);
o->T *= Sc(3.75);
o->T *= Tr(-5, -4.0, 4.5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Refract
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(0.0, 0.0, 1.0);
o->r_index = 1.54;
loadTexture(o, alphafile, 3, &t_list);  //alpha map
o->T *= RotY(PI);
o->T *= Sc(3.75);
o->T *= Tr(-5, 6.0, 8.5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
//insertObject(o, &(scene->object_list));

// Reflect
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(0.0, 1.0, 0.0);
o->refl_sig = .05;
o->r_index = 2.47;
strcpy(o->label, "Right Sphere");
//loadTexture(o, nfile, 2, &t_list);
o->T *= Sc(3.75);
o->T *= Tr(4, -3.75, 6.5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Jupiter
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 2.47;
strcpy(o->label, "Jupiter");
loadTexture(o, jupfile, 1, &t_list);
o->T *= Sc(2);
o->T *= Tr(7, 7, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &(scene->object_list));

// Planar light source at top
o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
strcpy(o->label, "Top Light");
o->T *= Sc(.5, 2.5, 1);
o->T *= RotX(PI / 2);
o->T *= Tr(0, 9.995, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
o->isLightSource = 1;
insertObject(o, &(scene->object_list));

p.x = 0;
p.y = 9.99;
p.z = 5;
p.w = 1;
l = newPLS(&p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));

// Planar light source at bottom
o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
strcpy(o->label, "bottom light");
o->T *= Sc(.5, 2.5, 1);
o->T *= RotX(-PI / 2);
o->T *= Tr(0, -9.995, 5);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
o->isLightSource = 1;
//insertObject(o, &(scene->object_list));
