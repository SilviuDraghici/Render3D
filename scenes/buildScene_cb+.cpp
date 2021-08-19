
Object *o;
point p;
PointLS *l;

struct textureNode *t_list = NULL;

const char *file = "scenes/marble.png";
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
o->name = "Left Wall";
setTexture(o, nfile, 2, scene->texture_list);
o->T *= RotY(PI / 2);
o->T *= Sc(25);
o->T *= Tr(-10, 0, 5);
o->invert_and_bound();
scene->insertObject(o);

// Right
o = new Plane(.25, .25, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
o->name = "Right Wall";
o->T *= RotY(PI / 2);
o->T *= Sc(25);
o->T *= Tr(10, 0, 5);
o->invert_and_bound();
scene->insertObject(o);

// Back
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
o->name = "Back Wall";
setTexture(o, file, 1, scene->texture_list);

//o->T *= RotateZ(o, PI/4);
o->T *= Sc(10);
o->T *= Tr(0, 0, 15);
o->invert_and_bound();
scene->insertObject(o);

// Bottom
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
o->name = "Bottom Wall";
o->T *= RotX(PI / 2);
o->T *= Sc(25);
o->T *= Tr(0, -10, 5);
o->invert_and_bound();
scene->insertObject(o);

// Top
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
o->name = "Top Wall";
o->T *= RotX(PI / 2);
o->T *= Sc(25);
o->T *= Tr(0, 10, 5);
o->invert_and_bound();
scene->insertObject(o);

// Two spheres scene
// Refract
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(0.0, 0.0, 1.0);
o->r_index = 1.54;
o->T *= RotY(PI);
o->T *= Sc(3.75);
o->T *= Tr(-5, -4.0, 4.5);
o->invert_and_bound();
scene->insertObject(o);

// Refract
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(0.0, 0.0, 1.0);
o->r_index = 1.54;
setTexture(o, alphafile, 3, scene->texture_list);  //alpha map
o->T *= RotY(PI);
o->T *= Sc(3.75);
o->T *= Tr(-5, 6.0, 8.5);
o->invert_and_bound();
//scene->insertObject(o);

// Reflect
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(0.0, 1.0, 0.0);
o->refl_sig = .05;
o->r_index = 2.47;
o->name = "Right Sphere";
//setTexture(o, nfile, 2, scene->texture_list);
o->T *= Sc(3.75);
o->T *= Tr(4, -3.75, 6.5);
o->invert_and_bound();
scene->insertObject(o);

// Jupiter
o = new Sphere(.99, .99, .99);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 2.47;
o->name = "Jupiter";
setTexture(o, jupfile, 1, scene->texture_list);
o->T *= Sc(2);
o->T *= Tr(7, 7, 5);
o->invert_and_bound();
scene->insertObject(o);

// Planar light source at top
o = new Plane(20.0, 20.0, 20.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->name = "Top Light";
o->T *= Sc(.5, 2.5, 1);
o->T *= RotX(PI / 2);
o->T *= Tr(0, 9.995, 5);
o->isLightSource = 1;
o->invert_and_bound();
scene->insertObject(o);

p.x = 0;
p.y = 9.99;
p.z = 5;
p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));

// Planar light source at bottom
o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
o->name = "bottom light";
o->T *= Sc(.5, 2.5, 1);
o->T *= RotX(-PI / 2);
o->T *= Tr(0, -9.995, 5);
o->isLightSource = 1;
o->invert_and_bound();
//scene->insertObject(o);
