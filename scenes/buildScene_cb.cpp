
Object *o;
point p;
PointLS *l;

struct textureNode *t_list = NULL;

// Cornell box

scene->cam_pos = point(0, 0, -15);
//scene->cam_gaze_point = point(0, 0, 0);
//scene->cam_gaze = cam_gaze_point - cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;

// Two spheres scene
// Refract
o = new Sphere(.99, .99, .99);
o->name = "Left Refract Shpere";
o->set_pathTrace_properties(0.0, 0.0, 1.0);
o->r_index = 1.54;
o->T *= RotY(PI);
o->T *= Sc(3.75);
o->T *= Tr(-5, -4.0, 4.5);
o->invert_and_bound();
scene->insertObject(o);

// Reflect
o = new Sphere(.99, .99, .99);
o->name = "Right Reflect Sphere";
o->set_pathTrace_properties(0.0, 1.0, 0.0);
o->refl_sig = 0.05;
o->r_index = 2.47;
o->T *= Sc(3.75);
o->T *= Tr(4, -3.75, 6.5);
o->invert_and_bound();
scene->insertObject(o);


// Left
o = new Plane(.75, .25, .25);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
o->name = "Left Wall";
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
o->T *= RotX(-PI / 2);
o->T *= Sc(25);
o->T *= Tr(0, 10, 5);
o->invert_and_bound();
scene->insertObject(o);

// Planar light source at top
o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
o->name = "Top Light";
o->T *= Sc(0.5, 2.5, 1);
o->T *= RotX(PI / 2);
o->T *= Tr(0, 9.995, 5);
o->isLightSource = 1;
o->pt.LSweight *= 0.5 * 2.5 * 1;  // <- scale weight by scale
o->invert_and_bound();
scene->insertObject(o);

o = new Plane(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
o->name = "Bottom Light";
o->T *= Sc(0.5, 2.5, 1);
o->T *= RotX(-PI / 2);
o->T *= Tr(0, -9.995, 5);
o->isLightSource = 1;
o->pt.LSweight *= 0.5 * 2.5 * 1;  // <- scale weight by scale
o->invert_and_bound();
scene->insertObject(o);

p.x = 0;
p.y = 9.9;
p.z = 5;
p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));
