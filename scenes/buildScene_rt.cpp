Object *o;
point p;
PointLS *l;

scene->cam_pos = point(0, 0, -1);

//scene->cam_gaze_point = point(0, 0, 0);
//scene->cam_gaze = scene->cam_gaze_point - scene->cam_pos;
//scene->cam_up = point(0, 1, 0);
//scene->cam_focal = -1;

o = new Sphere(1, .25, .25);  // Initialize a sphere
o->set_rayTrace_properties(.05, .95, .35, .35, 1, 6);
o->T = Tr(2, 2.5, 1.5) * RotZ(PI / 4) * Sc(1.5, .75, .75);
o->invert_and_bound();     // Compute the inverse transform * DON'T FORGET TO DO THIS! *
scene->insertObject(o);

o = new Sphere(.75, .95, .55);
o->set_rayTrace_properties(.05, .95, .95, .75, 1, 6);
o->T = Tr(-2.2, 1.75, 1.35) * RotZ(-PI / 1.5) * Sc(.95, 1.65, .65);
o->invert_and_bound();
scene->insertObject(o);

o = new Plane(.55, .8, .75);
o->set_rayTrace_properties(.05, .75, .05, .05, 1, 2);
//can still specify transforms line by line if I want
o->T *= Sc(11);
o->T *= RotZ(PI / 4);
o->T *= RotX(PI / 2);
o->T *= Tr(0, -4, 5);
o->invert_and_bound();
scene->insertObject(o);

// Insert a single point light source. We set up its position as a point structure, and specify its
// colour in terms of RGB (in [0,1]).
p.x = 0;
p.y = 25.5;
p.z = -3.5;
p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &scene->rt_point_light_list);

//sphere in case of path tracing
o = new Sphere(1.0, 1.0, 1.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
o->name = "Top Light";
o->T *= Tr(0, 30, -3.5);
o->isLightSource = 1;
o->invert_and_bound();
scene->insertObject(o);
