
Object *o;
point p;
PointLS *l;

struct textureNode *t_list = NULL;

const char *mesh_file = "scenes/wineglass_13k.obj";

// Cornell box

scene->cam_pos = point(0, 0, -15);
//scene->cam_gaze_point = point(0, 0, 0);
//scene->cam_gaze = cam_gaze_point - cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;


o = new Mesh(1.0, 1.0, 1.0);
((Mesh *)o)->setMesh(mesh_file);
strcpy(o->label, "Mesh test");
//o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->set_pathTrace_properties(0.0, 0.0, 1.0);
o->r_index = 1.54;
//o->T *= RotY(scene->frame * PI / 120);
//o->T *= RotX(scene->frame * PI / 120);
o->T *= Sc(13.0);
o->T *= Tr(0, -4, 2);
o->invert_and_bound();
scene->insertObject(o);

// Left
o = new Plane(.75, .25, .25);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Left Wall");
o->T *= RotY(PI / 2);
o->T *= Sc(25);
o->T *= Tr(-10, 0, 5);
o->invert_and_bound();
scene->insertObject(o);

// Right
o = new Plane(.25, .25, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Right Wall");
o->T *= RotY(PI / 2);
o->T *= Sc(25);
o->T *= Tr(10, 0, 5);
o->invert_and_bound();
scene->insertObject(o);

// Back
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Back Wall");
//o->T *= RotateZ(o, PI/4);
o->T *= Sc(10);
o->T *= Tr(0, 0, 15);
o->invert_and_bound();
scene->insertObject(o);

// Bottom
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Bottom Wall");
o->T *= RotX(PI / 2);
o->T *= Sc(25);
o->T *= Tr(0, -10, 5);
o->invert_and_bound();
scene->insertObject(o);

// Top
o = new Plane(.75, .75, .75);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.4;
strcpy(o->label, "Top Wall");
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
strcpy(o->label, "Top Light");
o->T *= Sc(0.5, 2.5, 1);
o->T *= RotX(PI / 2);
o->T *= Tr(0, 9.995, 5);
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
