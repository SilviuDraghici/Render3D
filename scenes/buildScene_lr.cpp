
Object *o;
point p;
PointLS *l;

struct textureNode *t_list = NULL;

// Cornell box with a living room inside

//const char *room_mesh = "scenes/living_room/tea_pot_N_992.obj";
const char *room_mesh = "scenes/living_room/living_room.obj";
const char *cube_mesh = "scenes/cube/cube.obj";
const char *sphere_mesh = "scenes/sphere/sphere-cylcoords-4k.obj";

scene->cam_pos = point(0, 0, -10);
scene->cam_pos = RotY(scene->frame * PI/120) * scene->cam_pos;
scene->cam_gaze_point = point(0, 0, 0);
scene->cam_gaze = scene->cam_gaze_point - scene->cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;


o = new Mesh(74 / 255.0, 255 / 255.0, 249 / 255.0);
((Mesh *)o)->setMesh(room_mesh);
strcpy(o->label, "room mesh");
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->T *= Sc(40.0);
//o->T *= RotX(scene->frame * PI / 120);
//o->T *= RotY(scene->frame * PI / 120);
o->T *= Tr(0,-5,0);
//o->T *= RotX(-PI / 2);
//o->T *= Tr(4, 0, 0);
o->invert_and_bound();
scene->insertObject(o);

o = new Mesh(255 / 255.0, 74 / 255.0, 249 / 255.0);
((Mesh *)o)->setMesh(cube_mesh);
strcpy(o->label, "cube mesh");
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->T *= Sc(1);
//o->T *= RotX(scene->frame * PI / 120);
//o->T *= RotY(scene->frame * PI / 120);
o->T *= Tr(2,0, 1);
//o->T *= RotX(-PI / 2);
//o->T *= Tr(4, 0, 0);
o->invert_and_bound();
scene->insertObject(o);


o = new Mesh(255 / 255.0, 249 / 255.0, 74 / 255.0);
((Mesh *)o)->setMesh(cube_mesh);
strcpy(o->label, "sphere mesh");
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->r_index = 1.54;
o->T *= Sc(1);
//o->T *= RotX(scene->frame * PI / 120);
//o->T *= RotY(scene->frame * PI / 120);
o->T *= Tr(-2,0,1);
//o->T *= RotX(-PI / 2);
//o->T *= Tr(4, 0, 0);
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
//scene->insertObject(o);

p.x = 0;
p.y = 0;
p.z = 0;
p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));
