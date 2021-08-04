
Object *o;
point p;
PointLS *l;
matrix t;

//living room scene
const std::string checkers_texture = "scenes/checkers.ppm";
const std::string diamond_texture = "scenes/cube/diamond.ppm";
const std::string lava_texture = "scenes/cube/lava.ppm";
const std::string world_mask_texture = "scenes/sphere/world-mask.ppm";

//const char *room_mesh = "scenes/living_room/tea_pot_N_992.obj";
const std::string room_mesh = "scenes/living_room/living_room.obj";
const std::string cube_mesh = "scenes/cube/cube.obj";
const std::string sphere_mesh = "scenes/sphere/sphere-cylcoords-4k.obj";

scene->cam_pos = point(0, 0, -10);
scene->cam_pos = RotY(-scene->frame * PI/120) * scene->cam_pos;
scene->cam_gaze_point = point(0, 0, 0);
scene->cam_gaze = scene->cam_gaze_point - scene->cam_pos;
//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -3;


scene->meshFactory.setDefaultColor(74 / 255.0, 255 / 255.0, 249 / 255.0);

scene->meshFactory.setTransform(Tr(0,-5,0) * Sc(40));
scene->meshFactory.loadMeshFile(room_mesh);

/*
scene->meshFactory.setDefaultColor(255 / 255.0, 74 / 255.0, 249 / 255.0);
t = Sc(1);
t *= Tr(2,0, 1);
scene->meshFactory.setTransform(t);
scene->meshFactory.loadMeshFile(cube_mesh);
*/

/*
scene->meshFactory.setDefaultColor(255 / 255.0, 249 / 255.0, 74 / 255.0);
t = Sc(1);
t *= RotX(scene->frame * PI / 120);
t *= RotY(scene->frame * PI / 120);
t *= Tr(-2,0,1);
scene->meshFactory.setTransform(t);
scene->meshFactory.loadMeshFile(sphere_mesh);
*/

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
//scene->insertObject(o);

p.x = 0;
p.y = 0;
p.z = 0;
p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));
