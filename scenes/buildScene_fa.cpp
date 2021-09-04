
Object *o;
point p;
PointLS *l;
matrix t;

//fancy appartment scene

const std::string room_mesh = "scenes/fancy_apartment/Svetaine.obj";

//default position for rotating the camera around the room
// scene->cam_pos = point(0, 0, -10);
// scene->cam_pos = RotY(-scene->frame * PI/120) * scene->cam_pos;
// scene->cam_gaze_point = point(0, 0, 0);

/* kitchen camera angle */
scene->cam_pos = point(5, -7.5, -11);
//scene->cam_pos = Tr(-2,0,0) * RotY(-15 * PI/120) * scene->cam_pos;
scene->cam_gaze = RotY(15 * PI/120) * point(0, 0, 1);


/* living room angle 
scene->cam_pos = point(22, -5.5, 15);
//scene->cam_pos = Tr(-2,0,0) * RotY(-15 * PI/120) * scene->cam_pos;
scene->cam_gaze = RotY(30 * PI/120) * point(0, 0, -1);
*/

/* upstairs camera angle 
scene->cam_pos = point(5, 7, 15);
//scene->cam_pos = Tr(-2,0,0) * RotY(-15 * PI/120) * scene->cam_pos;
scene->cam_gaze = RotY(0 * PI/120) * point(-1, 0, 0);
*/

//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -5;

scene->exposure = 1;

scene->meshFactory.setDefaultColor(74 / 255.0, 255 / 255.0, 249 / 255.0);

scene->meshFactory.setTransform(Tr(0,0,0) * RotX(PI/2) * RotZ(PI/2) * Sc(80));
scene->meshFactory.loadMeshFile(room_mesh);

// Planar light source outside window
o = new Plane(40, 40, 40);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
//setTexture(o, world_mask_texture, 1, scene->texture_list);
o->refl_sig = 0.0;
o->r_index = 1.54;
o->name = "Window Light";
o->T *= Sc(9.6, 5.4, 1);
o->T *= RotX(PI);
o->T *= Tr(-4.6, 0.8, 12);
o->isLightSource = 1;
o->pt.LSweight *= 6 * 1 * 1;  // <- scale weight by scale
o->invert_and_bound();
//scene->insertObject(o);


// Planar light source at top
o = new Plane(40.0, 40.0, 40.0);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
o->refl_sig = 0.0;
o->r_index = 1.54;
o->name = "Top Light";
o->T *= Sc(0.5, 2.5, 1);
o->T *= RotX(PI / 2);
o->T *= Tr(0, 9.995, 0);
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
