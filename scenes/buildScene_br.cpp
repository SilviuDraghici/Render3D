
Object *o;
point p;
PointLS *l;
matrix t;

//breakfast room
const std::string outside_window_light = "scenes/backyard.png";

const std::string room_mesh = "scenes/breakfast_room/breakfast_room.obj";

//default position for rotating the camera around the room
scene->cam_pos = point(0, 3, -29);
scene->cam_pos = scene->cam_pos;
scene->cam_gaze_point = point(0, 0, 0);
scene->cam_gaze = point(0, 0, 1);

//cam pos similar to refernce image
// scene->cam_pos = point(0, 0, -16);
// scene->cam_pos = Tr(-2,0,0) * RotY(-15 * PI/120) * scene->cam_pos;
// scene->cam_gaze_point = point(0, -4, 0);
// scene->cam_gaze = scene->cam_gaze_point - scene->cam_pos;
// scene->cam_gaze = {-6.8883, -4, 16.6298};
//std::cout << "cam_gaze: " << scene->cam_gaze << std::endl;

//scene->cam_up = point(0, 1, 0);
scene->cam_focal = -6;

scene->exposure = 1;

scene->meshFactory.setDefaultColor(74 / 255.0, 255 / 255.0, 249 / 255.0);

scene->meshFactory.setTransform(Tr(0,0,0) * Sc(40));
scene->meshFactory.loadMeshFile(room_mesh);

// Planar light source outside window
o = new Plane(20, 20, 20);
o->set_pathTrace_properties(1.0, 0.0, 0.0);
//setTexture(o, outside_window_light, 1, scene->texture_list);
o->refl_sig = 0.0;
o->r_index = 1.54;
o->name = "Window Light";
o->T *= Sc(15, 10, 1);
o->T *= RotX(3 * PI/4);
o->T *= RotY(PI/2);
o->T *= Tr(35, 30, 0);
o->isLightSource = 1;
o->pt.LSweight *= 6 * 1 * 1;  // <- scale weight by scale
o->invert_and_bound();
scene->insertObject(o);

p.x = 0;
p.y = 5;
p.z = 0;
p.w = 1;
l = new PointLS(p, .95, .95, .95);
insertPLS(l, &(scene->rt_point_light_list));