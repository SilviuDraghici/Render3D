struct object3D *o;
struct pointLS *l;
struct point p;

// Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)

o = newSphere(.05, .95, .35, .35, 1, .25, .25, 1, 1, 6);  // Initialize a sphere
o->T = Tr(2, 2.5, 1.5) * RotZ(PI / 4) * Sc(1.5, .75, .75);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);  // Compute the inverse transform * DON'T FORGET TO DO THIS! *
insertObject(o, &object_list);

o = newSphere(.05, .95, .95, .75, .75, .95, .55, 1, 1, 6);
o->T = Tr(-2.2, 1.75, 1.35) * RotZ(-PI / 1.5) * Sc(.95, 1.65, .65);
invert(&o->T.T[0][0], &o->Tinv.T[0][0]);
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .05, .55, .8, .75, 1, 1, 2);
//can still specify transforms line by line if I want
o->T = Sc(11);
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
