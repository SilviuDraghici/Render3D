struct object3D *o;
struct pointLS *l;
struct point3D p;

// Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)

o = newSphere(.05, .95, .35, .35, 1, .25, .25, 1, 1, 6); // Initialize a sphere
Scale(o, 1.5, .75, .75);                                 // Apply a few transforms (Translate * Rotate * Scale)
RotateZ(o, PI / 4);
Translate(o, 2.0, 2.5, 1.5);
invert(&o->T[0][0], &o->Tinv[0][0]); // Compute the inverse transform * DON'T FORGET TO DO THIS! *

// If needed, this is how you load a texture map
// loadTexture(o,"./Texture/mosaic2.ppm",1,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
// texture gets added to the texture list, and a
// pointer to it is stored within this object in the
// corresponding place. The '1' indicates this image
// will be used as a texture map. Use '2' to load
// an image as a normal map, and '3' to load an
// alpha map. Texture and normal maps are RGB .ppm
// files, alpha maps are grayscale .pgm files.
// * DO NOT * try to free image data loaded in this
// way, the cleanup function already provided will do
// this at the end.

insertObject(o, &object_list); // <-- If you don't insert the object into the object list,
                               //     nothing happens! your object won't be rendered.

// That's it for defining a single sphere... let's add a couple more objects
o = newSphere(.05, .95, .95, .75, .75, .95, .55, 1, 1, 6);
Scale(o, .95, 1.65, .65);
RotateZ(o, -PI / 1.5);
Translate(o, -2.2, 1.75, 1.35);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newPlane(.05, .75, .05, .05, .55, .8, .75, 1, 1, 2);
Scale(o, 11, 11, 11);
RotateZ(o, PI / 4);
RotateX(o, PI / 2);
Translate(o, 0, -4, 5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

// Insert a single point light source. We set up its position as a point structure, and specify its
// colour in terms of RGB (in [0,1]).
p.px = 0;
p.py = 25.5;
p.pz = -3.5;
p.pw = 1;
l = newPLS(&p, .95, .95, .95);
insertPLS(l, &light_list);
