#include "camera.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ray.h"
#include "utils.h"

struct view *setupView(struct point *e, struct point *g, struct point *up, double f, double wl, double wt, double wsize) {
    /*
    This function sets up the camera axes and viewing direction as discussed in the
    lecture notes.
    e - Camera center
    g - Gaze direction
    up - Up vector
    fov - Field of view in degrees
    f - focal length
    */
    struct view *c;
    struct point u, v;


    // Allocate space for the camera structure
    c = (struct view *)calloc(1, sizeof(struct view));
    if (c == NULL) {
        fprintf(stderr, "Out of memory setting up camera model!\n");
        return (NULL);
    }

    // Set up camera center and axes
    c->e.x = e->x;  // Copy camera center location, note we must make sure
    c->e.y = e->y;  // the camera center provided to this function has pw=1
    c->e.z = e->z;
    c->e.w = 1;

    // Set up w vector (camera's Z axis). w=-g/||g||
    c->w.x = -g->x;
    c->w.y = -g->y;
    c->w.z = -g->z;
    c->w.w = 1;
    normalize(&c->w);

    // Set up the horizontal direction, which must be perpenticular to w and up
    u = cross(&c->w, up);
    normalize(&u);
    c->u.x = u.x;
    c->u.y = u.y;
    c->u.z = u.z;
    c->u.w = 1;

    // Set up the remaining direction, v=(u x w)  - Mind the signs
    v = cross(&c->u, &c->w);
    normalize(&v);
    c->v.x = v.x;
    c->v.y = v.y;
    c->v.z = v.z;
    c->v.w = 1;

    // Copy focal length and window size parameters
    c->f = f;
    c->wl = wl;
    c->wt = wt;
    c->wsize = wsize;

    // Set up coordinate conversion matrices
    // Camera2World matrix (M_cw in the notes)
    // Mind the indexing convention [row][col]
    c->C2W.T[0][0] = c->u.x;
    c->C2W.T[1][0] = c->u.y;
    c->C2W.T[2][0] = c->u.z;
    c->C2W.T[3][0] = 0;

    c->C2W.T[0][1] = c->v.x;
    c->C2W.T[1][1] = c->v.y;
    c->C2W.T[2][1] = c->v.z;
    c->C2W.T[3][1] = 0;

    c->C2W.T[0][2] = c->w.x;
    c->C2W.T[1][2] = c->w.y;
    c->C2W.T[2][2] = c->w.z;
    c->C2W.T[3][2] = 0;

    c->C2W.T[0][3] = c->e.x;
    c->C2W.T[1][3] = c->e.y;
    c->C2W.T[2][3] = c->e.z;
    c->C2W.T[3][3] = 1;

    // World2Camera matrix (M_wc in the notes)
    // Mind the indexing convention [row][col]
    c->W2C.T[0][0] = c->u.x;
    c->W2C.T[1][0] = c->v.x;
    c->W2C.T[2][0] = c->w.x;
    c->W2C.T[3][0] = 0;

    c->W2C.T[0][1] = c->u.y;
    c->W2C.T[1][1] = c->v.y;
    c->W2C.T[2][1] = c->w.y;
    c->W2C.T[3][1] = 0;

    c->W2C.T[0][2] = c->u.z;
    c->W2C.T[1][2] = c->v.z;
    c->W2C.T[2][2] = c->w.z;
    c->W2C.T[3][2] = 0;

    c->W2C.T[0][3] = -dot(&c->u, &c->e);
    c->W2C.T[1][3] = -dot(&c->v, &c->e);
    c->W2C.T[2][3] = -dot(&c->w, &c->e);
    c->W2C.T[3][3] = 1;

    return (c);
}

void setPixelStep(Scene *scene, struct view *cam, double sx, double sy) {
    scene->du = cam->wsize / (sx - 1);  // du and dv. In the notes in terms of wl and wr, wt and wb,
    scene->dv = -cam->wsize / (sx - 1);
}

void getRayFromPixel(Scene *scene, struct Ray *ray, struct view *cam, double i, double j) {
    struct point pc;
    pc.x = cam->wl + i * scene->du;
    pc.y = cam->wt + j * scene->dv;
    pc.z = cam->f;
    pc.w = 1;

    // Convert image plane sample coordinates to world coordinates
    pc = cam->C2W * pc;

    // Now compute the ray direction
    ray->d = pc - cam->e;
    normalize(&ray->d);

    ray->p0 = pc;
}
