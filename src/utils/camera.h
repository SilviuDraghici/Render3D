#ifndef CAMERA_H
#define CAMERA_H

#include "utils.h"

/*
   The structure below is used to hold camera parameters. You will need
   to write code to initialize the camera position and orientation.
*/
struct view {
   struct point e;       // Location of the camera center
   struct point u;       // u vector
   struct point v;       // v vector
   struct point w;       // w vector
   double f;             // Focal length
   double wl;            // Left edge in camera coordinates
   double wt;            // Top edge in camera coordinates
   double wsize;         // Window size in distance units (not pixels!)
   struct matrix W2C; // World2Camera conversion matrix
   struct matrix C2W; // Camera2World conversion matrix
};

extern struct point cam_pos;
extern struct point cam_up;
extern struct point cam_gaze;
extern struct point cam_gaze_point;
extern double cam_focal; // should be negative
extern double du,dv;

struct view *setupView(struct point *e, struct point *g, struct point *up, double f, double wl, double wt, double wsize);
void setPixelStep(struct view * cam, double sx, double sy);
void getRayFromPixel(struct ray *ray, struct view* cam, double i, double j);
#endif