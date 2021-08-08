#ifndef IMAGEPROCESSOR_H
#define IMAGEPROCESSOR_H
/* The structure below is used to hold a single RGB image. 'rgbdata' points */
/* to an array of size sx*sy*[num layers] but note that your code must know */
/* the data type for the array. Within the raytracer, there will be double  */
/* floating point images, and also unsigned char images, and there will be  */
/* 3-layer images (your raytraced scene, texture and normal maps) as well   */
/* as 1-layer images (alpha maps)					      */
struct image {
    void *rgbdata;
    int sx;
    int sy;
};

extern struct image *outImage;
extern char output_name[1024];  // Name of the output file for the raytraced .ppm image

struct image *readPPMimage(const char *filename);
struct image *readPGMimage(const char *filename);
struct image *newImage(int size_x, int size_y, int pixel_size);


bool PNGImageOutput(image *im, const char *filename);
void PPMImageOutput(image *im, const char *filename);

void dataOutput(double *im, int sx, char *name);
void deleteImage(struct image *im);
#endif