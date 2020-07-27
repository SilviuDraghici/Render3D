#include "imageProcessor.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct image *outImage;
char output_name[1024];  // Name of the output file for the raytraced .ppm image

struct image *readPPMimage(const char *filename) {
    // Reads an image from a .ppm file. A .ppm file is a very simple image representation
    // format with a text header followed by the binary RGB data at 24bits per pixel.
    // The header has the following form:
    //
    // P6
    // # One or more comment lines preceded by '#'
    // 340 200
    // 255
    //
    // The first line 'P6' is the .ppm format identifier, this is followed by one or more
    // lines with comments, typically used to inidicate which program generated the
    // .ppm file.
    // After the comments, a line with two integer values specifies the image resolution
    // as number of pixels in x and number of pixels in y.
    // The final line of the header stores the maximum value for pixels in the image,
    // usually 255.
    // After this last header line, binary data stores the RGB values for each pixel
    // in row-major order. Each pixel requires 3 bytes ordered R, G, and B.
    //
    // NOTE: Windows file handling is rather crotchetty. You may have to change the
    //       way this file is accessed if the images are being corrupted on read
    //       on Windows.
    //
    // readPPMdata converts the image colour information to floating point. This is so that
    // the texture mapping function doesn't have to do the conversion every time
    // it is asked to return the colour at a specific location.
    //

    FILE *f;
    struct image *im;
    char line[1024];
    int sizx, sizy;
    int i;
    unsigned char *tmp;
    double *fRGB;
    int tmpi;
    char *tmpc;

    im = (struct image *)calloc(1, sizeof(struct image));
    if (im != NULL) {
        im->rgbdata = NULL;
        f = fopen(filename, "rb+");
        if (f == NULL) {
            fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
            free(im);
            return (NULL);
        }
        tmpc = fgets(&line[0], 1000, f);
        if (strcmp(&line[0], "P6\n") != 0) {
            fprintf(stderr, "Wrong file format, not a .ppm file or header end-of-line characters missing\n");
            free(im);
            fclose(f);
            return (NULL);
        }
        // fprintf(stderr,"%s\n",line);
        // Skip over comments
        tmpc = fgets(&line[0], 511, f);
        while (line[0] == '#') {
            //fprintf(stderr,"%s",line);
            tmpc = fgets(&line[0], 511, f);
        }
        sscanf(&line[0], "%d %d\n", &sizx, &sizy);  // Read file size
        //fprintf(stderr,"nx=%d, ny=%d\n\n",sizx,sizy);
        im->sx = sizx;
        im->sy = sizy;

        tmpc = fgets(&line[0], 9, f);  // Read the remaining header line
        //fprintf(stderr,"%s\n",line);
        tmp = (unsigned char *)calloc(sizx * sizy * 3, sizeof(unsigned char));
        fRGB = (double *)calloc(sizx * sizy * 3, sizeof(double));
        if (tmp == NULL || fRGB == NULL) {
            fprintf(stderr, "Out of memory allocating space for image\n");
            free(im);
            fclose(f);
            return (NULL);
        }

        tmpi = fread(tmp, sizx * sizy * 3 * sizeof(unsigned char), 1, f);
        fclose(f);

        // Conversion to floating point
        for (i = 0; i < sizx * sizy * 3; i++)
            *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
        free(tmp);
        im->rgbdata = (void *)fRGB;

        return (im);
    }

    fprintf(stderr, "Unable to allocate memory for image structure\n");
    return (NULL);
}

struct image *readPGMimage(const char *filename) {
    // Just like readPPMimage() except it is used to load grayscale alpha maps. In
    // alpha maps, a value of 255 corresponds to alpha=1 (fully opaque) and 0
    // correspondst to alpha=0 (fully transparent).
    // A .pgm header of the following form is expected:
    //
    // P5
    // # One or more comment lines preceded by '#'
    // 340 200
    // 255
    //
    // readPGMdata converts the image grayscale data to double floating point in [0,1].

    FILE *f;
    struct image *im;
    char line[1024];
    int sizx, sizy;
    int i;
    unsigned char *tmp;
    double *fRGB;
    int tmpi;
    char *tmpc;

    im = (struct image *)calloc(1, sizeof(struct image));
    if (im != NULL) {
        im->rgbdata = NULL;
        f = fopen(filename, "rb+");
        if (f == NULL) {
            fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
            free(im);
            return (NULL);
        }
        tmpc = fgets(&line[0], 1000, f);
        if (strcmp(&line[0], "P5\n") != 0) {
            fprintf(stderr, "Wrong file format, not a .pgm file or header end-of-line characters missing\n");
            free(im);
            fclose(f);
            return (NULL);
        }
        // Skip over comments
        tmpc = fgets(&line[0], 511, f);
        while (line[0] == '#')
            tmpc = fgets(&line[0], 511, f);
        sscanf(&line[0], "%d %d\n", &sizx, &sizy);  // Read file size
        im->sx = sizx;
        im->sy = sizy;

        tmpc = fgets(&line[0], 9, f);  // Read the remaining header line
        tmp = (unsigned char *)calloc(sizx * sizy, sizeof(unsigned char));
        fRGB = (double *)calloc(sizx * sizy, sizeof(double));
        if (tmp == NULL || fRGB == NULL) {
            fprintf(stderr, "Out of memory allocating space for image\n");
            free(im);
            fclose(f);
            return (NULL);
        }

        tmpi = fread(tmp, sizx * sizy * sizeof(unsigned char), 1, f);
        fclose(f);

        // Conversion to double floating point
        for (i = 0; i < sizx * sizy; i++)
            *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
        free(tmp);
        im->rgbdata = (void *)fRGB;

        return (im);
    }

    fprintf(stderr, "Unable to allocate memory for image structure\n");
    return (NULL);
}

struct image *newImage(int size_x, int size_y, int pixel_size) {
    // Allocates and returns a new image with all zeros. Assumes 24 bit per pixel,
    // unsigned char array.
    struct image *im;

    im = (struct image *)calloc(1, sizeof(struct image));
    if (im != NULL) {
        im->rgbdata = NULL;
        im->sx = size_x;
        im->sy = size_y;
        im->rgbdata = (void *)calloc(size_x * size_y * 3, pixel_size);
        if (im->rgbdata != NULL)
            return (im);
    }
    fprintf(stderr, "Unable to allocate memory for new image\n");
    return (NULL);
}

void imageOutput(struct image *im, const char *filename) {
    // Writes out a .ppm file from the image data contained in 'im'.
    // Note that Windows typically doesn't know how to open .ppm
    // images. Use Gimp or any other seious image processing
    // software to display .ppm images.
    // Also, note that because of Windows file format management,
    // you may have to modify this file to get image output on
    // Windows machines to work properly.
    //
    // Assumes a 24 bit per pixel image stored as unsigned chars
    //
    printf("Saving Image\n");
    FILE *f;

    if (im != NULL)
        if (im->rgbdata != NULL) {
            f = fopen(filename, "wb+");
            if (f == NULL) {
                fprintf(stderr, "Unable to open file %s for output! No image written\n", filename);
                return;
            }
            fprintf(f, "P6\n");
            fprintf(f, "# Output from RayTracer.c\n");
            fprintf(f, "%d %d\n", im->sx, im->sy);
            fprintf(f, "255\n");
            fwrite((unsigned char *)im->rgbdata, im->sx * im->sy * 3 * sizeof(unsigned char), 1, f);
            fclose(f);
            return;
        }
    fprintf(stderr, "imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im) {
    // De-allocates memory reserved for the image stored in 'im'
    if (im != NULL) {
        if (im->rgbdata != NULL)
            free(im->rgbdata);
        free(im);
    }
}

void dataOutput(double *im, int sx, char *name) {
    FILE *f;
    double *imT;
    double HDRhist[1000];
    int i, j;
    double mx, mi, biw, pct;
    unsigned char *bits24;
    int name_len = strlen(name);
    char pfmname[1024];

    imT = (double *)calloc(sx * sx * 3, sizeof(double));
    memcpy(imT, im, sx * sx * 3 * sizeof(double));
    strcpy(&pfmname[0], name);
    if (pfmname[(name_len - 1) - 3] == '.') {
        pfmname[(name_len - 1) - 2] = 'p';
        pfmname[(name_len - 1) - 1] = 'f';
        pfmname[(name_len - 1) - 0] = 'm';
    } else {
        strcat(&pfmname[0], ".pfm");
    }
    // Output the floating point data so we can post-process externally
    f = fopen(pfmname, "w");
    fprintf(f, "PF\n");
    fprintf(f, "%d %d\n", sx, sx);
    fprintf(f, "%1.1f\n", -1.0);
    fwrite(imT, sx * sx * 3 * sizeof(double), 1, f);
    fclose(f);

    // Post processing HDR map - find reasonable cutoffs for normalization
    for (j = 0; j < 1000; j++)
        HDRhist[j] = 0;

    mi = 10e6;
    mx = -10e6;
    for (i = 0; i < sx * sx * 3; i++) {
        if (*(imT + i) < mi)
            mi = *(imT + i);
        if (*(imT + i) > mx)
            mx = *(imT + i);
    }

    for (i = 0; i < sx * sx * 3; i++) {
        *(imT + i) = *(imT + i) - mi;
        *(imT + i) = *(imT + i) / (mx - mi);
    }
    fprintf(stderr, "\nSaving Image\n");
    //fprintf(stderr, "\nImage stats: Minimum=%f, maximum=%f\n", mi, mx);
    biw = 1.000001 / 1000.0;
    // Histogram
    for (i = 0; i < sx * sx * 3; i++) {
        for (j = 0; j < 1000; j++)
            if (*(imT + i) >= (biw * j) && *(imT + i) < (biw * (j + 1))) {
                HDRhist[j]++;
                break;
            }
    }

    pct = .005 * (sx * sx * 3);
    mx = 0;
    for (j = 5; j < 990; j++) {
        mx += HDRhist[j];
        if (HDRhist[j + 5] - HDRhist[j - 5] > pct)
            break;
        if (mx > pct)
            break;
    }
    mi = (biw * (.90 * j));

    for (j = 990; j > 5; j--) {
        if (HDRhist[j - 5] - HDRhist[j + 5] > pct)
            break;
    }
    mx = (biw * (j + (.25 * (999 - j))));

    //fprintf(stderr, "Limit values chosen at min=%f, max=%f... normalizing image\n", mi, mx);

    for (i = 0; i < sx * sx * 3; i++) {
        *(imT + i) = *(imT + i) - mi;
        *(imT + i) = *(imT + i) / (mx - mi);
        if (*(imT + i) < 0.0)
            *(imT + i) = 0.0;
        if (*(imT + i) > 1.0)
            *(imT + i) = 1.0;
        *(imT + i) = pow(*(imT + i), .75);
    }

    bits24 = (unsigned char *)calloc(sx * sx * 3, sizeof(unsigned char));
    for (int i = 0; i < sx * sx * 3; i++)
        *(bits24 + i) = (unsigned char)(255.0 * (*(imT + i)));
    f = fopen(name, "wb+");
    if (f == NULL) {
        fprintf(stderr, "Unable to open file %s for output! No image written\n", name);
        return;
    }
    fprintf(f, "P6\n");
    fprintf(f, "# Output from PathTracer.c\n");
    fprintf(f, "%d %d\n", sx, sx);
    fprintf(f, "255\n");
    fwrite(bits24, sx * sx * 3 * sizeof(unsigned char), 1, f);
    fclose(f);
    

    free(bits24);
    free(imT);
}