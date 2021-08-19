#include "imageProcessor.h"

#include "utils.h"

#include <iostream>

// on linux: sudo apt-get install libpng-dev
#include <png.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct image *outImage;
char output_name[1024];  // Name of the output file for the raytraced .ppm image

// Just ignore warnings when loading the PNG
void warn_fn(png_structp, png_const_charp) {  }

struct image *readPNGimage(const char *filename){
    FILE *fp = fopen(filename, "rb");
  if (!fp) {
    return NULL;
  }
  png_structp png_ptr =
      png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, warn_fn);
  if (!png_ptr) {
    fclose(fp);
    return NULL;
  }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_read_struct(&png_ptr, NULL, NULL);
    fclose(fp);
    return NULL;
  }
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    fclose(fp);
    return NULL;
  }
  png_init_io(png_ptr, fp);
  png_read_info(png_ptr, info_ptr);
  png_uint_32 width, height;
  int bit_depth, color_type, interlace_type;
  png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
               &interlace_type, NULL, NULL);
  if (color_type == PNG_COLOR_TYPE_PALETTE) {
    png_set_palette_to_rgb(png_ptr);
  }
  if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
    png_set_expand_gray_1_2_4_to_8(png_ptr);
  }
  if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) {
    png_set_tRNS_to_alpha(png_ptr);
  }
  if (bit_depth == 16) {
    png_set_strip_16(png_ptr);
  }
  if (bit_depth < 8) {
    png_set_packing(png_ptr);
  }
  if (color_type == PNG_COLOR_TYPE_GRAY ||
      color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
    png_set_gray_to_rgb(png_ptr);
  }
  png_read_update_info(png_ptr, info_ptr);
  int rowbytes = png_get_rowbytes(png_ptr, info_ptr);
  int channels = png_get_channels(png_ptr, info_ptr);
  int size = rowbytes * height;

  png_byte *image = new png_byte[size];
  png_bytep row_pointers[height];
  for (int i = 0; i < height; i++) {
    row_pointers[i] = &image[i * rowbytes];
  }
  png_read_image(png_ptr, row_pointers);
  fclose(fp);
  png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
  struct image *img = new struct image;
  img->sx = width;
  img->sy = height;
  double* rgbdata = new double[size];
  for (int i = 0; i < size; i++) {
    rgbdata[i] = (double)image[i] / 255.0;
  }
  delete[] image;
  img->rgbdata = rgbdata;
  return img;
}

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

    im = (image *)calloc(1, sizeof(image));
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

/**
 * Save the image stored in `img` into the given PNG file
 */
bool PNGImageOutput(image *img, const char *filename) {
    printf("\nSaving Image\n");
    double* color_data = (double*)img->rgbdata;
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        return false;
    }
    png_structp png_ptr =
        png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
        fclose(fp);
        return false;
    }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(fp);
        return false;
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp);
        return false;
    }
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, img->sx, img->sy, 8, PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
                PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);
    png_byte *image = new png_byte[img->sx * img->sy * 3];
    for (int i = 0; i < img->sx * img->sy * 3; i++) {
        image[i] = (png_byte)(color_data[i] * 255);
    }
    png_bytep row_pointers[img->sy];
    for (int i = 0; i < img->sy; i++) {
        row_pointers[i] = &image[i * img->sx * 3];
    }
    png_write_image(png_ptr, row_pointers);
    png_write_end(png_ptr, info_ptr);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);

    delete[] image;
    return true;
}

void PPMImageOutput(image *im, const char *filename) {
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
