#pragma once

#include "imageProcessor.h"

template <typename colorSize> class ColorTransform{
  public:
    virtual colorSize* operator()(image& im) = 0;
};

template <typename colorSize>
class LinerToSRGB: public ColorTransform<colorSize>{
  public:
    LinerToSRGB(double e = 1.0):exposure(e){
        std::cout << "Using Color Transform LinerToSRGB with exposure " << exposure << std::endl;
    }
    colorSize* operator()(image& im){
        colorSize* transformedColors  = new colorSize[im.sx * im.sy * 3];
        colorSize* im_data = (colorSize* )im.rgbdata;
        std::copy(im_data, im_data + im.sx * im.sy * 3, transformedColors);

        for (int i = 0; i < im.sx * im.sy * 3; i++) {
            transformedColors[i] = linearToSRGB(transformedColors[i] * exposure);
        }

        return transformedColors;
    }
  private:
    colorSize linearToSRGB(colorSize L){
        colorSize S = L * 12.92;
        if (L > 0.0031308) {
        S = 1.055 * pow(L, 1.0 / 2.4) - 0.055;
        }
        return clamp01(S);
    }
    double exposure;
};

template <typename colorSize>
class LinerToPacosFunction: public ColorTransform<colorSize>{
  public:
    LinerToPacosFunction(){
        std::cout << "Using Color Transform LinerToPacosFunction" << std::endl;
    }

    colorSize* operator()(image& im){
        FILE *f;
        double *imT;
        double HDRhist[1000];
        int i, j;
        double mx, mi, biw, pct;
        int sx = im.sx;

        imT = (double *)calloc(sx * sx * 3, sizeof(double));
        memcpy(imT, im.rgbdata, sx * sx * 3 * sizeof(double));
  
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
        return imT;
    }
  private:
};