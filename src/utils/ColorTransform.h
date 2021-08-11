#pragma once

#include "imageProcessor.h"

template <typename colorSize> class ColorTransform{
  public:
    virtual colorSize* operator()(image& im) = 0;
};

template <typename colorSize>
class LinerToSRGB: public ColorTransform<colorSize>{
  public:
    LinerToSRGB(double e = 1.0):exposure(e){}
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