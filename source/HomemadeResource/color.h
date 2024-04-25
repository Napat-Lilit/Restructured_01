// To output to a file, given a color object
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : None

#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"
#include <iostream>

void write_color(std::ostream& out, const color& pixel_color, int samples_per_pixel){

    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    auto scale = 1.0 / samples_per_pixel;

    if (b < 0 || b > 100000) {
        std :: cout << "pixel color value is weird b : " << b << std ::endl;
    }

    // Gamma-correct for gamma = 2.0
    r = std :: sqrt(scale * r);
    g = std :: sqrt(scale * g);
    b = std :: sqrt(scale * b);

    // Clamp to avoid the color data going over 256 limit
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << " "
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << " "
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << "\n";
}

#endif
