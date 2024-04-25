// This file is for the basic utilities underlining the whole system
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : Functions used to deal with float error analysis

#ifndef UTILITY_H
#define UTILITY_H

#include <cmath>
#include <limits>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <memory>

// Using
using std :: shared_ptr;
using std :: make_shared;

// Constants
static constexpr double Pi = 3.1415926535897932385;
static constexpr double InverseRandomMax = 1.0 / (RAND_MAX + 1.0);
static constexpr float Infinity = std :: numeric_limits<float> :: infinity();
static constexpr float MachineEpsilon= std :: numeric_limits<float> :: epsilon() * 0.5;

// Original utility functions
inline float degrees_to_radians (float degrees) {
    // If called frequently enough, converting the division into a constexpr maybe a good idea
    return degrees * Pi / 180.f;
}
inline float random_float (unsigned int &seed) {
    return rand_r(&seed) * InverseRandomMax;
}
inline float random_float (float min, float max, unsigned int &seed) {
    return min + (max - min) * random_float(seed);
}
inline int random_int (int min, int max, unsigned int &seed) {
    return static_cast<int>(random_float(min, max + 1, seed));
}
inline float clamp (float x, float min, float max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

// Utility functions from PBRT
inline uint32_t FloatToBits (float f) {
    uint32_t ui;
    memcpy(&ui, &f, sizeof(uint32_t));
    return ui;
}
inline float BitsToFloat (uint32_t ui){
    float f;
    memcpy(&f, &ui, sizeof(uint32_t));
    return f;
}
inline float NextFloatUp (float f) {
    if (std :: isinf(f) && f > 0.f)
        return f;
    if (f == -0.f)
        f = 0.f;
    uint32_t ui = FloatToBits(f);

    // By definition, plus and minus valued float is defined a little bit differently
    if (f >= 0) ++ui;
    else --ui;
    return BitsToFloat(ui);
}
inline float NextFloatDown (float f) {
    if (std :: isinf(f) && f< 0.f)
        return f;
    if (f == 0.f)
        f = -0.f;
    
    uint32_t ui = FloatToBits(f);
    // By definition, plus and minus valued float is defined a little bit differently
    if (f > 0) --ui;
    else ++ui;
    return BitsToFloat(ui);
}
inline constexpr float gamma (int n) {
    return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

#include "vec3.h"

#endif