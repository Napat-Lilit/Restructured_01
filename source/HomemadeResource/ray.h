// Ray object.
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : We don't need the time data, so it is taken out

#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray{
    public :
        ray() {}
        ray(const point3& origin, const vec3& direction) : orig(origin), dir(direction) {}

        point3 origin() const {return orig;}
        vec3 direction() const {return dir;} 

        vec3 at(float t) const {
            return orig + t * dir;
        }

    private:
        point3 orig;
        vec3 dir;
};

#endif