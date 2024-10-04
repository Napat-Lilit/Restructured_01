// Camera
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : Taken out the time

// Possible future adjustments : Stratified sampling

#ifndef CAMERA_H
#define CAMERA_H

#include "utility.h"
#include "vec3.h"

class camera{
    public:
        camera () {};
        camera (point3 lookfrom, point3 lookat, vec3 vup, float vfov, float aspect_ratio,
            float aperture, float focus_distance) {

            auto theta = degrees_to_radians(vfov);
            auto h = tan(theta / 2.0);
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;

            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);

            origin = lookfrom;
            horizontal = focus_distance * viewport_width * u;
            vertical = focus_distance * viewport_height * v;
            lower_left_corner = origin - horizontal/2.0 - vertical/2.0 - focus_distance * w;

            lens_radius = aperture / 2.0;
        }

        ray get_ray (float s, float t, unsigned int & seed){
            vec3 rd = lens_radius * random_in_unit_disk(seed);
            vec3 offset = u*rd.x() + v*rd.y();

            return ray(origin + offset, 
                lower_left_corner + s*horizontal + t*vertical - origin - offset);
        }

    public:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;

        vec3 u, v, w;
        float lens_radius;
};

#endif