#ifndef HIT_RECORD_H
#define HIT_RECORD_H

#include "utility.h"
#include "aabb.h"

class material;

struct hit_record{
    point3 p;
    shared_ptr<material> mat_ptr;
    float t_distance;

    // We need both of them since geometric normal is more reliable in some cases
    bool front_face;
    vec3 geometric_normal;
    vec3 interpolated_normal;

    // Error associated with this intersection calculation
    vec3 err;

    inline void set_err (const vec3& _err) {
        err = _err;
    }
    // Forcing the geometric_normal to always go against the ray
    inline void set_geometric_normal (const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0.f;
        geometric_normal = front_face ? outward_normal : -outward_normal;
    }
    inline void set_interpolated_normal (const vec3& _interpolated_normal) {
        interpolated_normal = _interpolated_normal;
    }
    inline void set_intersection (const vec3& _new_p) {
        p = _new_p;
    }
};

// Ensure that interpolated normal doesn't point to the opposite side of geometric one
// Beware that in cases where interpolated normals really are messed up, this alone is not enough
inline vec3 interpolated_allign_geometric (const hit_record& rec) {
    bool correct_side = dot(rec.interpolated_normal, rec.geometric_normal) > 0.f;
    if (correct_side) {
        return rec.interpolated_normal;
    }
    else {
        return -rec.interpolated_normal;
    }
}

#endif