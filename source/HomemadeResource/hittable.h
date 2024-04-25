// Used to deal with hittable class. Also hold an implementation of rotation. 
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : Minor adjustments in the bounding_box virtual func

// Possible bugs : Rotate normals back into world space surely suffered additional floating point error. How could we address them ?

#ifndef HITTABLE_H
#define HITTABLE_H

#include "hit_record.h"
#include "ray.h"

class hittable{
    public:
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
        virtual aabb WorldBound () const = 0;   // "Pure" virtual, forcing every hittable to implement this func
        virtual float pdf_value(const point3& origin, const vec3& v) const {
            return 0;
        }
        virtual vec3 random_dir_to_light(const vec3& o, unsigned int& seed) const {return vec3(1, 0, 0);}    // Literally have 0 idea why the default should be this
        virtual color emission() const {return color(0, 0, 0);}
};

class rotate_y : public hittable {
    public:
        rotate_y(shared_ptr<hittable> p, float angle);

        virtual bool hit (const ray& r, float t_min, float t_max, hit_record& rec) const override;
        virtual aabb WorldBound () const override {
            return bbox;
        }

    public:
        shared_ptr<hittable> ptr;
        float sin_theta;
        float cos_theta;
        aabb bbox;
};
rotate_y :: rotate_y (shared_ptr<hittable> p, float angle) : ptr(p) {
    auto radians = degrees_to_radians(angle);
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    bbox = ptr -> WorldBound();

    point3 min(Infinity, Infinity, Infinity);
    point3 max(-Infinity, -Infinity, -Infinity);

    // Check for min and max of each axis
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 2; k++){
                auto x = i*bbox.max().x() + (1-i)*bbox.min().x();
                auto y = j*bbox.max().y() + (1-j)*bbox.min().y();
                auto z = k*bbox.max().z() + (1-k)*bbox.min().z();

                auto newx =  cos_theta*x + sin_theta*z;
                auto newz = -sin_theta*x + cos_theta*z;

                vec3 tester(newx, y, newz);

                for (int c = 0; c < 3; c++) {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
    }
    bbox = aabb(min, max);
}
bool rotate_y :: hit (const ray& r, float t_min, float t_max, hit_record& rec) const {
    auto origin = r.origin();
    auto direction = r.direction();

    // The idea is to rotate the ray instead of the object
    origin[0] = cos_theta*r.origin()[0] - sin_theta*r.origin()[2];
    origin[2] = sin_theta*r.origin()[0] + cos_theta*r.origin()[2];

    direction[0] = cos_theta*r.direction()[0] - sin_theta*r.direction()[2];
    direction[2] = sin_theta*r.direction()[0] + cos_theta*r.direction()[2];

    ray rotated_r(origin, direction);

    if (!ptr->hit(rotated_r, t_min, t_max, rec))
        return false;

    auto p = rec.p;
    auto err = rec.err;
    auto geometric_normal = rec.geometric_normal;
    auto interpolated_normal = rec.interpolated_normal;

    p[0] =  cos_theta*rec.p[0] + sin_theta*rec.p[2];
    p[2] = -sin_theta*rec.p[0] + cos_theta*rec.p[2];
    // Accounting for additional error in p, caused by the above transformation
    float ab_cos = abs(cos_theta);
    float ab_sin = abs(sin_theta);
    err[0] = ab_cos * rec.err[0] + ab_sin * rec.err[2] + gamma(2) * (ab_cos * (abs(rec.p[0]) + rec.err[0]) + ab_sin * (abs(rec.p[2]) + rec.err[2]));
    err[2] = ab_sin * rec.err[0] + ab_cos * rec.err[2] + gamma(2) * (ab_sin * (abs(rec.p[0]) + rec.err[0]) + ab_cos * (abs(rec.p[2]) + rec.err[2]));
    geometric_normal[0] =  cos_theta*rec.geometric_normal[0] + sin_theta*rec.geometric_normal[2];
    geometric_normal[2] = -sin_theta*rec.geometric_normal[0] + cos_theta*rec.geometric_normal[2];
    interpolated_normal[0] = cos_theta*rec.interpolated_normal[0] + sin_theta*rec.interpolated_normal[2];
    interpolated_normal[2] = -sin_theta*rec.interpolated_normal[0] + cos_theta*rec.interpolated_normal[2];

    rec.p = p;
    rec.set_err(err);
    rec.set_geometric_normal(r, geometric_normal);
    rec.set_interpolated_normal(interpolated_normal);

    return true;
}

#endif