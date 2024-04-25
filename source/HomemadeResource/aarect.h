// Axis alligned rectangle
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : Just implementations of previously added data

#ifndef AARECT_H
#define AARECT_H

#include "utility.h"
#include "hittable.h"
#include "material.h"

class xy_rect : public hittable {
    public:
        xy_rect () {}
        xy_rect (float _x0, float _x1, float _y0, float _y1, float _k, shared_ptr<material> mat)
            : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {}

        virtual bool hit (const ray& r, float t_min, float t_max, hit_record& rec) const override;
        virtual aabb WorldBound() const override{
            return aabb(point3(x0, y0, k - 0.0001), point3(x1, y1, k + 0.0001));;
        }

        virtual float pdf_value (const point3& origin, const vec3& v) const override {
            hit_record temp_rec;
            if (!this -> hit(ray(origin, v), 0, Infinity, temp_rec)) return 0;

            auto area = (x1 - x0) * (y1 - y0);
            auto distance_squared = temp_rec.t_distance * temp_rec.t_distance * v.length_squared();
            auto cosine = fabs(dot(v, temp_rec.interpolated_normal)) / v.length();

            return distance_squared / (cosine * area);
            
        }
        virtual vec3 random_dir_to_light(const vec3& o, unsigned int& seed) const override {
            auto random_point = point3(random_float(x0, x1, seed), random_float(y0, y1, seed), k);
            vec3 random_dir_to_light = random_point - o;
            return random_dir_to_light;
        }

        virtual color emission() const {return mp -> emitted();}

    public:
        shared_ptr<material> mp;
        float x0, x1, y0, y1, k;
};

bool xy_rect :: hit (const ray& r, float t_min, float t_max, hit_record& rec) const {

    if (r.direction().z() == 0)
        return false;
    auto t = (k - r.origin().z()) / r.direction().z();
    if (t < t_min || t > t_max){
        return false;
    }
    
    auto x = r.origin().x() + t*r.direction().x();
    auto y = r.origin().y() + t*r.direction().y();

    if (x < x0 || x > x1 || y < y0 || y > y1) {
        return false;
    }

    auto outward_normal = vec3(0, 0, 1);
    rec.set_geometric_normal(r, outward_normal);
    rec.set_interpolated_normal(rec.geometric_normal);
    rec.mat_ptr = mp;

    // We don't need to be too precise with x and y, only z will be suffice to avoid self-intersection
    rec.set_err(vec3(0.f, 0.f, gamma(1)));
    rec.p = point3(x, y, k);

    rec.t_distance = t;

    return true;
}

class xz_rect : public hittable {
    public:
        xz_rect (float _x0, float _x1, float _z0, float _z1, double _k, shared_ptr<material> mat)
            : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {}

        xz_rect () {}

        virtual bool hit (const ray& r, float t_min, float t_max, hit_record& rec) const override;
        virtual aabb WorldBound() const override{
            return aabb(point3(x0, k - 0.0001, z0), point3(x1, k + 0.0001, z1));
        }
        virtual float pdf_value (const point3& origin, const vec3& v) const override {
            hit_record temp_rec;
            if (!this -> hit(ray(origin, v), 0, Infinity, temp_rec)) return 0;

            auto area = (x1 - x0) * (z1 - z0);
            auto distance_squared = temp_rec.t_distance * temp_rec.t_distance * v.length_squared();
            auto cosine = fabs(dot(v, temp_rec.interpolated_normal)) / v.length();

            // std :: cout << "distance squared : " << distance_squared << std :: endl;
            // std :: cout << "Interpolated normal : " << temp_rec.interpolated_normal << std :: endl;
            // std :: cout << "v : " << v << std :: endl;
            // std :: cout << "dot in cosine : " << fabs(dot(v, temp_rec.interpolated_normal)) << std :: endl;
            // std :: cout << "v length : " << v.length() << std :: endl;

            return distance_squared / (cosine * area);
            
        }
        virtual vec3 random_dir_to_light(const vec3& o, unsigned int& seed) const override {
            auto random_point = point3(random_float(x0, x1, seed), k, random_float(z0, z1, seed));
            vec3 random_dir_to_light = random_point - o;
            return random_dir_to_light;
        }

        virtual color emission() const {return mp -> emitted();}

    public:
        shared_ptr<material> mp;
        float x0, x1, z0, z1, k;
};

bool xz_rect :: hit (const ray& r, float t_min, float t_max, hit_record& rec) const {

    if (r.direction().y() == 0)
        return false;
    auto t = (k-r.origin().y()) / r.direction().y();
    if (t < t_min || t > t_max) {
        return false;
    }
        
    auto x = r.origin().x() + t*r.direction().x();
    auto z = r.origin().z() + t*r.direction().z();
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;

    auto outward_normal = vec3(0, 1, 0);
    rec.set_geometric_normal(r, outward_normal);
    rec.set_interpolated_normal(rec.geometric_normal);
    rec.mat_ptr = mp;

    rec.set_err(vec3(0.f, gamma(1), 0.f));
    rec.p = point3(x, k, z);

    rec.t_distance = t;
    // std :: cout << "Value of t : " << t << std :: endl;
    // std :: cout << "Value of t_max : " << t_max << std :: endl;

    return true;
}

class yz_rect : public hittable {
    public:
        yz_rect() {}

        yz_rect(float _y0, float _y1, float _z0, float _z1, float _k,
            shared_ptr<material> mat)
            : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const override;
        virtual aabb WorldBound() const override{
            return aabb(point3(k-0.0001, y0, z0), point3(k+0.0001, y1, z1));
        }
        virtual float pdf_value (const point3& origin, const vec3& v) const override {
            hit_record temp_rec;
            if (!this -> hit(ray(origin, v), 0, Infinity, temp_rec)) return 0;

            auto area = (y1 - y0) * (z1 - z0);
            auto distance_squared = temp_rec.t_distance * temp_rec.t_distance * v.length_squared();
            auto cosine = fabs(dot(v, temp_rec.interpolated_normal)) / v.length();

            return distance_squared / (cosine * area);
            
        }
        virtual vec3 random_dir_to_light(const vec3& o, unsigned int& seed) const override {
            auto random_point = point3(k, random_float(y0, y1, seed), random_float(z0, z1, seed));
            vec3 random_dir_to_light = random_point - o;
            return random_dir_to_light;
        }

        virtual color emission() const {return mp -> emitted();}

    public:
        shared_ptr<material> mp;
        float y0, y1, z0, z1, k;
};

bool yz_rect :: hit(const ray& r, float t_min, float t_max, hit_record& rec) const {

    if (r.direction().x() == 0)
        return false;
    auto t = (k-r.origin().x()) / r.direction().x();
    if (t < t_min || t > t_max)
        return false;

    auto y = r.origin().y() + t*r.direction().y();
    auto z = r.origin().z() + t*r.direction().z();
    if (y < y0 || y > y1 || z < z0 || z > z1)
        return false;

    auto outward_normal = vec3(1, 0, 0);
    rec.set_geometric_normal(r, outward_normal);
    rec.set_interpolated_normal(rec.geometric_normal);
    rec.mat_ptr = mp;
    
    rec.set_err(vec3(gamma(1), 0.f, 0.f));
    rec.p = point3(k, y, z);

    rec.t_distance = t;

    return true;
}

#endif
