#ifndef SPHERE_H
#define SPHERE_H

#include <iomanip>
#include "hittable.h"
#include "vec3.h"
#include "efloat.h"

class sphere : public hittable{
    public:
        sphere(point3 cen, float rad, shared_ptr<material> mat) : center(cen), radius(rad), mat_ptr(mat) {}
        bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const override;
        aabb WorldBound () const override {
            auto output_box = aabb(
            center - vec3(radius, radius, radius),
            center + vec3(radius, radius, radius));
            return output_box;
        }
        
    public:
        point3 center;
        float radius;
        shared_ptr<material> mat_ptr;
};

bool sphere :: hit(const ray& r, float t_min, float t_max, hit_record& rec) const {

    // std :: cout << "center x : " << center[0] << std :: endl;
    // std :: cout << "center y : " << center[1] << std :: endl;
    // std :: cout << "center z : " << center[2] << std :: endl;
    // std :: cout << "radiance : " << radius << std :: endl;
    // std :: cout << "ray origin : " << r.origin() << std :: endl;
    // std :: cout << "ray direction : " << r.direction() << std :: endl;

    EFloat transformed_origin[3];

    for (int i = 0; i < 3; i++) {
        transformed_origin[i] = EFloat(r.origin()[i]) - EFloat(center[i]);
    }
    
    EFloat ox = transformed_origin[0];
    EFloat oy = transformed_origin[1];
    EFloat oz = transformed_origin[2];
    EFloat dx(r.direction()[0]);
    EFloat dy(r.direction()[1]);
    EFloat dz(r.direction()[2]);

    // std :: cout << "ox : " << (float)ox << std :: endl;
    // std :: cout << "oy : " << (float)oy << std :: endl;
    // std :: cout << "dx : " << (float)dx << std :: endl;
    // std :: cout << "dy : " << (float)dy << std :: endl;

    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1))
        return false;

    if (t0.UpperBound() > t_max || t1.LowerBound() <= 0)
        return false;

    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;

        // std :: cout << "Lower bound lower than 0" << std :: endl;
        if (tShapeHit.UpperBound() > t_max)
            return false;
    }
    float t_hit = (float)tShapeHit;
    // std :: cout << "t_hit : " << t_hit << std :: endl;

    // Refine the t value to be of better seperation
    ray temp_r(point3(float(ox), float(oy), float(oz)), r.direction()); // temp_r, which is the transformed ray
    point3 p_hit = temp_r.at(t_hit);
    p_hit *=  std :: abs(radius) / p_hit.length();

    // std :: cout << "p_hit : " << p_hit << std :: endl;

    vec3 pError;
    for (int i = 0; i < 3; i++) {
        // One additional error bound from the final transformation back to world coordinate
        pError[i] = gamma(6) * std :: abs(p_hit[i]) + gamma(1) * std :: abs(center[i]);
        p_hit[i] += center[i];
    }

    rec.t_distance = t_hit;
    rec.p = p_hit;
    
    vec3 outward_normal = (rec.p - center) / radius;
    // std :: cout << "--------------------------outward : " << outward_normal << std :: endl;

    rec.set_geometric_normal(r, outward_normal);
    // rec.set_interpolated_normal(rec.geometric_normal);
    rec.set_interpolated_normal(outward_normal);

    rec.set_err(pError);
    rec.mat_ptr = mat_ptr;

    // std :: cout << "geo : " << rec.geometric_normal << std :: endl;
    // std :: cout << "interpo : " << rec.interpolated_normal << std :: endl;

    // std :: cout << "hit point : " << p_hit << std :: endl;
    // std :: cout << "center x : " << center[0] << std :: endl;
    // std :: cout << "center y : " << center[1] << std :: endl;
    // std :: cout << "center z : " << center[2] << std :: endl;
    // std :: cout << "radiance : " << radius << std :: endl;
    // std :: cout << "ray origin : " << r.origin() << std :: endl;
    // std :: cout << "ray direction : " << r.direction() << std :: endl;

    if (outward_normal[0] != outward_normal[0]) {
        std :: cout << "Got you" << std :: endl;
        return false;
    }

    return true;
}

#endif
