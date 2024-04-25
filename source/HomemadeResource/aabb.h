// This file is for the axis-alligned bounding box
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : Addition of centroid and other functions useful mainly for acceleration structures

#ifndef AABB_H
#define AABB_H

#include "utility.h"
#include "ray.h"

class aabb {
    public :
        aabb() {}
        aabb(const point3& min, const point3& max) : minimum(min), maximum(max), centroid(0.5 * min + 0.5 * max) {}

        point3 min() const {return minimum;}
        point3 max() const {return maximum;}
        point3 cen() const {return centroid;}

        int MaximumExtent() const;
        float SurfaceArea() const;
        vec3 Diagonal() const;

        inline bool hit(const ray& r, float t_min, float t_max) const {
            for (int a = 0; a < 3; a++) {
                auto invD = 1.0 / r.direction()[a];
                auto t0 = (min()[a] - r.origin()[a]) * invD;
                auto t1 = (max()[a] - r.origin()[a]) * invD;
                if (invD < 0.f) {
                    std :: swap(t0, t1);
                }
                t_min = t0 > t_min ? t0 : t_min; // Use the bigger one
                t_max = t1 < t_max ? t1 : t_max; // Use the smaller one
                if (t_max <= t_min) return false;
            }
            return true;
        }
    
    private:
        point3 minimum;
        point3 maximum;
        point3 centroid;
};

// Surrounding box of either 2 smaller boxes, or a box and a point
aabb surrounding_box (aabb box0, aabb box1) {
    point3 small(fmin(box0.min().x(), box1.min().x()),
        fmin(box0.min().y(), box1.min().y()), fmin(box0.min().z(), box1.min().z()));
    point3 big(fmax(box0.max().x(), box1.max().x()),
        fmax(box0.max().y(), box1.max().y()), fmax(box0.max().z(), box1.max().z()));

    return aabb(small, big);
}
aabb surrounding_box (aabb box, point3 point) {
    point3 small(fmin(box.min().x(), point.x()),
        fmin(box.min().y(), point.y()), fmin(box.min().z(), point.z()));
    point3 big(fmax(box.max().x(), point.x()),
        fmax(box.max().y(), point.y()), fmax(box.max().z(), point.z()));

    return aabb(small, big);
}

vec3 aabb :: Diagonal() const {
    return maximum - minimum;
}
int aabb :: MaximumExtent() const {
    vec3 d = Diagonal();
    return MaxDimension(d);
} 
float aabb :: SurfaceArea() const {
    vec3 d = Diagonal();
    return 2.f * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
}

#endif