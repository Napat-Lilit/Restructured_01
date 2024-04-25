// 
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : Minor changes to bounding box function

#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"
#include "ray.h"
#include <vector>

class hittable_list : public hittable{
    public:
        hittable_list () {}
        hittable_list (shared_ptr<hittable> object) { add(object); }
        void clear() { objects.clear(); }

        void add(shared_ptr<hittable> object) { objects.push_back(object); }

        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const override;
        virtual aabb WorldBound () const override{
            if (objects.empty()) {
                throw std :: invalid_argument("Give empthy vector to the hittable list");
            }

            aabb temp_box;
            bool first_box = true;

            for (const auto& object : objects) {
                if (first_box) {
                    first_box = false;
                    temp_box = object -> WorldBound();
                }
                else {
                    temp_box = surrounding_box(temp_box, (object -> WorldBound()));
                }
            }
            return temp_box;
        }
    
    public:
        std :: vector<shared_ptr<hittable>> objects;
};
bool hittable_list :: hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_any = false;
    auto closest_so_far = t_max;

    for (const auto& object : objects){
        if (object -> hit(r, t_min, closest_so_far, temp_rec)){
            hit_any = true;
            closest_so_far = temp_rec.t_distance;
            rec = temp_rec;
        }
    }
    return hit_any;
}

#endif
