#ifndef BVH_H
#define BVH_H

#include "aabb.h"
#include "hittable.h"
#include "hittable_list.h"

#include <algorithm>

inline bool cen_compare(const shared_ptr<hittable> a, const shared_ptr<hittable> b, int axis) {
    // Of course, this works on the assumption that every hittable properly implemented the WorldBound() func.
    aabb box_a = a -> WorldBound();
    aabb box_b = b -> WorldBound();

    return box_a.cen()[axis] < box_b.cen()[axis];
}
bool cen_x_compare (const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return cen_compare(a, b, 0);
}
bool cen_y_compare (const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return cen_compare(a, b, 1);
}
bool cen_z_compare (const shared_ptr<hittable> a, const shared_ptr<hittable> b) {
    return cen_compare(a, b, 2);
}

size_t split_point_sah (std :: vector<shared_ptr<hittable>>& objects, 
size_t start, size_t end, int axis, aabb bounds, aabb cen_bounds);

class bvh_node : public hittable {
    public:
        // Delegating constructor
        bvh_node (hittable_list& list, unsigned int& seed)
            : bvh_node (list.objects, 0, list.objects.size(), seed) {}
        bvh_node (std :: vector<shared_ptr<hittable>>& src_objects, 
            size_t start, size_t end, unsigned int& seed);
        
        bool hit (const ray& r, float t_min, float t_max, hit_record& rec) const override;

        aabb WorldBound() const override {
            return box;
        }
        point3 Centroid() const {
            return box.cen();
        }
        int CutAxis() const {
            return cutAxis;
        }

    private:
        shared_ptr<hittable> left;
        shared_ptr<hittable> right;
        aabb box;
        int cutAxis;
};

bool bvh_node :: hit (const ray& r, float t_min, float t_max, hit_record& rec) const {
    // Don't hit the big box, nothing to look more
    if (box.hit(r, t_min, t_max) == false) {
        return false;
    }

    // std :: cout << "At least it hit the box" << std :: endl; 

    if (r.direction()[cutAxis] > 0.f) {
        // Come from left to right
        bool hit_left = left -> hit(r, t_min, t_max, rec);
        if (hit_left) {
            if (rec.p[cutAxis] < right -> WorldBound().min()[cutAxis]) return true;
            // std :: cout << "Interesting intersection" << std :: endl;
            bool hit_right = right -> hit(r, t_min, rec.t_distance, rec);
            return true;
        }

        bool hit_right = right -> hit(r, t_min, t_max, rec);
        return hit_right;
    }
    else if (r.direction()[cutAxis] < 0.f) {
        // Come from right to right
        bool hit_right = right -> hit(r, t_min, t_max, rec);
        if (hit_right) {
            if (rec.p[cutAxis] > left -> WorldBound().max()[cutAxis]) return true;
            // std :: cout << "Interesting intersection" << std :: endl;
            bool hit_left = left -> hit(r, t_min, rec.t_distance, rec);
            return true;
        }

        bool hit_left = left -> hit(r, t_min, t_max, rec);
        return hit_left;
    }
    else {
        bool hit_left = left -> hit(r, t_min, t_max, rec);
        bool hit_right = right -> hit(r, t_min, hit_left ? rec.t_distance : t_max, rec);
        return hit_left || hit_right;
    }
}

bvh_node :: bvh_node (std :: vector<shared_ptr<hittable>>& objects, 
    size_t start, size_t end, unsigned int& seed) {

    // BVH construction methods. 1 is random axis bisection. 2 is largest axis bisection. 3(default) is SAH.
    int method = 3;

    int axis = 0;
    auto comparator = cen_x_compare;    // Both will be overwrtitten later. These are just place holders

    // Dont forget that "start" begins with 0 but "end" begins with 1
    size_t object_span = end - start;
    size_t seperation_point = 0;
    aabb bounds;
    aabb centroid_bounds;

    switch (method)
    {
    case 1:
        axis = random_int(0, 2, seed);
        comparator = (axis == 0) ? cen_x_compare : (axis == 1) ? cen_y_compare : cen_z_compare;
        cutAxis = axis;

        if (object_span == 1) {
            left = right = objects[start];
        }
        else if (object_span == 2) {
            if (comparator(objects[start], objects[start + 1])){
                left = objects[start];
                right = objects[start + 1];
            }
            else {
                left = objects[start + 1];
                right = objects[start];
            }
        }
        else {
            std :: sort(objects.begin() + start, objects.begin() + end, comparator);
            seperation_point = start + object_span / 2;
            left = make_shared<bvh_node>(objects, start, seperation_point, seed);
            right = make_shared<bvh_node>(objects, seperation_point, end, seed);
        }
        break;

    case 2:
        centroid_bounds = aabb(objects[start] -> WorldBound().cen(), objects[start] -> WorldBound().cen());
        for (int i = start + 1; i < end; i++) {
            centroid_bounds = surrounding_box(centroid_bounds, objects[i] -> WorldBound().cen()); 
        }
        axis = centroid_bounds.MaximumExtent();
        comparator = (axis == 0) ? cen_x_compare : (axis == 1) ? cen_y_compare : cen_z_compare;
        cutAxis = axis;

        if (object_span == 1) {
            left = right = objects[start];
        }
        else if (object_span == 2) {
            if (comparator(objects[start], objects[start + 1])){
                left = objects[start];
                right = objects[start + 1];
            }
            else {
                left = objects[start + 1];
                right = objects[start];
            }
        }
        else {
            std :: sort(objects.begin() + start, objects.begin() + end, comparator);
            seperation_point = start + object_span / 2;
            left = make_shared<bvh_node>(objects, start, seperation_point, seed);
            right = make_shared<bvh_node>(objects, seperation_point, end, seed);
        }
        break;
    default:
        bounds = objects[start] -> WorldBound();
        centroid_bounds = aabb(objects[start] -> WorldBound().cen(), objects[start] -> WorldBound().cen());
        for (int i = start + 1; i < end; i++) {
            centroid_bounds = surrounding_box(centroid_bounds, objects[i] -> WorldBound().cen()); 
            bounds = surrounding_box(bounds, objects[i] -> WorldBound());
        }
        axis = centroid_bounds.MaximumExtent();
        comparator = (axis == 0) ? cen_x_compare : (axis == 1) ? cen_y_compare : cen_z_compare;
        cutAxis = axis;

        if (object_span <= 4) {
            if (object_span == 1) {
                left = right = objects[start];
            }
            else if (object_span == 2) {
                if (comparator(objects[start], objects[start + 1])){
                    left = objects[start];
                    right = objects[start + 1];
                }
                else {
                    left = objects[start + 1];
                    right = objects[start];
                }
            }
            else {
                std :: sort(objects.begin() + start, objects.begin() + end, comparator);

                seperation_point = start + object_span / 2;
                left = make_shared<bvh_node>(objects, start, seperation_point, seed);
                right = make_shared<bvh_node>(objects, seperation_point, end, seed);
            }
        }
        else {
            // Real SAH
            std :: sort(objects.begin() + start, objects.begin() + end, comparator);
            
            seperation_point = start + split_point_sah(objects, start, end, axis, bounds, centroid_bounds);
            left = make_shared<bvh_node>(objects, start, seperation_point, seed);
            right = make_shared<bvh_node>(objects, seperation_point, end, seed);
        }
        break;
    }

    aabb box_left, box_right;
    box_left = left -> WorldBound();
    box_right = right -> WorldBound();
    box = surrounding_box(box_left, box_right);
}

struct  BucketInfo
{
    int count = 0;
    aabb bounds;
};

size_t split_point_sah (std :: vector<shared_ptr<hittable>>& objects, 
size_t start, size_t end, int axis, aabb bounds, aabb cen_bounds) {

    constexpr int nBuckets = 12;
    BucketInfo buckets[nBuckets];

    // For each object in the range, determine which bucket it belongs to
    for (int i = start; i < end; i++) {
        int b = nBuckets * (objects[i] -> WorldBound().cen()[axis] - cen_bounds.min()[axis]) 
        / (cen_bounds.max()[axis] - cen_bounds.min()[axis]);
        if (b == nBuckets) b--;

        buckets[b].count++;
        if (buckets[b].count == 1) {
            // The first element
            buckets[b].bounds = objects[i] -> WorldBound();
        }
        else {
            buckets[b].bounds = surrounding_box(buckets[b].bounds, objects[i] -> WorldBound());
        }
    }

    // Compute cost of splitting after each bucket
    float cost[nBuckets - 1];
    for (int i = 0; i < nBuckets - 1; i++){
        aabb b0, b1;
        int count0 = 0;
        int count1 = 0;
        for (int j = 0; j <= i; j++) {
            if (count0 == 0 && buckets[j].count != 0) {
                b0 = buckets[j].bounds;
                count0 = buckets[j].count;
                continue;
            }
            if (buckets[j].count != 0) {
                b0 = surrounding_box(b0, buckets[j].bounds);
                count0 += buckets[j].count;
            }
        }
        for (int j = i + 1; j < nBuckets; j++) {
            if (count1 == 0 && buckets[j].count != 0) {
                b1 = buckets[j].bounds;
                count1 = buckets[j].count;
                continue;
            }
            if (buckets[j].count != 0) {
                b1 = surrounding_box(b1, buckets[j].bounds);
                count1 += buckets[j].count;
            }
        }
        cost[i] = .125f + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
    }

    //Find bucket to split that minimize SAH
    float minCost = cost[0];
    int minCostSplitBucket = 0;
    for (int i = 0; i < nBuckets - 1; i++) {
        if (cost[i] < minCost) {
            minCost = cost[i];
            minCostSplitBucket = i;
        }
    }
    //Find specific object to split that minimize SAH
    size_t splitPoint = 0;
    for (int i = 0; i <= minCostSplitBucket; i++) {
        splitPoint += buckets[i].count;
    }

    return splitPoint;
}

#endif