#ifndef TRACER_H
#define TRACER_H

#include "utility.h"
#include "ray.h"
#include "hittable_list.h"
#include "hit_record.h"
#include "material.h"

// #define DYNAMIC_BACKGROUND
bool DynamicLight;

inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
    float f = nf * fPdf, g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

color naive_tracer(const ray& r, const color& background, const hittable_list& world, int depth, std :: vector<shared_ptr<hittable>>& lights, bool include_emitted, unsigned int& seed){
    hit_record rec;
    scatter_record srec;

    if (depth <= 0) {
        // std :: cout << "Reach max depth" << std :: endl;
        return color(0, 0, 0);
    }

    // Debug
    // std :: cout << "the result of the hit test being" <<  world.hit(r, 0.000001f, Infinity, rec) << std :: endl;
    // std :: cout << "Max of world : " << world.WorldBound().max() << std :: endl;
    // std :: cout << "Min of world : " << world.WorldBound().min() << std :: endl;

    if (!world.hit(r, 0.f, Infinity, rec)) {

        // #ifdef DYNAMIC_BACKGROUND
        // vec3 unit_direction = unit_vector(r.direction());
        // auto t = 0.5 * (unit_direction.z() + 1.0);
        // return ((1.0 - t)*color(1.0, 1.0, 1.0) + t*color(0.3, 0.5, 0.8));
        // #endif

        // std :: cout << "hit nothing at depth : " << depth << std :: endl;
        if (DynamicLight) {
            vec3 unit_direction = unit_vector(r.direction());
            auto t = 0.5 * (unit_direction.z() + 1.0);
            return ((1.0 - t)*color(1.0, 1.0, 1.0) + t*color(0.3, 0.5, 0.8));
        }

        return background;
    }

    // std :: cout << "hit something at depth : " << depth << std :: endl;
    color emitted = rec.mat_ptr -> emitted();

    if (!rec.mat_ptr -> scatter(r, rec, srec, seed)) {
        return emitted;
    }
    if (srec.is_specular) {
        return srec.attenuation * naive_tracer(srec.specular_ray, background, world, depth-1, lights, include_emitted, seed);
    }

    ray scattered(rec.p, srec.pdf_ptr -> generate(seed));
    auto pdf_val = srec.pdf_ptr -> value(scattered.direction());
    if (pdf_val < 0.000001) pdf_val = 0.000001;

    float cosine_term = dot(rec.interpolated_normal, unit_vector(scattered.direction()));
    if (cosine_term < 0) {
        std :: cout << "Interpolated_normal wrong side" << std :: endl;
        return emitted;
    }
    color scattering_pdf = cosine_term * rec.mat_ptr -> scattering_bsdf(r, rec, scattered);

    return emitted + scattering_pdf * naive_tracer(scattered, background, world, depth - 1, lights, include_emitted, seed) / pdf_val;
}

color explicit_tracer(const ray& r, const color& background, const hittable_list& world, int depth, std :: vector<shared_ptr<hittable>>& lights, bool include_emitted, unsigned int& seed){
    hit_record rec;
    scatter_record srec;

    if (depth <= 0) return color(0, 0, 0);
    if (!world.hit(r, 0.f, Infinity, rec)) {

        // #ifdef DYNAMIC_BACKGROUND
        // vec3 unit_direction = unit_vector(r.direction());
        // auto t = 0.5 * (unit_direction.z() + 1.0);
        // return (1.0 - t)*color(1.0, 1.0, 1.0) + t*color(0.3, 0.5, 0.8);
        // #endif
        if (DynamicLight) {
            vec3 unit_direction = unit_vector(r.direction());
            auto t = 0.5 * (unit_direction.z() + 1.0);
            return ((1.0 - t)*color(1.0, 1.0, 1.0) + t*color(0.3, 0.5, 0.8));
        }

        // std :: cout << "hit nothing at depth : " << depth << std :: endl;
        return background;
    }
    
    // std :: cout << "hit somthing at depth : " << depth << std :: endl;
    color emitted;
    if  (include_emitted == true) {
        emitted = rec.mat_ptr -> emitted();
    }
    else {
        emitted = color(0, 0, 0);
    }

    // Can not scatter signal hitting a light source, cut the path here
    bool can_scatter = rec.mat_ptr -> scatter(r, rec, srec, seed);
    if (!can_scatter)
        return emitted;
    
    // In case of specular, just go on. Remember that since we haven't sample lights in specular case, if we hit light we have to add it manually
    if (srec.is_specular) {
        return srec.attenuation * explicit_tracer(srec.specular_ray, background, world, depth-1, lights, true, seed);
    }

    // Direct Lighting
    auto mat_ptr = srec.pdf_ptr;

    for (auto& light : lights) {

        // Sample light source with MIS
        float light_pdf = 0, scattering_pdf = 0;

        vec3 random_vec_to_light = light -> random_dir_to_light(rec.p, seed);
        vec3 f = rec.mat_ptr -> scattering_bsdf(r, rec, ray(rec.p, random_vec_to_light)) * dot(rec.interpolated_normal, unit_vector(random_vec_to_light));
        scattering_pdf = mat_ptr -> value(random_vec_to_light);

        // Ad Hoc
        if (dot(rec.interpolated_normal, unit_vector(random_vec_to_light)) > 0.0001) {
            float distance_just_before_light = random_vec_to_light.length() - 0.0001;
            hit_record temp_rec;
            ray shadow_ray(rec.p, random_vec_to_light);
            world.hit(shadow_ray, 0.f, distance_just_before_light, temp_rec);
            bool unoccluded = (temp_rec.t_distance * random_vec_to_light).length() > distance_just_before_light;

            if (unoccluded) {
                // Found nothing obstructing the light, calculating light pdf
                auto light_pdf = light -> pdf_value(rec.p, random_vec_to_light);

                // Debug
                if (light_pdf > 0.0001) {
                    float weigth = PowerHeuristic(1, light_pdf, 1, scattering_pdf);
                    emitted += f * (light -> emission()) * weigth / light_pdf;
                }
            }
        }

        // Sample BSDF
        ray scattered(rec.p, mat_ptr -> generate(seed));
        scattering_pdf = mat_ptr -> value(scattered.direction());
        f = rec.mat_ptr -> scattering_bsdf(r, rec, ray(rec.p, scattered.direction())) * dot(rec.interpolated_normal, unit_vector(scattered.direction()));

        light_pdf = light -> pdf_value(rec.p, scattered.direction());
        if (light_pdf == 0) {
            continue;
        }
        else {
            float weigth = PowerHeuristic(1, scattering_pdf, 1, light_pdf);
            
            hit_record temp_rec;
            world.hit(scattered,  0.f, Infinity, temp_rec);

            // Ad Hoc, but probably ok as long as the light placements are not too complicated. In such cases, just adjust the emission slightly between light objects
            if ((temp_rec.mat_ptr -> emitted()).x() == (light -> emission()).x()) {

                // Debug
                if (scattering_pdf > 0.0001) {
                    emitted += f * (light -> emission()) * weigth / scattering_pdf;
                }
            }
        }
        
    }

    // auto cosine_ptr = make_shared<cosine_pdf>(rec.interpolated_normal);
    ray scattered(rec.p, mat_ptr -> generate(seed));
    auto pdf_val = mat_ptr -> value(scattered.direction());
    if (pdf_val < 0.0001) pdf_val = 0.0001;

    float cosine_term = dot(rec.interpolated_normal, unit_vector(scattered.direction()));
    if (cosine_term < 0.0001) {
        return emitted;
    }
    color scattering_pdf = cosine_term * rec.mat_ptr -> scattering_bsdf(r, rec, scattered);

    return emitted + scattering_pdf * explicit_tracer(scattered, background, world, depth - 1, lights, false, seed) / pdf_val;
}

#endif