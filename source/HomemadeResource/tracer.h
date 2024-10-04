#ifndef TRACER_H
#define TRACER_H

#include "utility.h"
#include "ray.h"
#include "hittable_list.h"
#include "hit_record.h"
#include "material.h"

#include "exrHandler.h"

extern bool hasOverlayBackgroundColor;
extern color BackgroundColorForOutput;

inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
    float f = nf * fPdf, g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

color naive_tracer(const ray& r, const color& BackgroundColorForIllumination, const hittable_list& world, int depth, bool firstShot, unsigned int& seed){
    hit_record rec;
    scatter_record srec;

    if (depth <= 0) {
        // std :: cout << "Reach max depth" << std :: endl;
        return color(0, 0, 0);
    }

    if (!world.hit(r, 0.f, Infinity, rec)) {
        if (firstShot && hasOverlayBackgroundColor) return BackgroundColorForOutput;
        else return BackgroundColorForIllumination;
    }

    // std :: cout << "hit something at depth : " << depth << std :: endl;
    color emitted = rec.mat_ptr -> emitted();

    if (!rec.mat_ptr -> scatter(r, rec, srec, seed)) {
        return emitted;
    }
    if (srec.is_specular) {
        return srec.attenuation * naive_tracer(srec.specular_ray, BackgroundColorForIllumination, world, depth-1, false, seed);
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

    return emitted + scattering_pdf * naive_tracer(scattered, BackgroundColorForIllumination, world, depth - 1, false, seed) / pdf_val;
}
color naive_tracer_exr(const ray& r, const exrHandler& exrObject ,const hittable_list& world, int depth, bool firstShot, unsigned int& seed){
    hit_record rec;
    scatter_record srec;

    if (depth <= 0) {
        // std :: cout << "Reach max depth" << std :: endl;
        return color(0, 0, 0);
    }

    if (!world.hit(r, 0.f, Infinity, rec)) {
        if (firstShot && hasOverlayBackgroundColor) return BackgroundColorForOutput;
        else return exrObject.getSampleFromRay(r);
    }

    // std :: cout << "hit something at depth : " << depth << std :: endl;
    color emitted = rec.mat_ptr -> emitted();

    if (!rec.mat_ptr -> scatter(r, rec, srec, seed)) {
        return emitted;
    }
    if (srec.is_specular) {
        return srec.attenuation * naive_tracer_exr(srec.specular_ray, exrObject, world, depth-1, false, seed);
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

    return emitted + scattering_pdf * naive_tracer_exr(scattered, exrObject, world, depth-1, false, seed) / pdf_val;
}

// These will be used for rendering normal and positional mapping, necessary for denoising purpose
color normal_color(const ray& r, const hittable_list& world){
    hit_record rec;

    if (!world.hit(r, 0.f, Infinity, rec)) {
        // Just return black if hit nothing
        return color();
    }

    // Prevent minus value
    color baseColor = rec.interpolated_normal;
    float red = (baseColor.x() + 1.0) / 2.0;
    float green = (baseColor.y() + 1.0) / 2.0;
    float blue = (baseColor.z() + 1.0) / 2.0;
    return color(red, green, blue);
}
color position_color(const ray& r, const hittable_list& world, const vec3& min_point_ref, 
    const vec3& max_point_ref, const vec3& non_hit_ref){
    hit_record rec;

    if (!world.hit(r, 0.f, Infinity, rec)) {
        return non_hit_ref;
    }

    // Prevent minus value
    float x_scale = max_point_ref.x() - min_point_ref.x();
    float y_scale = max_point_ref.y() - min_point_ref.y();
    float z_scale = max_point_ref.z() - min_point_ref.z();

    color baseColor = rec.p;
    float red = (baseColor.x() - min_point_ref.x()) / x_scale;
    float green = (baseColor.y() - min_point_ref.y()) / y_scale;
    float blue = (baseColor.z() - min_point_ref.z()) / z_scale;
    return color(red, green, blue);
}

#endif