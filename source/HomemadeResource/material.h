#ifndef MATERIAL_H
#define MATERIAL_H

#include "hit_record.h"
#include "utility.h"
#include "pdf.h"
#include "onb.h"

struct scatter_record {
    ray specular_ray;
    bool is_specular;
    color attenuation;
    shared_ptr<pdf> pdf_ptr;
};

class material{
    public:
        virtual bool scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed)
        const {return false;};
        virtual bool scatter () const = 0;

        virtual color scattering_bsdf (const ray& r_in, const hit_record& rec, const ray& scattered) const {return color();}

        virtual color emitted () const {
            // Impure virtual function with black as default return value
            return color(0, 0, 0);
        }
};

class lambertian : public material{
    public:
        lambertian(const color& a) : albedo(a) {}
        virtual bool scatter () const override {return true;}
        bool scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const override;
        // Pretty sure this thing is BSDF in disguise
        color scattering_bsdf (const ray& r_in, const hit_record& rec, const ray& scattered) const override {
            // auto cosine = dot(rec.interpolated_normal, unit_vector(scattered.direction()));
            // return cosine < 0 ? 0 : cosine / Pi;
            return albedo / Pi;
        }
    
    public:
        color albedo;
};
bool lambertian :: scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const {

    srec.is_specular = false;
    srec.attenuation = albedo;
    vec3 interpolated_normal = interpolated_allign_geometric(rec);
    rec.set_interpolated_normal(interpolated_normal);

    // std :: cout << "Inisde lambertian scatter. The value of interpolated normal feed to cosine_pdf is : " << interpolated_normal << std :: endl;

    srec.pdf_ptr = make_shared<cosine_pdf>(interpolated_normal);

    point3 spawn_orig = OffsetRayOrigin(rec.p, rec.err, rec.geometric_normal, rec.geometric_normal);
    rec.p = spawn_orig;

    return true;
}

class metal : public material{
    public:
        metal (const color& a, float f) : albedo(a), fuzz(f < 1 ? f : 1.f) {}
        virtual bool scatter () const override {return true;}
        bool scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const override;

    public:
        color albedo;
        float fuzz;
};
bool metal :: scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const{
    srec.is_specular = true;
    srec.attenuation = albedo;
    srec.pdf_ptr = 0;
    vec3 interpolated_normal = interpolated_allign_geometric(rec);
    vec3 reflected = reflect(unit_vector(r_in.direction()), rec.interpolated_normal);
    point3 spawn_orig = OffsetRayOrigin(rec.p, rec.err, rec.geometric_normal, rec.geometric_normal);
    srec.specular_ray = ray(spawn_orig, reflected+fuzz*random_in_unit_sphere(seed));
    rec.p = spawn_orig;

    return true;
}

class dielectric : public material{
    public:
        dielectric (double index_of_refraction) : ir(index_of_refraction) {}
        virtual bool scatter () const override {return true;}
        bool scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const override;

    private:
        static float reflectance (float cosine, float ref_index) {
            // Schlick's appoximation
            auto r0 = (1 - ref_index) / (1 + ref_index);
            r0 = r0 * r0;
            return r0 + (1 - r0) * pow((1 - cosine), 5);
        }

    public:
        float ir;
};
bool dielectric :: scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const {
    srec.is_specular = true;
    srec.pdf_ptr = 0;
    srec.attenuation = color(1.0, 1.0, 1.0);
    vec3 interpolated_normal = interpolated_allign_geometric(rec);
    float refraction_ratio = rec.front_face? (1.0 / ir) : ir;

    vec3 unit_direction = unit_vector(r_in.direction());
    float cos_theta = fmin(dot(-unit_direction, rec.geometric_normal), 1.0);    // Maybe interpolated_nomal will be more correct here
    float sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    bool cannot_refract = refraction_ratio * sin_theta > 1.0;
    vec3 direction;

    // For the new ray origin, which side it should be on depend on weather we choose to refract or reflect
    point3 spawn_orig;

    if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_float(seed)) {
        direction = reflect(unit_direction, interpolated_normal);

        // Should offset to the same side as normal
        spawn_orig = OffsetRayOrigin(rec.p, rec.err, rec.geometric_normal, rec.geometric_normal);
    }
    else {
        direction = refract(unit_direction, interpolated_normal, refraction_ratio);

        // Should offset to the opposite side of normal
        spawn_orig = OffsetRayOrigin(rec.p, rec.err, rec.geometric_normal, -rec.geometric_normal);
    }

    rec.p = spawn_orig;
    srec.specular_ray = ray(spawn_orig, direction);
    return true;
}

class beer_lambert_dielectric : public material {
    public:
        beer_lambert_dielectric (const float& index_of_refraction, const float& den, const color& volColor) :
        ir(index_of_refraction), inverse_ir(1.0 / index_of_refraction), density(den), volumeColor(volColor) {}

        virtual bool scatter () const override {return true;}
        bool scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const override;

    private:
        static float reflectance (float cosine, float ref_index) {
            // Schlick's appoximation
            auto r0 = (1 - ref_index) / (1 + ref_index);
            r0 = r0 * r0;
            return r0 + (1 - r0) * pow((1 - cosine), 5);
        }

    public:
        float density;
        color volumeColor;
        float ir; // Index of refraction
        float inverse_ir;
};
bool beer_lambert_dielectric :: scatter (const ray& r_in, hit_record& rec, scatter_record& srec, unsigned int & seed) const {

    float fuzz = 0.001;  // Manipultaing this value can give us new-looking material as well as more realistic glass
    
    float refraction_ratio;

    srec.is_specular = true;
    // srec.attenuation = volumeColor;
    srec.attenuation = color(1.0, 1.0, 1.0);

    // Potential bug. Must believe that rec.interpolated_normal indeed point to the correct "outside" of the object
    bool entering = dot(r_in.direction(), rec.interpolated_normal) < 0;
    vec3 ref_normal = entering ? rec.interpolated_normal : -rec.interpolated_normal;

    if (entering) {
        refraction_ratio = inverse_ir;
    }
    else {
        refraction_ratio = ir;

        // Exiting the object. Applying Beer-Lambert attenuation
        float t_distance = (rec.p - r_in.origin()).length();

        const color absorb = t_distance * density * volumeColor;
        // std :: cout << absorb << std :: endl;
        const color transparency = color(exp(-absorb[0]), exp(-absorb[1]), exp(-absorb[2]));
        srec.attenuation[0] *= transparency[0];
        srec.attenuation[1] *= transparency[1];
        srec.attenuation[2] *= transparency[2];
    }

    vec3 unit_direction = unit_vector(r_in.direction());
    float cos_theta = fmin(dot(-unit_direction, ref_normal), 1.0);
    float sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    bool cannot_refract = refraction_ratio * sin_theta > 1.0;
    vec3 direction;

    // For the new ray origin, which side it should be on depend on weather we choose to refract or reflect
    point3 spawn_orig;

    if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_float(seed)) {
        direction = reflect(unit_direction, ref_normal);

        // Should offset to the same side as normal
        spawn_orig = OffsetRayOrigin(rec.p, rec.err, rec.geometric_normal, rec.geometric_normal);
    }
    else {
        direction = refract(unit_direction, ref_normal, refraction_ratio);
        direction += random_in_unit_sphere(seed) * fuzz;

        // Should offset to the opposite side of normal
        spawn_orig = OffsetRayOrigin(rec.p, rec.err, rec.geometric_normal, -rec.geometric_normal);
    }

    // With this, we can reliably use t_distance for Beer-Lambert
    direction = unit_vector(direction);
    rec.p = spawn_orig;
    srec.specular_ray = ray(spawn_orig, direction);
    return true;
}

class diffuse_light : public material {
    public:
        diffuse_light (color a) : emit(a) {}
        virtual bool scatter () const override {return false;}

        // virtual bool scatter (const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered, unsigned int & seed) 
        //     const override {return false;}

        virtual color emitted () const override {
            return emit;
        }

    public:
        color emit;
};

#endif
