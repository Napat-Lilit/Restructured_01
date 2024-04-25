#ifndef PDF_H
#define PDF_H

#include "utility.h"
#include "onb.h"
#include "hittable.h"

class pdf {
    public:
        virtual ~pdf() {}
        virtual float value(const vec3& direction) const = 0;
        virtual vec3 generate(unsigned int& seed) const = 0;
};

// For cosine sampling on a hemisphere as relative to z axis
inline vec3 random_cosine_direction(unsigned int &seed) {
    auto r1 = random_float(seed);
    auto r2 = random_float(seed);
    auto z = std :: sqrt(1-r2);

    auto phi = 2 * Pi * r1;
    auto x = cos(phi) * std :: sqrt(r2);
    auto y = sin(phi) * std :: sqrt(r2);

    return vec3(x, y, z);
}
class cosine_pdf : public pdf {
    public:
        cosine_pdf(const vec3& w) { uvw.build_from_w(w); }

        virtual float value(const vec3& direction) const override {
            auto cosine = dot(unit_vector(direction), uvw.w());
            return (cosine <= 0) ? 0 : cosine / Pi;
        }

        virtual vec3 generate(unsigned int& seed) const override {
            return uvw.local(random_cosine_direction(seed));
        }

    public:
        onb uvw;
};

#endif