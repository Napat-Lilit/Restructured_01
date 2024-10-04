// Vector3 class. This same class also aliasing as both point3 and color
// Originally written by Peter Shirley in Ray tracing in one weekend serie
// Adjustments from the original include : Apart from just one abs func, practically no change

#ifndef VEC3_H
#define VEC3_H

#include <iostream>
#include <cmath>

#include "utility.h"

class vec3 {
    public:
        vec3 () : e{0, 0, 0} {}
        vec3 (float e0, float e1, float e2) : e{e0, e1, e2} {}

        float x() const {return e[0];}
        float y() const {return e[1];}
        float z() const {return e[2];}

        void x(float newVal) {e[0] = newVal;}
        void y(float newVal) {e[1] = newVal;}
        void z(float newVal) {e[2] = newVal;}

        // Operator overload
        vec3 operator-() const {return vec3(-e[0], -e[1], -e[2]);}
        float operator[](int i) const {return e[i];}
        float& operator[](int i) {return e[i];}
        vec3& operator+=(const vec3 &v) {
            // We return this for other implementations such as chaining
            e[0] += v[0];
            e[1] += v[1];
            e[2] += v[2];
            return *this;
        }
        vec3& operator*=(const float t){
            e[0] *= t;
            e[1] *= t;
            e[2] *= t;
            return *this;
        }
        vec3& operator/=(const float t){
            return *this *= 1/t;
        }

        // Tools to deal with length-related questions
        float length_squared () const {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }
        float length () const{
            return std :: sqrt(length_squared()); // ---> std :: is important
        }
        bool near_zero () {
            const auto s = 1e-8;
            return ((fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s));
        }

        // Tools to deal with random vec3
        inline static vec3 random(unsigned int& seed) {
            return vec3(random_float(seed), random_float(seed), random_float(seed));
        }
        inline static vec3 random(float min, float max, unsigned int& seed) {
            return vec3(random_float(min, max, seed), random_float(min, max, seed), random_float(min, max, seed));
        }

    private:
        float e[3];
};

// Type aliasing
using point3 = vec3;
using color = vec3;

// Original utility functions
inline std::ostream& operator<<(std::ostream& out, const vec3& v){
    return out << v[0] << " " << v[1] << " " << v[2];
}
inline vec3 operator+(const vec3& a, const vec3& b){
    return vec3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
inline vec3 operator-(const vec3& a, const vec3& b){
    return vec3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
inline vec3 operator*(const vec3& a, const vec3& b){
    return vec3(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
inline vec3 operator*(const vec3& a, double t){
    return vec3(a[0]*t, a[1]*t, a[2]*t);
}
inline vec3 operator*(double t, const vec3& a){
    return a*t;
}
inline vec3 operator/(const vec3& a, double t){
    float invert_t = 1.f / t;
    return a*invert_t;
}
inline vec3 abs(const vec3& a) {
    return vec3(fabs(a.x()), fabs(a.y()), fabs(a.z()));
}
inline float dot(const vec3& a, const vec3& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline vec3 cross(const vec3& a, const vec3& b){
    return vec3(a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0]);
}
inline vec3 unit_vector(const vec3& a){
    return a / a.length();
}
inline vec3 random_in_unit_sphere(unsigned int & seed) {
    while (true)
    {
        auto p = vec3 :: random(-1.f, 1.f, seed);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}
inline vec3 random_unit_vector(unsigned int & seed) {
    return unit_vector(random_in_unit_sphere(seed));
}

// For defocus blur
inline vec3 random_in_unit_disk(unsigned int & seed){
    while (true)
    {
        auto p = vec3(random_float(-1, 1, seed), random_float(-1, 1, seed), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

// Tools to deal with reflect and refract
inline vec3 reflect(const vec3& v, const vec3& n){
    return v - 2*dot(v, n)*n;
}
inline vec3 refract (const vec3& u, const vec3& n, double etai_over_etat) {
    auto cos_theta = fmin(dot(-u, n), 1.f);
    vec3 r_out_perp = etai_over_etat * (u + cos_theta * n);
    vec3 r_out_parallel = - (std :: sqrt(fabs(1.f - r_out_perp.length_squared()))) * n;

    return r_out_perp + r_out_parallel;
}

// Origin spawning. Normal is the outward normal, while w is the side of the light that we want the new origin to be in
inline point3 OffsetRayOrigin (const point3& p, const vec3& p_error, const vec3& normal, const vec3& w) {
    float d = dot(p_error, abs(normal));
    // std :: cout << "In OffsetRayOrigin d : " << d << std :: endl;
    vec3 offset = d * normal;
    // std :: cout << "In OffsetRayOrigin offset : " << offset << std :: endl;
    if (dot(w, normal) < 0.f) {
        offset = -offset;
    }
    
    point3 po = p + offset;
    // std :: cout << "In OffsetRayOrigin po before rounding away : " << po << std :: endl;

    // Rounding away from p
    for (int i = 0; i < 3; i++) {
        if (offset[i] > 0)      po[i] = NextFloatUp(po[i]);
        else if (offset[i] < 0) po[i] = NextFloatDown(po[i]);
    }
    return po;
}
inline point3 Permute(const point3& p, int x, int y, int z) {
    return point3(p[x], p[y], p[z]);
}
inline int MaxDimension(const point3& p) {
    int i = 0;
    if (p[1] > p[0]) 
        i = 1;
    if (p[2] > p[i])
        i = 2;
    return i;
}
inline float MaxComponent(const point3& p) {
    float i = p[0];
    if (p[1] > p[0]) 
        i = p[1];
    if (p[2] > p[i])
        i = p[2];
    return i;
}

#endif