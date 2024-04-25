#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <vector>
#include "vec3.h"
#include "hittable.h"
#include "transform.h"

extern bool TrustNormal;

struct TriangleMesh {
    TriangleMesh(const Transform& tran, int nTriangles, const int *vertexIndices, int nVertices, const vec3 *P, const vec3 *N);
    // TriangleMesh(const Transform& tran, int nTriangles, const std :: vector <int> *vertexIndices, int nVertices, const std :: vector <vec3>  *P, const std :: vector <vec3> *N);
    // TriangleMesh(const Transform& tran, int nTriangles, const int *vertexIndices, int nVertices, const vec3  *P, const vec3 *N);

    const int nTriangles, nVertices;
    std :: vector<int> vertexIndices;
    std::unique_ptr<vec3[]> p;
    std::unique_ptr<vec3[]> n;
};

TriangleMesh :: TriangleMesh(const Transform& tran, int nTriangles, const int *vertexIndices, int nVertices, const vec3 *P, const vec3 *N)
// TriangleMesh :: TriangleMesh(const Transform& tran, int nTriangles, const std :: vector <int> *vertexIndices, int nVertices, const std :: vector <vec3> *P, const std :: vector <vec3> *N)
// TriangleMesh :: TriangleMesh (const Transform& tran, int nTriangles, const int *vertexIndices, int nVertices, const vec3  *P, const vec3 *N)
: nTriangles(nTriangles), nVertices(nVertices), vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles) {
    // Transform mesh vertices
    p.reset(new vec3[nVertices]);
    for (int i = 0; i < nVertices; ++i) {
        p[i] = tran.point_transform(P[i]);
        // p[i] = tran.point_transform(P -> at(i));
    }
    // N if present
    if (N) {
        n.reset(new vec3[nVertices]);
        for (int i = 0; i < nVertices; ++i) {
            n[i] = tran.norm_transform(N[i]);
            // n[i] = tran.norm_transform(N -> at(i));
        }
    }
}

class Triangle : public hittable {
  public:
    Triangle(const std::shared_ptr<TriangleMesh> &mesh, int triNumber, shared_ptr<material> mat) : mesh(mesh), mat_ptr(mat) {
        v = &mesh->vertexIndices[3 * triNumber];
    }

    aabb WorldBound() const;
    bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;

  private:
    std::shared_ptr<TriangleMesh> mesh;
    const int *v; // Very weird using pointer here

    shared_ptr<material> mat_ptr;
};

aabb Triangle :: WorldBound () const {
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const point3 &p0 = mesh->p[v[0]]; // Works like an array
    const point3 &p1 = mesh->p[v[1]];
    const point3 &p2 = mesh->p[v[2]];

    float min_vx = fmin(p0.x(), p1.x());
    float min_vy = fmin(p0.y(), p1.y());
    float min_vz = fmin(p0.z(), p1.z());
    float max_vx = fmax(p0.x(), p1.x());
    float max_vy = fmax(p0.y(), p1.y());
    float max_vz = fmax(p0.z(), p1.z());

    point3 small(min_vx, min_vy, min_vz);
    point3 large(max_vx, max_vy, max_vz);
    
    return surrounding_box(aabb(small, large), p2);  
};
bool Triangle ::  hit (const ray& r, float t_min, float t_max, hit_record& rec) const {

    // std :: cout << "<<< Someone trying to hit a triangle" << std :: endl;
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const point3 &p0 = mesh->p[v[0]];
    const point3 &p1 = mesh->p[v[1]];
    const point3 &p2 = mesh->p[v[2]];
    
    // std :: cout << "p0t[2] : " << p0[2] << ", p1t[2] : " << p1[2] << ", p2t[2] : " << p2[2] << std :: endl;

    // Translate vertices based on ray origin
    point3 p0t = p0 - r.origin();
    point3 p1t = p1 - r.origin();
    point3 p2t = p2 - r.origin();

    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;

    int kz = MaxDimension(abs(r.direction()));
    int kx = kz + 1; if (kx == 3) kx = 0;
    int ky = kx + 1; if (ky == 3) ky = 0;
    
    vec3 new_d = Permute(r.direction(), kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);
    // std :: cout << "Permutation ok" << std :: endl;
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;

    // Apply shear transformation to translated vertex positions
    float sx = -new_d[0] / new_d[2];
    float sy = -new_d[1] / new_d[2];
    float sz = 1.f / new_d[2];
    p0t[0] += sx * p0t[2];
    p0t[1] += sy * p0t[2];
    p1t[0] += sx * p1t[2];
    p1t[1] += sy * p1t[2];
    p2t[0] += sx * p2t[2];
    p2t[1] += sy * p2t[2];
    // std :: cout << "Transformation ok" << std :: endl;
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;

    // Edge functions
    float e0 = p1t[0] * p2t[1] - p1t[1] * p2t[0];
    float e1 = p2t[0] * p0t[1] - p2t[1] * p0t[0];
    float e2 = p0t[0] * p1t[1] - p0t[1] * p1t[0];

    // Fall back to double precision test at triangle edges
    if (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f) {
        double p2txp1ty = (double)p2t[0] * (double)p1t[1];
        double p2typ1tx = (double)p2t[1] * (double)p1t[0];
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t[0] * (double)p2t[1];
        double p0typ2tx = (double)p0t[1] * (double)p2t[0];
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t[0] * (double)p0t[1];
        double p1typ0tx = (double)p1t[1] * (double)p0t[0];
        e2 = (float)(p1typ0tx - p1txp0ty);
    }
    // std :: cout << "Fall back to double ok" << std :: endl;

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    float det = e0 + e1 + e2;
    if (det == 0) return false; // The perfect side hit
    
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;
    // Pay the rest of the transformation price that we have withholded
    p0t[2] *= sz;
    p1t[2] *= sz;
    p2t[2] *= sz;
    float tScaled = e0 * p0t[2] + e1 * p1t[2] + e2 * p2t[2];
    if (det < 0 && (tScaled >= 0 || tScaled < t_max * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det))
        return false;

    // std :: cout << "Just before entering Bary ok" << std :: endl;
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;
    // Barycentric coordinates
    float invDet = 1 / det;
    float b0 = e0 * invDet; 
    float b1 = e1 * invDet; 
    float b2 = e2 * invDet; 
    float t = tScaled * invDet;

    // std :: cout << "Just before entering t check ok" << std :: endl;
    // Ensure that this t is indeed greater than 0
    // Delta_Z
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;
    float maxZ = MaxComponent(abs(point3(p0t[2], p1t[2], p2t[2])));
    // std :: cout << "maxZ : " << maxZ << std :: endl;
    float deltaZ = gamma(3) * maxZ;
    // std :: cout << "deltaZ : " << deltaZ << std :: endl;
    // Delta_X and Y
    float maxX = MaxComponent(abs(point3(p0t[0], p1t[0], p2t[0])));
    float maxY = MaxComponent(abs(point3(p0t[1], p1t[1], p2t[1])));
    float deltaX = gamma(5) * (maxX + maxZ);
    float deltaY = gamma(5) * (maxY + maxZ);
    // Delta_e
    // std :: cout << "Just before delta_e check ok" << std :: endl;
    float deltaE = 2 * (gamma(2) * maxX * maxY + deltaY * maxX + deltaX * maxY);

    float maxE = MaxComponent(abs(point3(e0, e1, e2)));
    float deltaT = 3 * (gamma(3) * maxE * maxZ + deltaE * maxZ + deltaZ * maxE) * std :: abs(invDet);
    if (t <= deltaT)
        return false;
    // std :: cout << "Just after entering t check ok" << std :: endl;

    // A hit is officially recognized, calculate associated error bounds
    float xAbsSum = (std :: abs(b0 * p0[0]) + std :: abs(b1 * p1[0]) + std :: abs(b2 * p2[0]));
    float yAbsSum = (std :: abs(b0 * p0[1]) + std :: abs(b1 * p1[1]) + std :: abs(b2 * p2[1]));
    float zAbsSum = (std :: abs(b0 * p0[2]) + std :: abs(b1 * p1[2]) + std :: abs(b2 * p2[2]));
    vec3 pError = gamma(7) * vec3(xAbsSum, yAbsSum, zAbsSum);

    // Using barycentric to find hit point
    point3 p_Hit = b0 * p0 + b1 * p1 + b2 * p2;

    // Fill in hit_record from triangle hit
    rec.t_distance = t;
    rec.p = p_Hit;

    vec3 v0v1 = p1 - p0;
    vec3 v0v2 = p2 - p0;
    vec3 geometric_normal = cross(v0v1, v0v2);
    geometric_normal = geometric_normal / geometric_normal.length();

    vec3 interpolated_normal;

    // std :: cout << "Before assign interpolated normal" << std :: endl;

    // if (mesh -> n) {
    if (TrustNormal) {
        interpolated_normal = unit_vector(b0 * mesh->n[v[0]] + b1 * mesh->n[v[1]] + b2 * mesh->n[v[2]]);
    }
    else {
        interpolated_normal = geometric_normal;
    }

    rec.set_geometric_normal(r, geometric_normal);
    rec.set_interpolated_normal(interpolated_normal);
    rec.set_err(pError);
    rec.mat_ptr = mat_ptr;    
    return true;
};

// Probably a very bad idea
class triangle : public hittable {
    public:
        triangle () {}
        triangle (point3 _v0, point3 _v1, point3 _v2, shared_ptr<material> mat) 
            : v0(_v0), v1(_v1), v2(_v2), mat_ptr(mat) {
                vec3 v0v1 = v1 - v0;
                vec3 v0v2 = v2 - v0;
                normal = cross(v0v1, v0v2);
                normal = normal / normal.length();
            }
        triangle (point3 _v0, point3 _v1, point3 _v2, vec3 _norm, shared_ptr<material> mat) 
            : v0(_v0), v1(_v1), v2(_v2), normal(_norm), mat_ptr(mat) {}

        bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
        aabb WorldBound() const;

    public:
        // const float epsilon = 0.0000001;
        point3 v0;
        point3 v1;
        point3 v2;
        vec3 normal;

        shared_ptr<material> mat_ptr;
};
aabb triangle :: WorldBound () const {
    // Get triangle vertices in _p0_, _p1_, and _p2_
    const point3 &p0 = v0; // Works like an array
    const point3 &p1 = v1;
    const point3 &p2 = v2;

    float min_vx = fmin(p0.x() - 0.00001f, p1.x() - 0.00001f);
    float min_vy = fmin(p0.y() - 0.00001f, p1.y() - 0.00001f);
    float min_vz = fmin(p0.z() - 0.00001f, p1.z() - 0.00001f);
    float max_vx = fmax(p0.x() + 0.00001f, p1.x() + 0.00001f);
    float max_vy = fmax(p0.y() + 0.00001f, p1.y() + 0.00001f);
    float max_vz = fmax(p0.z() + 0.00001f, p1.z() + 0.00001f);

    point3 small(min_vx, min_vy, min_vz);
    point3 large(max_vx, max_vy, max_vz);
    
    return surrounding_box(aabb(small, large), p2);  
};
bool triangle :: hit (const ray& r, float t_min, float t_max, hit_record& rec) const {

    point3 p0t = v0 - r.origin();
    point3 p1t = v1 - r.origin();
    point3 p2t = v2 - r.origin();

    int kz = MaxDimension(abs(r.direction()));
    int kx = kz + 1; if (kx == 3) kx = 0;
    int ky = kx + 1; if (ky == 3) ky = 0;
    
    vec3 new_d = Permute(r.direction(), kx, ky, kz);
    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);
    // std :: cout << "Permutation ok" << std :: endl;
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;

    // Apply shear transformation to translated vertex positions
    float sx = -new_d[0] / new_d[2];
    float sy = -new_d[1] / new_d[2];
    float sz = 1.f / new_d[2];
    p0t[0] += sx * p0t[2];
    p0t[1] += sy * p0t[2];
    p1t[0] += sx * p1t[2];
    p1t[1] += sy * p1t[2];
    p2t[0] += sx * p2t[2];
    p2t[1] += sy * p2t[2];
    // std :: cout << "Transformation ok" << std :: endl;
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;

    // Edge functions
    float e0 = p1t[0] * p2t[1] - p1t[1] * p2t[0];
    float e1 = p2t[0] * p0t[1] - p2t[1] * p0t[0];
    float e2 = p0t[0] * p1t[1] - p0t[1] * p1t[0];

    // Fall back to double precision test at triangle edges
    if (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f) {
        double p2txp1ty = (double)p2t[0] * (double)p1t[1];
        double p2typ1tx = (double)p2t[1] * (double)p1t[0];
        e0 = (float)(p2typ1tx - p2txp1ty);
        double p0txp2ty = (double)p0t[0] * (double)p2t[1];
        double p0typ2tx = (double)p0t[1] * (double)p2t[0];
        e1 = (float)(p0typ2tx - p0txp2ty);
        double p1txp0ty = (double)p1t[0] * (double)p0t[1];
        double p1typ0tx = (double)p1t[1] * (double)p0t[0];
        e2 = (float)(p1typ0tx - p1txp0ty);
    }
    // std :: cout << "Fall back to double ok" << std :: endl;

    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return false;
    float det = e0 + e1 + e2;
    if (det == 0) return false; // The perfect side hit
    
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;
    // Pay the rest of the transformation price that we have withholded
    p0t[2] *= sz;
    p1t[2] *= sz;
    p2t[2] *= sz;
    float tScaled = e0 * p0t[2] + e1 * p1t[2] + e2 * p2t[2];
    if (det < 0 && (tScaled >= 0 || tScaled < t_max * det))
        return false;
    else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det))
        return false;

    // std :: cout << "Just before entering Bary ok" << std :: endl;
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;
    // Barycentric coordinates
    float invDet = 1 / det;
    float b0 = e0 * invDet; 
    float b1 = e1 * invDet; 
    float b2 = e2 * invDet; 
    float t = tScaled * invDet;

    // std :: cout << "Just before entering t check ok" << std :: endl;
    // Ensure that this t is indeed greater than 0
    // Delta_Z
    // std :: cout << "p0t[2] : " << p0t[2] << ", p1t[2] : " << p1t[2] << ", p2t[2] : " << p2t[2] << std :: endl;
    float maxZ = MaxComponent(abs(point3(p0t[2], p1t[2], p2t[2])));
    // std :: cout << "maxZ : " << maxZ << std :: endl;
    float deltaZ = gamma(3) * maxZ;
    // std :: cout << "deltaZ : " << deltaZ << std :: endl;
    // Delta_X and Y
    float maxX = MaxComponent(abs(point3(p0t[0], p1t[0], p2t[0])));
    float maxY = MaxComponent(abs(point3(p0t[1], p1t[1], p2t[1])));
    float deltaX = gamma(5) * (maxX + maxZ);
    float deltaY = gamma(5) * (maxY + maxZ);
    // Delta_e
    // std :: cout << "Just before delta_e check ok" << std :: endl;
    float deltaE = 2 * (gamma(2) * maxX * maxY + deltaY * maxX + deltaX * maxY);

    float maxE = MaxComponent(abs(point3(e0, e1, e2)));
    float deltaT = 3 * (gamma(3) * maxE * maxZ + deltaE * maxZ + deltaZ * maxE) * std :: abs(invDet);
    if (t <= deltaT)
        return false;
    // std :: cout << "Just after entering t check ok" << std :: endl;

    // A hit is officially recognized, calculate associated error bounds
    float xAbsSum = (std :: abs(b0 * v0[0]) + std :: abs(b1 * v1[0]) + std :: abs(b2 * v2[0]));
    float yAbsSum = (std :: abs(b0 * v0[1]) + std :: abs(b1 * v1[1]) + std :: abs(b2 * v2[1]));
    float zAbsSum = (std :: abs(b0 * v0[2]) + std :: abs(b1 * v1[2]) + std :: abs(b2 * v2[2]));
    vec3 pError = gamma(7) * vec3(xAbsSum, yAbsSum, zAbsSum);

    // Using barycentric to find hit point
    point3 p_Hit = b0 * v0 + b1 * v1 + b2 * v2;

    // Fill in hit_record from triangle hit
    rec.t_distance = t;
    rec.p = p_Hit;

    vec3 v0v1 = v1 - v0;
    vec3 v0v2 = v2 - v0;
    vec3 geometric_normal = cross(v0v1, v0v2);
    geometric_normal = geometric_normal / geometric_normal.length();

    vec3 interpolated_normal;
    interpolated_normal = geometric_normal;

    // std :: cout << "Before assign interpolated normal" << std :: endl;

    // if (mesh -> n) {
    // if (TrustNormal) {
    //     interpolated_normal = unit_vector(b0 * mesh->n[v[0]] + b1 * mesh->n[v[1]] + b2 * mesh->n[v[2]]);
    // }
    // else {
    //     interpolated_normal = geometric_normal;
    // }

    rec.set_geometric_normal(r, geometric_normal);
    rec.set_interpolated_normal(interpolated_normal);
    rec.set_err(pError);
    rec.mat_ptr = mat_ptr;    
    return true;
}

#endif