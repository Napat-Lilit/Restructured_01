#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "vec3.h"

struct  Matrix4x4
{
    Matrix4x4() : Matrix4x4 (1, 0, 0, 0,
    0, 1, 0, 0, 
    0, 0, 1, 0, 
    0, 0, 0, 1) {};
    Matrix4x4(float t00, float t01, float t02, float t03,
          float t10, float t11, float t12, float t13,
          float t20, float t21, float t22, float t23,
          float t30, float t31, float t32, float t33);
    float m[4][4];
};
Matrix4x4 :: Matrix4x4 (float t00, float t01, float t02, float t03,
float t10, float t11, float t12, float t13,
float t20, float t21, float t22, float t23,
float t30, float t31, float t32, float t33) {
    m[0][0] = t00, m[0][1] = t01, m[0][2] = t02, m[0][3] = t03,
    m[1][0] = t10, m[1][1] = t11, m[1][2] = t12, m[1][3] = t13,
    m[2][0] = t20, m[2][1] = t21, m[2][2] = t22, m[2][3] = t23,
    m[3][0] = t30, m[3][1] = t31, m[3][2] = t32, m[3][3] = t33;
}
Matrix4x4 Transpose(const Matrix4x4 &m) {
    return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0],
    m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1],
    m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2],
    m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]);
}
Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
    Matrix4x4 r;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            r.m[i][j] = m1.m[i][0] * m2.m[0][j] + 
                        m1.m[i][1] * m2.m[1][j] + 
                        m1.m[i][2] * m2.m[2][j] + 
                        m1.m[i][3] * m2.m[3][j];
    return r;
}

class Transform {
    public :
        Transform() {}
        Transform(const Matrix4x4 &m, const Matrix4x4 &mInv) : m(m), mInv(mInv) {}  // Force every transformation to provide inv when initialize
        Transform operator* (const Transform &t2) const;
        inline point3 point_transform(const point3& p) const;
        inline vec3 vec_transform(const vec3& p) const;
        inline vec3 norm_transform(const vec3& p) const;

    private :
        Matrix4x4 m, mInv;
};
Transform Translate(const vec3 &delta) {
    Matrix4x4 m(1, 0, 0, delta.x(),
                0, 1, 0, delta.y(),
                0, 0, 1, delta.z(), 
                0, 0, 0,       1);
    Matrix4x4 minv(1, 0, 0, -delta.x(),
                0, 1, 0, -delta.y(),
                0, 0, 1, -delta.z(), 
                0, 0, 0,        1);
    return Transform(m, minv);
}
Transform Scale(float x, float y, float z) {
    Matrix4x4 m(x, 0, 0, 0,
                0, y, 0, 0,
                0, 0, z, 0,
                0, 0, 0, 1);
    Matrix4x4 minv(1/x,   0,   0, 0,
                   0,   1/y,   0, 0,
                   0,     0, 1/z, 0,
                   0,     0,   0, 1);
    return Transform(m, minv);
}
Transform RotateX(float theta) {
    float sinTheta = std::sin(degrees_to_radians(theta));
    float cosTheta = std::cos(degrees_to_radians(theta));
    Matrix4x4 m(1,        0,         0, 0, 
                0, cosTheta, -sinTheta, 0,
                0, sinTheta,  cosTheta, 0,
                0,        0,         0, 1);
    return Transform(m, Transpose(m));
}
// For composite transformation
Transform Transform :: operator*(const Transform &t2) const {
    return Transform(Mul(m, t2.m),
                     Mul(t2.mInv, mInv));
}

inline point3 Transform :: point_transform(const point3& p) const {
    float x = p.x(), y = p.y(), z = p.z();
    float xp = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
    float yp = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
    float zp = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
    float wp = m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + m.m[3][3];
    if (wp == 1) return point3(xp, yp, zp);
    else         return point3(xp, yp, zp) / wp;
}
inline vec3 Transform :: vec_transform(const vec3& v) const {
    float x = v.x(), y = v.y(), z = v.z();
    return vec3(m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z,
                m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z,
                m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z);
}
inline vec3 Transform :: norm_transform(const vec3& n) const {
    float x = n.x(), y = n.y(), z = n.z();
    return vec3(mInv.m[0][0]*x + mInv.m[1][0]*y + mInv.m[2][0]*z,
                mInv.m[0][1]*x + mInv.m[1][1]*y + mInv.m[2][1]*z,
                mInv.m[0][2]*x + mInv.m[1][2]*y + mInv.m[2][2]*z);
}

#endif