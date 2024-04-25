// Used to deal with error caused by the natural imprecision of floating point calculation
// Originally written in Physical Based Ray Tracing (PBRT)
// Adjustments from the original include : None

#ifndef EFLOAT_H
#define EFLOAT_H

#include "utility.h"

class EFloat {
    public:
        EFloat() {}
        EFloat(float v, float err = 0.f) : v(v) {
            if (err == 0.)
                low = high = v;
            else {
                low = NextFloatDown(v - err);
                high = NextFloatUp(v + err);
            }
        }

        EFloat operator+(EFloat ef) const {
            EFloat r;
            r.v = v + ef.v;
            r.low = NextFloatDown(LowerBound() + ef.LowerBound());
            r.high = NextFloatUp(UpperBound() + ef.UpperBound());
            return r;
        }

        explicit operator float() const { return v; }
        explicit operator double() const { return v; }
        float GetAbsoluteError() const { return NextFloatUp(std::max(std::abs(high - v), std::abs(v - low))); }
        float UpperBound() const { return high; }
        float LowerBound() const { return low; }

        EFloat operator-(EFloat ef) const {
            EFloat r;
            r.v = v - ef.v;
            r.low = NextFloatDown(LowerBound() - ef.UpperBound());
            r.high = NextFloatUp(UpperBound() - ef.LowerBound());
            return r;
        }

        EFloat operator*(EFloat ef) const {
            EFloat r;
            r.v = v * ef.v;
            float prod[4] = {
            LowerBound() * ef.LowerBound(), UpperBound() * ef.LowerBound(),
            LowerBound() * ef.UpperBound(), UpperBound() * ef.UpperBound()};
            r.low = NextFloatDown(
                std::min(std::min(prod[0], prod[1]), std::min(prod[2], prod[3])));
            r.high = NextFloatUp(
                std::max(std::max(prod[0], prod[1]), std::max(prod[2], prod[3])));
            return r;
        }

        EFloat operator/(EFloat ef) const {
            EFloat r;
            r.v = v / ef.v;
            if (ef.low < 0 && ef.high > 0) {
                // Bah. The interval we're dividing by straddles zero, so just
                // return an interval of everything.
                r.low = -Infinity;
                r.high = Infinity;
            }
            else {
                float div[4] = {
                    LowerBound() / ef.LowerBound(), UpperBound() / ef.LowerBound(),
                    LowerBound() / ef.UpperBound(), UpperBound() / ef.UpperBound()};
                r.low = NextFloatDown(
                    std::min(std::min(div[0], div[1]), std::min(div[2], div[3])));
                r.high = NextFloatUp(
                    std::max(std::max(div[0], div[1]), std::max(div[2], div[3])));
            }
            return r;
        }

        EFloat operator-() const {
            EFloat r;
            r.v = -v;
            r.low = -high;
            r.high = -low;
            return r;
        }

        inline bool operator==(EFloat fe) const { return v == fe.v; }

        EFloat &operator=(const EFloat &ef) {
            if (&ef != this) {
                v = ef.v;
                low = ef.low;
                high = ef.high;
            }
            return *this;
        }

    private:
        float v, low, high;

        // friend here to allow outside access to private data
        friend inline EFloat sqrt(EFloat fe);
        friend inline EFloat abs(EFloat fe);
        friend inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat *t0, EFloat *t1);
};

// EFloat Inline Functions
inline EFloat operator*(float f, EFloat fe) { return EFloat(f) * fe; }
inline EFloat operator/(float f, EFloat fe) { return EFloat(f) / fe; }
inline EFloat operator+(float f, EFloat fe) { return EFloat(f) + fe; }
inline EFloat operator-(float f, EFloat fe) { return EFloat(f) - fe; }

inline EFloat sqrt(EFloat fe) {
    EFloat r;
    r.v = std::sqrt(fe.v);
    r.low = NextFloatDown(std::sqrt(fe.low));
    r.high = NextFloatUp(std::sqrt(fe.high));
    return r;
}
inline EFloat abs(EFloat fe) {
    if (fe.low >= 0)
        // The entire interval is greater than zero, so we're all set.
        return fe;
    else if (fe.high <= 0) {
        // The entire interval is less than zero.
        EFloat r;
        r.v = -fe.v;
        r.low = -fe.high;
        r.high = -fe.low;
        return r;
    } else {
        // The interval straddles zero.
        EFloat r;
        r.v = std::abs(fe.v);
        r.low = 0;
        r.high = std::max(-fe.low, fe.high);
        return r;
    }
}

inline bool Quadratic(EFloat A, EFloat B, EFloat C, EFloat *t0, EFloat *t1) {
    // Find quadratic discriminant
    double discrim = (double)B.v * (double)B.v - 4. * (double)A.v * (double)C.v;
    if (discrim < 0.) return false;
    double rootDiscrim = std::sqrt(discrim);

    EFloat floatRootDiscrim(rootDiscrim, MachineEpsilon * rootDiscrim);
    EFloat q;

    if ((float)B < 0)
        q = -.5 * (B - floatRootDiscrim);
    else
        q = -.5 * (B + floatRootDiscrim);
    
    *t0 = q / A;
    *t1 = C / q;

    if ((float)*t0 > (float)*t1) std::swap(*t0, *t1);
    return true;
}

#endif