#ifndef SIREMATHS_MATHS_H
#define SIREMATHS_MATHS_H

#include "sireglobal.h"

#include <cmath>

#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <gsl/gsl_sys.h>

#include "constants.h"

#ifdef SIRE_MAX
#undef SIRE_MAX
#endif

#ifdef SIRE_MIN
#undef SIRE_MIN
#endif

#define SIRE_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define SIRE_MIN(x, y) (((x) < (y)) ? (x) : (y))

namespace SireMaths
{

typedef boost::variate_generator<boost::mt19937, boost::uniform_real<double> > UniformRand;
typedef boost::variate_generator<boost::mt19937, boost::uniform_int<qint32> > UniformRanInt;

/** Function to return a uniform double random number generator with a range
    from 'min' to 'max' */
inline UniformRand uniformRanGenerator(double min=0.0, double max=1.0)
{
    return UniformRand( boost::mt19937(), boost::uniform_real<double>(min,max) );
}

/** Function to return a uniform double random number generator with a range
    from 'min' to 'max' and seeded with 'seed' */
inline UniformRand uniformRanGenerator(double min, double max, uint32_t seed)
{
    return UniformRand( boost::mt19937(seed), boost::uniform_real<double>(min,max) );
}

/** Function to return a uniform int random number generate with a range from
    'min' to 'max' */
inline UniformRanInt uniformRanIntGenerator(qint32 min=0, qint32 max=100)
{
    return UniformRanInt( boost::mt19937(), boost::uniform_int<qint32>(min,max) );
}

/** Function to return a uniform int random number generate with a range from
    'min' to 'max' and seeded with 'seed' */
inline UniformRanInt uniformRanIntGenerator(qint32 min, qint32 max, uint32_t seed)
{
    return UniformRanInt( boost::mt19937(seed), boost::uniform_int<qint32>(min,max) );
}

/** Return x raised to the power 2 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_2(const T &x) { return x*x;   }

/** Return x raised to the power 3 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_3(const T &x) { return x*x*x; }

/** Return x raised to the power 4 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_4(const T &x) { T x2 = x*x;   return x2*x2;    }

/** Return x raised to the power 5 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_5(const T &x) { T x2 = x*x;   return x2*x2*x;  }

/** Return x raised to the power 6 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_6(const T &x) { T x2 = x*x;   return x2*x2*x2; }

/** Return x raised to the power 7 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_7(const T &x) { T x3 = x*x*x; return x3*x3*x;  }

/** Return x raised to the power 8 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_8(const T &x) { T x2 = x*x;   T x4 = x2*x2; return x4*x4; }

/** Return x raised to the power 9 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_9(const T &x) { T x3 = x*x*x; return x3*x3*x3; }

/** Return x raised to the power 10 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_10(const T &x) { T x2 = x*x; T x5 = x2*x2*x; return x5*x5; }

/** Return x raised to the power 11 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_11(const T &x) { T x2 = x*x; T x5 = x2*x2*x; return x5*x5*x; }

/** Return x raised to the power 12 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_12(const T &x) { T x2 = x*x; T x6 = x2*x2*x2; return x6*x6; }

/** Return x raised to the power 13 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_13(const T &x) { T x2 = x*x; T x6 = x2*x2*x2; return x6*x6*x; }

/** Return x raised to the power 14 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_14(const T &x) { T x2 = x*x; T x7 = x2*x2*x2*x; return x7*x7; }

/** Return x raised to the power 15 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_15(const T &x) { T x2 = x*x; T x5 = x2*x2*x; return x5*x5*x5; }

/** Return x raised to the power 16 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_16(const T &x) { T x2 = x*x; T x4 = x2*x2; return x4*x4*x4*x4; }

/** Return x raised to the power 17 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_17(const T &x) { T x2 = x*x; T x4 = x2*x2; return x4*x4*x4*x4*x; }

/** Return x raised to the power 18 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_18(const T &x) { T x3 = x*x*x; T x9 = x3*x3*x3; return x9*x9; }

/** Return x raised to the power 19 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_19(const T &x) { T x3 = x*x*x; T x9 = x3*x3*x3; return x9*x9*x; }

/** Return x raised to the power 20 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_20(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x8=x4*x4; return x8*x8*x4; }

/** Return x raised to the power 21 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_21(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x7=x4*x2*x; return x7*x7*x7; }

/** Return x raised to the power 22 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_22(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x6=x4*x2; return x6*x6*x4; }

/** Return x raised to the power 23 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_23(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x6=x4*x2; return x6*x6*x4*x; }

/** Return x raised to the power 24 */
template<typename T>
SIRE_INLINE_TEMPLATE
T pow_24(const T &x)
    { T x2 = x*x; T x4 = x2*x2; T x8=x4*x4; return x8*x8*x8; }

/** Internal function used by the 'pow' function - raise 'x' to
    the positive power 'n' */
inline double pow_pvt(double x, int n)
{
    switch(n)
    {
        case 1:
            return x;
        case 2:
            return pow_2(x);
        case 3:
            return pow_3(x);
        case 4:
            return pow_4(x);
        case 5:
            return pow_5(x);
        case 6:
            return pow_6(x);
        case 7:
            return pow_7(x);
        case 8:
            return pow_8(x);
        case 9:
            return pow_9(x);
        case 10:
            return pow_10(x);
        case 11:
            return pow_11(x);
        case 12:
            return pow_12(x);
        default:
            return std::pow(x, n);
    }
}

/** Return x raised to the power n */
inline double pow(double x, int n)
{
    if (n > 0)
        return pow_pvt(x,n);
    else if (n == 0)
        return double(1.0);
    else
        return double(1.0) / pow_pvt(x, -1*n);
}

/** Return whether or not an integer is even */
inline bool isEven(int val)
{
    return val % 2 == 0;
}

/** Return whether or not an integer is odd */
inline bool isOdd(int val)
{
    return val % 2 != 0;
}

/** Return whether or not this is an integer */
inline bool isInteger(double val)
{
    return (smallest_int <= val) and (val <= largest_int) and (val == double(int(val)));
}

/** Return pow(x,n) for non-integer n */
inline double pow(double x, double n)
{
    if (isInteger(n))
        return pow(x, int(n));
    else
        return std::pow(x,n);
}

/** Here is a non-fast version of fastInvSqrt */
inline double invSqrt(double x)
{
    return double(1.0) / std::sqrt(x);
}

/** Here is a non-fast version of fastInvSqrt */
inline float invSqrt(float x)
{
    return float(1.0) / std::sqrt(x);
}

/** Return true if two numbers are equal. This uses gsl_fcmp
    for the comparison, and the values must be equal to within
    a range of 2*delta, where delta = 2^k * epsilon, where
    k equals the maximum base-2 exponent of val0 and val1 as
    calculated via frexp(). This should normally be a value
    around 1e-6 */
inline bool areEqual(double val0, double val1, double epsilon)
{
    return not ( gsl_fcmp(val0,val1,epsilon) );
}

/** Return true if two numbers are equal */
inline bool areEqual(double val0, double val1)
{
    return not ( gsl_fcmp(val0,val1,1e-6) );
}

/** Return true if this number is equal to zero */
template<typename T>
SIRE_INLINE_TEMPLATE
bool isZero(const T &val)
{
    return (std::abs(val) < std::numeric_limits<T>::epsilon());
}

}

#endif
