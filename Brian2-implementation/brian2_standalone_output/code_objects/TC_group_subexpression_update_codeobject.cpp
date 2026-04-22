#include "code_objects/TC_group_subexpression_update_codeobject.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<chrono>
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>

////// SUPPORT CODE ///////
namespace {
        
    template <typename T>
    static inline T _clip(const T value, const double a_min, const double a_max)
    {
        if (value < a_min)
            return a_min;
        if (value > a_max)
            return a_max;
        return value;
    }
    template < typename T1, typename T2 > struct _higher_type;
    template < > struct _higher_type<int32_t,int32_t> { typedef int32_t type; };
    template < > struct _higher_type<int32_t,int64_t> { typedef int64_t type; };
    template < > struct _higher_type<int32_t,float> { typedef float type; };
    template < > struct _higher_type<int32_t,double> { typedef double type; };
    template < > struct _higher_type<int32_t,long double> { typedef long double type; };
    template < > struct _higher_type<int64_t,int32_t> { typedef int64_t type; };
    template < > struct _higher_type<int64_t,int64_t> { typedef int64_t type; };
    template < > struct _higher_type<int64_t,float> { typedef float type; };
    template < > struct _higher_type<int64_t,double> { typedef double type; };
    template < > struct _higher_type<int64_t,long double> { typedef long double type; };
    template < > struct _higher_type<float,int32_t> { typedef float type; };
    template < > struct _higher_type<float,int64_t> { typedef float type; };
    template < > struct _higher_type<float,float> { typedef float type; };
    template < > struct _higher_type<float,double> { typedef double type; };
    template < > struct _higher_type<float,long double> { typedef long double type; };
    template < > struct _higher_type<double,int32_t> { typedef double type; };
    template < > struct _higher_type<double,int64_t> { typedef double type; };
    template < > struct _higher_type<double,float> { typedef double type; };
    template < > struct _higher_type<double,double> { typedef double type; };
    template < > struct _higher_type<double,long double> { typedef long double type; };
    template < > struct _higher_type<long double,int32_t> { typedef long double type; };
    template < > struct _higher_type<long double,int64_t> { typedef long double type; };
    template < > struct _higher_type<long double,float> { typedef long double type; };
    template < > struct _higher_type<long double,double> { typedef long double type; };
    template < > struct _higher_type<long double,long double> { typedef long double type; };
    // General template, used for floating point types
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_mod(T1 x, T2 y)
    {
        return x-y*floor(1.0*x/y);
    }
    // Specific implementations for integer types
    // (from Cython, see LICENSE file)
    template <>
    inline int32_t _brian_mod(int32_t x, int32_t y)
    {
        int32_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    template <>
    inline int64_t _brian_mod(int32_t x, int64_t y)
    {
        int64_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    template <>
    inline int64_t _brian_mod(int64_t x, int32_t y)
    {
        int64_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    template <>
    inline int64_t _brian_mod(int64_t x, int64_t y)
    {
        int64_t r = x % y;
        r += ((r != 0) & ((r ^ y) < 0)) * y;
        return r;
    }
    // General implementation, used for floating point types
    template < typename T1, typename T2 >
    static inline typename _higher_type<T1,T2>::type
    _brian_floordiv(T1 x, T2 y)
    {{
        return floor(1.0*x/y);
    }}
    // Specific implementations for integer types
    // (from Cython, see LICENSE file)
    template <>
    inline int32_t _brian_floordiv<int32_t, int32_t>(int32_t a, int32_t b) {
        int32_t q = a / b;
        int32_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    template <>
    inline int64_t _brian_floordiv<int32_t, int64_t>(int32_t a, int64_t b) {
        int64_t q = a / b;
        int64_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    template <>
    inline int64_t _brian_floordiv<int64_t, int>(int64_t a, int32_t b) {
        int64_t q = a / b;
        int64_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    template <>
    inline int64_t _brian_floordiv<int64_t, int64_t>(int64_t a, int64_t b) {
        int64_t q = a / b;
        int64_t r = a - q*b;
        q -= ((r != 0) & ((r ^ b) < 0));
        return q;
    }
    #ifdef _MSC_VER
    #define _brian_pow(x, y) (pow((double)(x), (y)))
    #else
    #define _brian_pow(x, y) (pow((x), (y)))
    #endif

}

////// HASH DEFINES ///////



void _run_TC_group_subexpression_update_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double A_TC = 5.1819378374737015;
const size_t _numCa_TC = 20;
const int64_t N = 20;
const size_t _numdCa_extra = 20;
const double gT_TC = 20.0;
const size_t _numhT_TC = 20;
const double inf = INFINITY;
const double mM = 1.0;
const double mV = 0.001;
const double second = 1.0;
const size_t _numv = 20;
const int64_t vShiftT_TC = 2;
    ///// POINTERS ////////////
        
    double* __restrict  _ptr_array_TC_group_Ca_TC = _array_TC_group_Ca_TC;
    double* __restrict  _ptr_array_TC_group_dCa_extra = _array_TC_group_dCa_extra;
    double* __restrict  _ptr_array_TC_group_hT_TC = _array_TC_group_hT_TC;
    double* __restrict  _ptr_array_TC_group_v = _array_TC_group_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double _lio_1 = 1.0f*(8314.41 * (309.15 * mV))/193692.0;
    const double _lio_2 = 2.0 * mM;
    const double _lio_3 = 57.0 + vShiftT_TC;
    const double _lio_4 = 1.0f*1.0/mV;
    const double _lio_5 = - A_TC;
    const double _lio_6 = 1.0f*(INFINITY * mM)/second;


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const double Ca_TC = _ptr_array_TC_group_Ca_TC[_idx];
        const double hT_TC = _ptr_array_TC_group_hT_TC[_idx];
        const double v = _ptr_array_TC_group_v[_idx];
        double dCa_extra;
        const double ET_TC = _lio_1 * log(1.0f*_lio_2/Ca_TC);
        const double Minf_T_TC = 1.0f*1.0/(1.0 + exp(0.16129032258064516 * (- (_lio_3 + (_lio_4 * v)))));
        const double iT_TC = gT_TC * (((_brian_pow(Minf_T_TC, 2)) * hT_TC) * (v - ET_TC));
        dCa_extra = _clip(_lio_5 * iT_TC, 0.0, _lio_6);
        _ptr_array_TC_group_dCa_extra[_idx] = dCa_extra;

    }

}


