#include "code_objects/IN_group_stateupdater_codeobject.h"
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



void _run_IN_group_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double Cm_IN = 0.009999999999999998;
const double EK_IN = - 0.09;
const double ELeak_IN = - 0.0638;
const double ENa_IN = 0.055;
const int64_t N = 20;
const size_t _numdt = 1;
const double gK_IN = 90.0;
const double gLeak_IN = 1.025;
const double gNa_IN = 350.0;
const size_t _numhNa_IN = 20;
const size_t _numiAMPA_PYso_IN = 20;
const size_t _numiAMPA_TC_IN = 20;
const size_t _numiGABAA_IN_IN = 20;
const size_t _numiNMDA_PYso_IN = 20;
const double mV = 0.001;
const double ms = 0.001;
const size_t _numnK_IN = 20;
const size_t _numnot_refractory = 20;
const size_t _numv = 20;
    ///// POINTERS ////////////
        
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_IN_group_hNa_IN = _array_IN_group_hNa_IN;
    double* __restrict  _ptr_array_IN_group_iAMPA_PYso_IN = _array_IN_group_iAMPA_PYso_IN;
    double* __restrict  _ptr_array_IN_group_iAMPA_TC_IN = _array_IN_group_iAMPA_TC_IN;
    double* __restrict  _ptr_array_IN_group_iGABAA_IN_IN = _array_IN_group_iGABAA_IN_IN;
    double* __restrict  _ptr_array_IN_group_iNMDA_PYso_IN = _array_IN_group_iNMDA_PYso_IN;
    double* __restrict  _ptr_array_IN_group_nK_IN = _array_IN_group_nK_IN;
    char* __restrict  _ptr_array_IN_group_not_refractory = _array_IN_group_not_refractory;
    double* __restrict  _ptr_array_IN_group_v = _array_IN_group_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = (- 25.0) * mV;
    const double _lio_2 = 1.0f*(- 5.0)/ms;
    const double _lio_3 = 1.0f*0.1/mV;
    const double _lio_4 = 1.0f*1.0/ms;
    const double _lio_5 = 1.0f*0.05/mV;
    const double _lio_6 = 1.0f*(- 0.625)/ms;
    const double _lio_7 = 1.0f*0.0125/mV;
    const double _lio_8 = 1.0f*dt/Cm_IN;
    const double _lio_9 = - gK_IN;
    const double _lio_10 = - EK_IN;
    const double _lio_11 = - ELeak_IN;
    const double _lio_12 = 1.0f*gNa_IN/(_brian_pow(ms, 3));
    const double _lio_13 = 1.0f*0.5/mV;
    const double _lio_14 = - ENa_IN;
    const double _lio_15 = 1.0f*20.0/ms;
    const double _lio_16 = 1.0f*0.05555555555555555/mV;


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double hNa_IN = _ptr_array_IN_group_hNa_IN[_idx];
        const double iAMPA_PYso_IN = _ptr_array_IN_group_iAMPA_PYso_IN[_idx];
        const double iAMPA_TC_IN = _ptr_array_IN_group_iAMPA_TC_IN[_idx];
        const double iGABAA_IN_IN = _ptr_array_IN_group_iGABAA_IN_IN[_idx];
        const double iNMDA_PYso_IN = _ptr_array_IN_group_iNMDA_PYso_IN[_idx];
        double nK_IN = _ptr_array_IN_group_nK_IN[_idx];
        char not_refractory = _ptr_array_IN_group_not_refractory[_idx];
        double v = _ptr_array_IN_group_v[_idx];
        if(!not_refractory)
            not_refractory = false || (! (v > _lio_1));
        else 
            not_refractory = true || (! (v > _lio_1));
        const double _hNa_IN = (dt * ((1.0f*(_lio_2 * hNa_IN)/(1.0 + exp((- 2.8000000000000003) - (_lio_3 * v)))) + (_lio_4 * ((0.35 - (0.35 * hNa_IN)) * exp((- 2.9000000000000004) - (_lio_5 * v)))))) + hNa_IN;
        const double _nK_IN = (dt * ((_lio_6 * (nK_IN * exp((- 0.55) - (_lio_7 * v)))) + (1.0f*(_lio_4 * ((1.0 - nK_IN) * (1.7 + (_lio_5 * v))))/(1.0 - exp((- 3.4000000000000004) - (_lio_3 * v)))))) + nK_IN;
        const double _v = v + (_lio_8 * ((((((_lio_9 * ((_brian_pow(nK_IN, 4)) * (_lio_10 + v))) + iAMPA_PYso_IN) + iAMPA_TC_IN) + iGABAA_IN_IN) + iNMDA_PYso_IN) - ((gLeak_IN * (_lio_11 + v)) + (1.0f*(_lio_12 * ((hNa_IN * (_brian_pow(17.5 + (_lio_13 * v), 3))) * (_lio_14 + v)))/((_brian_pow(1.0 - exp((- 3.5) - (_lio_3 * v)), 3)) * (_brian_pow((_lio_15 * exp((- 3.333333333333333) - (_lio_16 * v))) + (1.0f*(_lio_4 * (17.5 + (_lio_13 * v)))/(1.0 - exp((- 3.5) - (_lio_3 * v)))), 3)))))));
        hNa_IN = _hNa_IN;
        nK_IN = _nK_IN;
        v = _v;
        _ptr_array_IN_group_hNa_IN[_idx] = hNa_IN;
        _ptr_array_IN_group_nK_IN[_idx] = nK_IN;
        _ptr_array_IN_group_not_refractory[_idx] = not_refractory;
        _ptr_array_IN_group_v[_idx] = v;

    }

}


