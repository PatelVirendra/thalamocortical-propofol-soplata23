#include "code_objects/IAMPA_PYso_RE_stateupdater_codeobject.h"
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



void _run_IAMPA_PYso_RE_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numN = 1;
double* const _array_IAMPA_PYso_RE_P_AMPA = _dynamic_array_IAMPA_PYso_RE_P_AMPA.empty()? 0 : &_dynamic_array_IAMPA_PYso_RE_P_AMPA[0];
const size_t _numP_AMPA = _dynamic_array_IAMPA_PYso_RE_P_AMPA.size();
const size_t _numdt = 1;
const double mV = 0.001;
double* const _array_IAMPA_PYso_RE_sAMPA_PYso_RE = _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.empty()? 0 : &_dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE[0];
const size_t _numsAMPA_PYso_RE = _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.size();
double* const _array_IAMPA_PYso_RE_tauAMPA = _dynamic_array_IAMPA_PYso_RE_tauAMPA.empty()? 0 : &_dynamic_array_IAMPA_PYso_RE_tauAMPA[0];
const size_t _numtauAMPA = _dynamic_array_IAMPA_PYso_RE_tauAMPA.size();
const size_t _numv_pre = 100;
int32_t* const _array_IAMPA_PYso_RE__synaptic_pre = _dynamic_array_IAMPA_PYso_RE__synaptic_pre.empty()? 0 : &_dynamic_array_IAMPA_PYso_RE__synaptic_pre[0];
const size_t _num_presynaptic_idx = _dynamic_array_IAMPA_PYso_RE__synaptic_pre.size();
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_IAMPA_PYso_RE_N = _array_IAMPA_PYso_RE_N;
    double* __restrict  _ptr_array_IAMPA_PYso_RE_P_AMPA = _array_IAMPA_PYso_RE_P_AMPA;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_IAMPA_PYso_RE_sAMPA_PYso_RE = _array_IAMPA_PYso_RE_sAMPA_PYso_RE;
    double* __restrict  _ptr_array_IAMPA_PYso_RE_tauAMPA = _array_IAMPA_PYso_RE_tauAMPA;
    double* __restrict  _ptr_array_PYso_group_v = _array_PYso_group_v;
    int32_t* __restrict  _ptr_array_IAMPA_PYso_RE__synaptic_pre = _array_IAMPA_PYso_RE__synaptic_pre;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = 1.0f*0.25/mV;


    const int _N = _array_IAMPA_PYso_RE_N[0];
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _presynaptic_idx = _ptr_array_IAMPA_PYso_RE__synaptic_pre[_idx];
        const double P_AMPA = _ptr_array_IAMPA_PYso_RE_P_AMPA[_idx];
        double sAMPA_PYso_RE = _ptr_array_IAMPA_PYso_RE_sAMPA_PYso_RE[_idx];
        const double tauAMPA = _ptr_array_IAMPA_PYso_RE_tauAMPA[_idx];
        const double v_pre = _ptr_array_PYso_group_v[_presynaptic_idx];
        const double _sAMPA_PYso_RE = (dt * (((P_AMPA * (1.0 - sAMPA_PYso_RE)) * (1.0 + tanh(_lio_1 * v_pre))) - (1.0f*sAMPA_PYso_RE/tauAMPA))) + sAMPA_PYso_RE;
        sAMPA_PYso_RE = _sAMPA_PYso_RE;
        _ptr_array_IAMPA_PYso_RE_sAMPA_PYso_RE[_idx] = sAMPA_PYso_RE;

    }

}


