#include "code_objects/IKNa_PYso_PYdr_stateupdater_codeobject.h"
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



void _run_IKNa_PYso_PYdr_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double ENaP_local = 0.055;
const double ENa_local = 0.055;
const size_t _numN = 1;
const double RPump = 18.0;
const double alphaNa = 10000000000.0;
const double areaDR = 3.5e-08;
const double areaSO = 1.5e-08;
double* const _array_IKNa_PYso_PYdr_concNa = _dynamic_array_IKNa_PYso_PYdr_concNa.empty()? 0 : &_dynamic_array_IKNa_PYso_PYdr_concNa[0];
const size_t _numconcNa = _dynamic_array_IKNa_PYso_PYdr_concNa.size();
const size_t _numdt = 1;
const double eqNaPumpTerm = 0.2025753861602528;
const double gNaP_local = 0.6859999999999999;
const double gNa_local = 500.0;
double* const _array_IKNa_PYso_PYdr_hNa_local = _dynamic_array_IKNa_PYso_PYdr_hNa_local.empty()? 0 : &_dynamic_array_IKNa_PYso_PYdr_hNa_local[0];
const size_t _numhNa_local = _dynamic_array_IKNa_PYso_PYdr_hNa_local.size();
const double mM = 1.0;
const double mV = 0.001;
const double ms = 0.001;
const int64_t phi_Na = 4;
const size_t _numv_post = 100;
const size_t _numv_pre = 100;
int32_t* const _array_IKNa_PYso_PYdr__synaptic_post = _dynamic_array_IKNa_PYso_PYdr__synaptic_post.empty()? 0 : &_dynamic_array_IKNa_PYso_PYdr__synaptic_post[0];
const size_t _num_postsynaptic_idx = _dynamic_array_IKNa_PYso_PYdr__synaptic_post.size();
int32_t* const _array_IKNa_PYso_PYdr__synaptic_pre = _dynamic_array_IKNa_PYso_PYdr__synaptic_pre.empty()? 0 : &_dynamic_array_IKNa_PYso_PYdr__synaptic_pre[0];
const size_t _num_presynaptic_idx = _dynamic_array_IKNa_PYso_PYdr__synaptic_pre.size();
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_IKNa_PYso_PYdr_N = _array_IKNa_PYso_PYdr_N;
    double* __restrict  _ptr_array_IKNa_PYso_PYdr_concNa = _array_IKNa_PYso_PYdr_concNa;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_IKNa_PYso_PYdr_hNa_local = _array_IKNa_PYso_PYdr_hNa_local;
    double* __restrict  _ptr_array_PYso_group_v = _array_PYso_group_v;
    double* __restrict  _ptr_array_PYdr_group_v = _array_PYdr_group_v;
    int32_t* __restrict  _ptr_array_IKNa_PYso_PYdr__synaptic_post = _array_IKNa_PYso_PYdr__synaptic_post;
    int32_t* __restrict  _ptr_array_IKNa_PYso_PYdr__synaptic_pre = _array_IKNa_PYso_PYdr__synaptic_pre;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = - RPump;
    const double _lio_2 = 0.0 - eqNaPumpTerm;
    const double _lio_3 = 3375.0 * (_brian_pow(mM, 3));
    const double _lio_4 = areaDR * gNaP_local;
    const double _lio_5 = - ENaP_local;
    const double _lio_6 = 1.0f*(- 0.12987012987013)/mV;
    const double _lio_7 = 1.0f*(areaSO * gNa_local)/(_brian_pow(ms, 3));
    const double _lio_8 = 1.0f*0.1/mV;
    const double _lio_9 = - ENa_local;
    const double _lio_10 = 1.0f*0.0455608884980535/ms;
    const double _lio_11 = 1.0f*0.08333333333333333/mV;
    const double _lio_12 = 1.0f*1.0/ms;
    const double _lio_13 = dt * phi_Na;


    const int _N = _array_IKNa_PYso_PYdr_N[0];
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _postsynaptic_idx = _ptr_array_IKNa_PYso_PYdr__synaptic_post[_idx];
        const int32_t _presynaptic_idx = _ptr_array_IKNa_PYso_PYdr__synaptic_pre[_idx];
        double concNa = _ptr_array_IKNa_PYso_PYdr_concNa[_idx];
        double hNa_local = _ptr_array_IKNa_PYso_PYdr_hNa_local[_idx];
        const double v_post = _ptr_array_PYso_group_v[_postsynaptic_idx];
        const double v_pre = _ptr_array_PYdr_group_v[_presynaptic_idx];
        const double _concNa = concNa + (dt * ((_lio_1 * (_lio_2 + (1.0f*(_brian_pow(concNa, 3))/(_lio_3 + (_brian_pow(concNa, 3)))))) - (alphaNa * ((1.0f*(_lio_4 * (_lio_5 + v_pre))/(_brian_pow(1.0 + (0.000721797280255039 * exp(_lio_6 * v_pre)), 3))) + (1.0f*(_lio_7 * ((hNa_local * (_brian_pow(3.3 + (_lio_8 * v_post), 3))) * (_lio_9 + v_post)))/((_brian_pow(1.0 - exp((- 3.3000000000000003) - (_lio_8 * v_post)), 3)) * (_brian_pow((_lio_10 * exp(_lio_11 * (- v_post))) + (1.0f*(_lio_12 * (3.3 + (_lio_8 * v_post)))/(1.0 - exp((- 3.3000000000000003) - (_lio_8 * v_post)))), 3))))))));
        const double _hNa_local = (_lio_13 * ((1.0f*(_lio_12 * (- hNa_local))/(1.0 + exp((- 2.0) - (_lio_8 * v_post)))) + (_lio_12 * ((0.07 - (0.07 * hNa_local)) * exp((- 5.0) - (_lio_8 * v_post)))))) + hNa_local;
        concNa = _concNa;
        hNa_local = _hNa_local;
        _ptr_array_IKNa_PYso_PYdr_concNa[_idx] = concNa;
        _ptr_array_IKNa_PYso_PYdr_hNa_local[_idx] = hNa_local;

    }

}


