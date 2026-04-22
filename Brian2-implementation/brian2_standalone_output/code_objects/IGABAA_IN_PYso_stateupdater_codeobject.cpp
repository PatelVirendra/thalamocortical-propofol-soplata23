#include "code_objects/IGABAA_IN_PYso_stateupdater_codeobject.h"
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



void _run_IGABAA_IN_PYso_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numN = 1;
double* const _array_IGABAA_IN_PYso_alpha_GABAA = _dynamic_array_IGABAA_IN_PYso_alpha_GABAA.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_alpha_GABAA[0];
const size_t _numalpha_GABAA = _dynamic_array_IGABAA_IN_PYso_alpha_GABAA.size();
double* const _array_IGABAA_IN_PYso_deprFactor = _dynamic_array_IGABAA_IN_PYso_deprFactor.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_deprFactor[0];
const size_t _numdeprFactor = _dynamic_array_IGABAA_IN_PYso_deprFactor.size();
const size_t _numdt = 1;
const double mV = 0.001;
double* const _array_IGABAA_IN_PYso_propoTauMult = _dynamic_array_IGABAA_IN_PYso_propoTauMult.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_propoTauMult[0];
const size_t _numpropoTauMult = _dynamic_array_IGABAA_IN_PYso_propoTauMult.size();
double* const _array_IGABAA_IN_PYso_res_GABAA = _dynamic_array_IGABAA_IN_PYso_res_GABAA.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_res_GABAA[0];
const size_t _numres_GABAA = _dynamic_array_IGABAA_IN_PYso_res_GABAA.size();
double* const _array_IGABAA_IN_PYso_sGABAA = _dynamic_array_IGABAA_IN_PYso_sGABAA.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_sGABAA[0];
const size_t _numsGABAA = _dynamic_array_IGABAA_IN_PYso_sGABAA.size();
double* const _array_IGABAA_IN_PYso_spike_activity = _dynamic_array_IGABAA_IN_PYso_spike_activity.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_spike_activity[0];
const size_t _numspike_activity = _dynamic_array_IGABAA_IN_PYso_spike_activity.size();
double* const _array_IGABAA_IN_PYso_tauGABAA = _dynamic_array_IGABAA_IN_PYso_tauGABAA.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_tauGABAA[0];
const size_t _numtauGABAA = _dynamic_array_IGABAA_IN_PYso_tauGABAA.size();
double* const _array_IGABAA_IN_PYso_tauRes_GABAA = _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso_tauRes_GABAA[0];
const size_t _numtauRes_GABAA = _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.size();
const size_t _numv_pre = 20;
int32_t* const _array_IGABAA_IN_PYso__synaptic_pre = _dynamic_array_IGABAA_IN_PYso__synaptic_pre.empty()? 0 : &_dynamic_array_IGABAA_IN_PYso__synaptic_pre[0];
const size_t _num_presynaptic_idx = _dynamic_array_IGABAA_IN_PYso__synaptic_pre.size();
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_IGABAA_IN_PYso_N = _array_IGABAA_IN_PYso_N;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_alpha_GABAA = _array_IGABAA_IN_PYso_alpha_GABAA;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_deprFactor = _array_IGABAA_IN_PYso_deprFactor;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_propoTauMult = _array_IGABAA_IN_PYso_propoTauMult;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_res_GABAA = _array_IGABAA_IN_PYso_res_GABAA;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_sGABAA = _array_IGABAA_IN_PYso_sGABAA;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_spike_activity = _array_IGABAA_IN_PYso_spike_activity;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_tauGABAA = _array_IGABAA_IN_PYso_tauGABAA;
    double* __restrict  _ptr_array_IGABAA_IN_PYso_tauRes_GABAA = _array_IGABAA_IN_PYso_tauRes_GABAA;
    double* __restrict  _ptr_array_IN_group_v = _array_IN_group_v;
    int32_t* __restrict  _ptr_array_IGABAA_IN_PYso__synaptic_pre = _array_IGABAA_IN_PYso__synaptic_pre;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = 1.0f*1.0/dt;
    const double _lio_2 = 1.0f*0.5/mV;
    const double _lio_3 = 20.0 * mV;


    const int _N = _array_IGABAA_IN_PYso_N[0];
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _presynaptic_idx = _ptr_array_IGABAA_IN_PYso__synaptic_pre[_idx];
        const double alpha_GABAA = _ptr_array_IGABAA_IN_PYso_alpha_GABAA[_idx];
        const double deprFactor = _ptr_array_IGABAA_IN_PYso_deprFactor[_idx];
        const double propoTauMult = _ptr_array_IGABAA_IN_PYso_propoTauMult[_idx];
        double res_GABAA = _ptr_array_IGABAA_IN_PYso_res_GABAA[_idx];
        double sGABAA = _ptr_array_IGABAA_IN_PYso_sGABAA[_idx];
        const double spike_activity = _ptr_array_IGABAA_IN_PYso_spike_activity[_idx];
        const double tauGABAA = _ptr_array_IGABAA_IN_PYso_tauGABAA[_idx];
        const double tauRes_GABAA = _ptr_array_IGABAA_IN_PYso_tauRes_GABAA[_idx];
        const double v_pre = _ptr_array_IN_group_v[_presynaptic_idx];
        const double _res_GABAA = (dt * ((spike_activity * ((1.0f*((- 1.0) + res_GABAA)/tauRes_GABAA) + (_lio_1 * ((deprFactor * res_GABAA) - res_GABAA)))) + (1.0f*(1.0 - res_GABAA)/tauRes_GABAA))) + res_GABAA;
        const double _sGABAA = (dt * ((1.0f*alpha_GABAA/(1.0 + exp(_lio_2 * (_lio_3 - v_pre)))) - (1.0f*sGABAA/(propoTauMult * tauGABAA)))) + sGABAA;
        res_GABAA = _res_GABAA;
        sGABAA = _sGABAA;
        _ptr_array_IGABAA_IN_PYso_res_GABAA[_idx] = res_GABAA;
        _ptr_array_IGABAA_IN_PYso_sGABAA[_idx] = sGABAA;

    }

}


