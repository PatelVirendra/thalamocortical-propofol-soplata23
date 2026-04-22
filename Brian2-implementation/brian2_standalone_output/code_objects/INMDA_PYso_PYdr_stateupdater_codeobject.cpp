#include "code_objects/INMDA_PYso_PYdr_stateupdater_codeobject.h"
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



void _run_INMDA_PYso_PYdr_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numN = 1;
double* const _array_INMDA_PYso_PYdr_alphaS_NMDA = _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA[0];
const size_t _numalphaS_NMDA = _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.size();
double* const _array_INMDA_PYso_PYdr_alphaX_NMDA = _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA[0];
const size_t _numalphaX_NMDA = _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.size();
double* const _array_INMDA_PYso_PYdr_deprFactor = _dynamic_array_INMDA_PYso_PYdr_deprFactor.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_deprFactor[0];
const size_t _numdeprFactor = _dynamic_array_INMDA_PYso_PYdr_deprFactor.size();
const size_t _numdt = 1;
const double mV = 0.001;
double* const _array_INMDA_PYso_PYdr_res_NMDA = _dynamic_array_INMDA_PYso_PYdr_res_NMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_res_NMDA[0];
const size_t _numres_NMDA = _dynamic_array_INMDA_PYso_PYdr_res_NMDA.size();
double* const _array_INMDA_PYso_PYdr_sNMDA = _dynamic_array_INMDA_PYso_PYdr_sNMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_sNMDA[0];
const size_t _numsNMDA = _dynamic_array_INMDA_PYso_PYdr_sNMDA.size();
double* const _array_INMDA_PYso_PYdr_spike_activity = _dynamic_array_INMDA_PYso_PYdr_spike_activity.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_spike_activity[0];
const size_t _numspike_activity = _dynamic_array_INMDA_PYso_PYdr_spike_activity.size();
double* const _array_INMDA_PYso_PYdr_tauRes_NMDA = _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA[0];
const size_t _numtauRes_NMDA = _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.size();
double* const _array_INMDA_PYso_PYdr_tauS_NMDA = _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_tauS_NMDA[0];
const size_t _numtauS_NMDA = _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.size();
double* const _array_INMDA_PYso_PYdr_tauX_NMDA = _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_tauX_NMDA[0];
const size_t _numtauX_NMDA = _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.size();
const size_t _numv_pre = 100;
double* const _array_INMDA_PYso_PYdr_xNMDA = _dynamic_array_INMDA_PYso_PYdr_xNMDA.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_xNMDA[0];
const size_t _numxNMDA = _dynamic_array_INMDA_PYso_PYdr_xNMDA.size();
int32_t* const _array_INMDA_PYso_PYdr__synaptic_pre = _dynamic_array_INMDA_PYso_PYdr__synaptic_pre.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr__synaptic_pre[0];
const size_t _num_presynaptic_idx = _dynamic_array_INMDA_PYso_PYdr__synaptic_pre.size();
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_INMDA_PYso_PYdr_N = _array_INMDA_PYso_PYdr_N;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_alphaS_NMDA = _array_INMDA_PYso_PYdr_alphaS_NMDA;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_alphaX_NMDA = _array_INMDA_PYso_PYdr_alphaX_NMDA;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_deprFactor = _array_INMDA_PYso_PYdr_deprFactor;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_res_NMDA = _array_INMDA_PYso_PYdr_res_NMDA;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_sNMDA = _array_INMDA_PYso_PYdr_sNMDA;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_spike_activity = _array_INMDA_PYso_PYdr_spike_activity;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_tauRes_NMDA = _array_INMDA_PYso_PYdr_tauRes_NMDA;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_tauS_NMDA = _array_INMDA_PYso_PYdr_tauS_NMDA;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_tauX_NMDA = _array_INMDA_PYso_PYdr_tauX_NMDA;
    double* __restrict  _ptr_array_PYso_group_v = _array_PYso_group_v;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_xNMDA = _array_INMDA_PYso_PYdr_xNMDA;
    int32_t* __restrict  _ptr_array_INMDA_PYso_PYdr__synaptic_pre = _array_INMDA_PYso_PYdr__synaptic_pre;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = 1.0f*1.0/dt;
    const double _lio_2 = 1.0f*0.5/mV;
    const double _lio_3 = 20.0 * mV;


    const int _N = _array_INMDA_PYso_PYdr_N[0];
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _presynaptic_idx = _ptr_array_INMDA_PYso_PYdr__synaptic_pre[_idx];
        const double alphaS_NMDA = _ptr_array_INMDA_PYso_PYdr_alphaS_NMDA[_idx];
        const double alphaX_NMDA = _ptr_array_INMDA_PYso_PYdr_alphaX_NMDA[_idx];
        const double deprFactor = _ptr_array_INMDA_PYso_PYdr_deprFactor[_idx];
        double res_NMDA = _ptr_array_INMDA_PYso_PYdr_res_NMDA[_idx];
        double sNMDA = _ptr_array_INMDA_PYso_PYdr_sNMDA[_idx];
        const double spike_activity = _ptr_array_INMDA_PYso_PYdr_spike_activity[_idx];
        const double tauRes_NMDA = _ptr_array_INMDA_PYso_PYdr_tauRes_NMDA[_idx];
        const double tauS_NMDA = _ptr_array_INMDA_PYso_PYdr_tauS_NMDA[_idx];
        const double tauX_NMDA = _ptr_array_INMDA_PYso_PYdr_tauX_NMDA[_idx];
        const double v_pre = _ptr_array_PYso_group_v[_presynaptic_idx];
        double xNMDA = _ptr_array_INMDA_PYso_PYdr_xNMDA[_idx];
        const double _res_NMDA = (dt * ((spike_activity * ((1.0f*((- 1.0) + res_NMDA)/tauRes_NMDA) + (_lio_1 * ((deprFactor * res_NMDA) - res_NMDA)))) + (1.0f*(1.0 - res_NMDA)/tauRes_NMDA))) + res_NMDA;
        const double _sNMDA = (dt * (((alphaS_NMDA * xNMDA) * (1.0 - sNMDA)) - (1.0f*sNMDA/tauS_NMDA))) + sNMDA;
        const double _xNMDA = (dt * ((1.0f*alphaX_NMDA/(1.0 + exp(_lio_2 * (_lio_3 - v_pre)))) - (1.0f*xNMDA/tauX_NMDA))) + xNMDA;
        res_NMDA = _res_NMDA;
        sNMDA = _sNMDA;
        xNMDA = _xNMDA;
        _ptr_array_INMDA_PYso_PYdr_res_NMDA[_idx] = res_NMDA;
        _ptr_array_INMDA_PYso_PYdr_sNMDA[_idx] = sNMDA;
        _ptr_array_INMDA_PYso_PYdr_xNMDA[_idx] = xNMDA;

    }

}


