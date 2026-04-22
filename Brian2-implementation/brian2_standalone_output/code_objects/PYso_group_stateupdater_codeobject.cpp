#include "code_objects/PYso_group_stateupdater_codeobject.h"
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



void _run_PYso_group_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double Cm_PYso = 0.009999999999999998;
const double EA_PYso = - 0.1;
const double EKS_PYso = - 0.1;
const double EK_PYso = - 0.1;
const double ELeak_PYso = - 0.060950000000000004;
const double ENa_PYso = 0.055;
const int64_t N = 100;
const size_t _numdt = 1;
const double gA_PYso = 10.0;
const double gKS_PYso = 5.76;
const double gK_PYso = 105.0;
const double gLeak_PYso = 0.6669999999999999;
const double gNa_PYso = 500.0;
const size_t _numhA = 100;
const size_t _numhNa = 100;
const size_t _numiCOM = 100;
const size_t _numiGABAA_IN_PYso = 100;
const size_t _numiKNa = 100;
const size_t _nummKS = 100;
const double mV = 0.001;
const double ms = 0.001;
const size_t _numnK = 100;
const size_t _numnot_refractory = 100;
const int64_t phi_A_PYso = 1;
const int64_t phi_KS_PYso = 1;
const int64_t phi_K_PYso = 4;
const int64_t phi_Na_PYso = 4;
const double tauH_A_PYso = 0.015;
const size_t _numv = 100;
    ///// POINTERS ////////////
        
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_PYso_group_hA = _array_PYso_group_hA;
    double* __restrict  _ptr_array_PYso_group_hNa = _array_PYso_group_hNa;
    double* __restrict  _ptr_array_PYso_group_iCOM = _array_PYso_group_iCOM;
    double* __restrict  _ptr_array_PYso_group_iGABAA_IN_PYso = _array_PYso_group_iGABAA_IN_PYso;
    double* __restrict  _ptr_array_PYso_group_iKNa = _array_PYso_group_iKNa;
    double* __restrict  _ptr_array_PYso_group_mKS = _array_PYso_group_mKS;
    double* __restrict  _ptr_array_PYso_group_nK = _array_PYso_group_nK;
    char* __restrict  _ptr_array_PYso_group_not_refractory = _array_PYso_group_not_refractory;
    double* __restrict  _ptr_array_PYso_group_v = _array_PYso_group_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = (- 25.0) * mV;
    const double _lio_2 = 1.0f*(dt * phi_A_PYso)/tauH_A_PYso;
    const double _lio_3 = 1.0f*0.16666666666666666/mV;
    const double _lio_4 = dt * phi_Na_PYso;
    const double _lio_5 = 1.0f*1.0/ms;
    const double _lio_6 = 1.0f*0.1/mV;
    const double _lio_7 = 1.0f*(0.125 * (dt * phi_KS_PYso))/ms;
    const double _lio_8 = 1.0f*(- 0.153846153846154)/mV;
    const double _lio_9 = 1.0f*0.03333333333333333/mV;
    const double _lio_10 = dt * phi_K_PYso;
    const double _lio_11 = 1.0f*(- 0.125)/ms;
    const double _lio_12 = 1.0f*0.04/mV;
    const double _lio_13 = 1.0f*0.01/mV;
    const double _lio_14 = 1.0f*dt/Cm_PYso;
    const double _lio_15 = - gA_PYso;
    const double _lio_16 = - EA_PYso;
    const double _lio_17 = 1.0f*0.05/mV;
    const double _lio_18 = - EKS_PYso;
    const double _lio_19 = - EK_PYso;
    const double _lio_20 = - ELeak_PYso;
    const double _lio_21 = 1.0f*gNa_PYso/(_brian_pow(ms, 3));
    const double _lio_22 = - ENa_PYso;
    const double _lio_23 = 1.0f*0.0455608884980535/ms;
    const double _lio_24 = 1.0f*0.08333333333333333/mV;


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double hA = _ptr_array_PYso_group_hA[_idx];
        double hNa = _ptr_array_PYso_group_hNa[_idx];
        const double iCOM = _ptr_array_PYso_group_iCOM[_idx];
        const double iGABAA_IN_PYso = _ptr_array_PYso_group_iGABAA_IN_PYso[_idx];
        const double iKNa = _ptr_array_PYso_group_iKNa[_idx];
        double mKS = _ptr_array_PYso_group_mKS[_idx];
        double nK = _ptr_array_PYso_group_nK[_idx];
        char not_refractory = _ptr_array_PYso_group_not_refractory[_idx];
        double v = _ptr_array_PYso_group_v[_idx];
        if(!not_refractory)
            not_refractory = false || (! (v > _lio_1));
        else 
            not_refractory = true || (! (v > _lio_1));
        const double _hA = (_lio_2 * ((- hA) + (1.0f*1.0/(1.0 + exp(13.333333333333334 + (_lio_3 * v)))))) + hA;
        const double _hNa = (_lio_4 * ((1.0f*(_lio_5 * (- hNa))/(1.0 + exp((- 2.0) - (_lio_6 * v)))) + (_lio_5 * ((0.07 - (0.07 * hNa)) * exp((- 5.0) - (_lio_6 * v)))))) + hNa;
        const double _mKS = (_lio_7 * (((- mKS) + (1.0f*1.0/(1.0 + (0.00534940878975042 * exp(_lio_8 * v))))) * (exp((- 1.8333333333333333) - (_lio_9 * v)) + exp(1.8333333333333333 + (_lio_9 * v))))) + mKS;
        const double _nK = (_lio_10 * ((_lio_11 * (nK * exp((- 1.76) - (_lio_12 * v)))) + (1.0f*(_lio_5 * ((0.34 + (_lio_13 * v)) * (1.0 - nK)))/(1.0 - exp((- 3.4000000000000004) - (_lio_6 * v)))))) + nK;
        const double _v = v + (_lio_14 * (((((1.0f*(_lio_15 * (hA * (_lio_16 + v)))/(_brian_pow(1.0 + exp((- 2.5) - (_lio_17 * v)), 3))) + iCOM) + iGABAA_IN_PYso) + iKNa) - ((((gKS_PYso * ((_brian_pow(mKS, 3)) * (_lio_18 + v))) + (gK_PYso * ((_brian_pow(nK, 4)) * (_lio_19 + v)))) + (gLeak_PYso * (_lio_20 + v))) + (1.0f*(_lio_21 * ((hNa * (_brian_pow(3.3 + (_lio_6 * v), 3))) * (_lio_22 + v)))/((_brian_pow(1.0 - exp((- 3.3000000000000003) - (_lio_6 * v)), 3)) * (_brian_pow((_lio_23 * exp(_lio_24 * (- v))) + (1.0f*(_lio_5 * (3.3 + (_lio_6 * v)))/(1.0 - exp((- 3.3000000000000003) - (_lio_6 * v)))), 3)))))));
        hA = _hA;
        hNa = _hNa;
        mKS = _mKS;
        nK = _nK;
        v = _v;
        _ptr_array_PYso_group_hA[_idx] = hA;
        _ptr_array_PYso_group_hNa[_idx] = hNa;
        _ptr_array_PYso_group_mKS[_idx] = mKS;
        _ptr_array_PYso_group_nK[_idx] = nK;
        _ptr_array_PYso_group_not_refractory[_idx] = not_refractory;
        _ptr_array_PYso_group_v[_idx] = v;

    }

}


