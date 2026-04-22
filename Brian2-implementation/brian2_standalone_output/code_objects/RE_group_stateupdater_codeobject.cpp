#include "code_objects/RE_group_stateupdater_codeobject.h"
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



void _run_RE_group_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double Cm_RE = 0.009999999999999998;
const double EKLeak_RE = - 0.095;
const double EK_RE = - 0.1;
const double ELeak_RE = - 0.09;
const double ENa_RE = 0.05;
const double ET_RE = 0.12;
const double Eext_RE = 0.0;
const int64_t N = 20;
const size_t _numdt = 1;
const double gKLeak_RE = 0.0;
const double gK_RE = 200.0;
const double gLeak_RE = 0.5;
const double gNa_RE = 2000.0;
const size_t _numgPoisson_RE = 20;
const double gT_RE = 30.0;
const size_t _numhNa_RE = 20;
const size_t _numhT_RE = 20;
const size_t _numiAMPA_PYso_RE = 20;
const size_t _numiAMPA_TC_RE = 20;
const size_t _numiGABAA_RE_RE = 20;
const size_t _nummNa_RE = 20;
const size_t _nummT_RE = 20;
const double mV = 0.001;
const double ms = 0.001;
const size_t _numnK_RE = 20;
const size_t _numnot_refractory = 20;
const double phiH_RE = 3.73;
const double phiM_RE = 6.81;
const double tau_RE = 0.002;
const size_t _numv = 20;
const int64_t vShiftK_RE = - 55;
const int64_t vShiftNa_RE = - 55;
const int64_t vShiftT_RE = 4;
    ///// POINTERS ////////////
        
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_RE_group_gPoisson_RE = _array_RE_group_gPoisson_RE;
    double* __restrict  _ptr_array_RE_group_hNa_RE = _array_RE_group_hNa_RE;
    double* __restrict  _ptr_array_RE_group_hT_RE = _array_RE_group_hT_RE;
    double* __restrict  _ptr_array_RE_group_iAMPA_PYso_RE = _array_RE_group_iAMPA_PYso_RE;
    double* __restrict  _ptr_array_RE_group_iAMPA_TC_RE = _array_RE_group_iAMPA_TC_RE;
    double* __restrict  _ptr_array_RE_group_iGABAA_RE_RE = _array_RE_group_iGABAA_RE_RE;
    double* __restrict  _ptr_array_RE_group_mNa_RE = _array_RE_group_mNa_RE;
    double* __restrict  _ptr_array_RE_group_mT_RE = _array_RE_group_mT_RE;
    double* __restrict  _ptr_array_RE_group_nK_RE = _array_RE_group_nK_RE;
    char* __restrict  _ptr_array_RE_group_not_refractory = _array_RE_group_not_refractory;
    double* __restrict  _ptr_array_RE_group_v = _array_RE_group_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = (- 25.0) * mV;
    const double _lio_2 = 1.0f*(- dt)/tau_RE;
    const double _lio_3 = 1.0f*(- 4.0)/ms;
    const double _lio_4 = 8.0 + (0.2 * vShiftNa_RE);
    const double _lio_5 = 1.0f*0.2/mV;
    const double _lio_6 = 1.0f*1.0/ms;
    const double _lio_7 = (0.05555555555555555 * vShiftNa_RE) + 0.9444444444444444;
    const double _lio_8 = 1.0f*0.05555555555555555/mV;
    const double _lio_9 = 1.0f*(dt * phiH_RE)/ms;
    const double _lio_10 = (0.2 * vShiftT_RE) + 15.6;
    const double _lio_11 = (0.02 * (- vShiftT_RE)) - 8.1;
    const double _lio_12 = 1.0f*0.02/mV;
    const double _lio_13 = (0.25 * vShiftT_RE) + 11.5;
    const double _lio_14 = 1.0f*0.25/mV;
    const double _lio_15 = (- 11.2) + ((- 0.28) * vShiftNa_RE);
    const double _lio_16 = 1.0f*0.28/mV;
    const double _lio_17 = (- 8.0) + (0.2 * (- vShiftNa_RE));
    const double _lio_18 = 4.16 + (0.32 * vShiftNa_RE);
    const double _lio_19 = 1.0f*0.32/mV;
    const double _lio_20 = (0.25 * vShiftNa_RE) + 3.25;
    const double _lio_21 = 1.0f*(dt * phiM_RE)/ms;
    const double _lio_22 = (- 0.135135135135135) * vShiftT_RE;
    const double _lio_23 = 1.0f*0.135135135135135/mV;
    const double _lio_24 = (0.06666666666666667 * (- vShiftT_RE)) - 6.666666666666667;
    const double _lio_25 = 1.0f*0.06666666666666667/mV;
    const double _lio_26 = (0.1 * vShiftT_RE) + 2.5;
    const double _lio_27 = 1.0f*0.1/mV;
    const double _lio_28 = 1.0f*(- 0.5)/ms;
    const double _lio_29 = (0.025 * vShiftK_RE) + 0.25;
    const double _lio_30 = 1.0f*0.025/mV;
    const double _lio_31 = 0.48 + (0.032 * vShiftK_RE);
    const double _lio_32 = 1.0f*0.032/mV;
    const double _lio_33 = 3.0 + (0.2 * vShiftK_RE);
    const double _lio_34 = 1.0f*dt/Cm_RE;
    const double _lio_35 = - gKLeak_RE;
    const double _lio_36 = - EKLeak_RE;
    const double _lio_37 = - EK_RE;
    const double _lio_38 = - ELeak_RE;
    const double _lio_39 = - ENa_RE;
    const double _lio_40 = - Eext_RE;
    const double _lio_41 = - ET_RE;


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double gPoisson_RE = _ptr_array_RE_group_gPoisson_RE[_idx];
        double hNa_RE = _ptr_array_RE_group_hNa_RE[_idx];
        double hT_RE = _ptr_array_RE_group_hT_RE[_idx];
        const double iAMPA_PYso_RE = _ptr_array_RE_group_iAMPA_PYso_RE[_idx];
        const double iAMPA_TC_RE = _ptr_array_RE_group_iAMPA_TC_RE[_idx];
        const double iGABAA_RE_RE = _ptr_array_RE_group_iGABAA_RE_RE[_idx];
        double mNa_RE = _ptr_array_RE_group_mNa_RE[_idx];
        double mT_RE = _ptr_array_RE_group_mT_RE[_idx];
        double nK_RE = _ptr_array_RE_group_nK_RE[_idx];
        char not_refractory = _ptr_array_RE_group_not_refractory[_idx];
        double v = _ptr_array_RE_group_v[_idx];
        if(!not_refractory)
            not_refractory = false || (! (v > _lio_1));
        else 
            not_refractory = true || (! (v > _lio_1));
        const double _gPoisson_RE = (_lio_2 * gPoisson_RE) + gPoisson_RE;
        const double _hNa_RE = (dt * ((1.0f*(_lio_3 * hNa_RE)/(1.0 + exp(_lio_4 - (_lio_5 * v)))) + (_lio_6 * ((0.128 - (0.128 * hNa_RE)) * exp(_lio_7 - (_lio_8 * v)))))) + hNa_RE;
        const double _hT_RE = (1.0f*(_lio_9 * ((- hT_RE) + (1.0f*1.0/(1.0 + exp(_lio_10 + (_lio_5 * v))))))/(85.0 + (1.0f*1.0/(exp(_lio_11 - (_lio_12 * v)) + exp(_lio_13 + (_lio_14 * v)))))) + hT_RE;
        const double _mNa_RE = (dt * ((1.0f*(_lio_6 * ((- mNa_RE) * (_lio_15 + (_lio_16 * v))))/((- 1.0) + exp(_lio_17 + (_lio_5 * v)))) + (1.0f*(_lio_6 * ((1.0 - mNa_RE) * (_lio_18 - (_lio_19 * v))))/((- 1.0) + exp(_lio_20 - (_lio_14 * v)))))) + mNa_RE;
        const double _mT_RE = (1.0f*(_lio_21 * ((- mT_RE) + (1.0f*1.0/(1.0 + (0.00116299493943615 * exp(_lio_22 - (_lio_23 * v)))))))/(3.0 + (1.0f*1.0/(exp(_lio_24 - (_lio_25 * v)) + exp(_lio_26 + (_lio_27 * v)))))) + mT_RE;
        const double _nK_RE = (dt * ((_lio_28 * (nK_RE * exp(_lio_29 - (_lio_30 * v)))) + (1.0f*(_lio_6 * ((1.0 - nK_RE) * (_lio_31 - (_lio_32 * v))))/((- 1.0) + exp(_lio_33 - (_lio_5 * v)))))) + nK_RE;
        const double _v = v + (_lio_34 * (((((_lio_35 * (_lio_36 + v)) + iAMPA_PYso_RE) + iAMPA_TC_RE) + iGABAA_RE_RE) - (((((gK_RE * ((_brian_pow(nK_RE, 4)) * (_lio_37 + v))) + (gLeak_RE * (_lio_38 + v))) + (gNa_RE * ((hNa_RE * (_brian_pow(mNa_RE, 3))) * (_lio_39 + v)))) + (gPoisson_RE * (_lio_40 + v))) + (gT_RE * ((hT_RE * (_brian_pow(mT_RE, 2))) * (_lio_41 + v))))));
        gPoisson_RE = _gPoisson_RE;
        hNa_RE = _hNa_RE;
        hT_RE = _hT_RE;
        mNa_RE = _mNa_RE;
        mT_RE = _mT_RE;
        nK_RE = _nK_RE;
        v = _v;
        _ptr_array_RE_group_gPoisson_RE[_idx] = gPoisson_RE;
        _ptr_array_RE_group_hNa_RE[_idx] = hNa_RE;
        _ptr_array_RE_group_hT_RE[_idx] = hT_RE;
        _ptr_array_RE_group_mNa_RE[_idx] = mNa_RE;
        _ptr_array_RE_group_mT_RE[_idx] = mT_RE;
        _ptr_array_RE_group_nK_RE[_idx] = nK_RE;
        _ptr_array_RE_group_not_refractory[_idx] = not_refractory;
        _ptr_array_RE_group_v[_idx] = v;

    }

}


