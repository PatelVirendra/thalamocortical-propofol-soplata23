#include "code_objects/TC_group_stateupdater_codeobject.h"
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



void _run_TC_group_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numCa_TC = 20;
const double Ca_inf_TC = 0.00024;
const double Cac_TC = 0.002;
const double Cm_TC = 0.009999999999999998;
const double EH_TC = - 0.04;
const double EKLeak_TC = - 0.1;
const double EK_TC = - 0.1;
const double ELeak_TC = - 0.07;
const double ENa_TC = 0.05;
const double Eext_TC = 0.0;
const int64_t N = 20;
const size_t _numOpen = 20;
const size_t _numOpenLocked = 20;
const size_t _numPone = 20;
const size_t _numdCa_extra = 20;
const size_t _numdt = 1;
const double gH_TC = 0.05;
const int64_t gInc_TC = 2;
const double gKLeak_TC = 0.17200000000000001;
const double gK_TC = 100.0;
const double gLeak_TC = 0.1;
const double gNa_TC = 899.9999999999999;
const size_t _numgPoisson_TC = 20;
const double gT_TC = 20.0;
const size_t _numhNa_TC = 20;
const size_t _numhT_TC = 20;
const size_t _numiAMPA_PYso_TC = 20;
const size_t _numiGABAA_RE_TC = 20;
const size_t _numiGABAB_RE_TC = 20;
const double k2_TC = 0.4;
const double k4_TC = 1.0;
const double mM = 1.0;
const size_t _nummNa_TC = 20;
const double mV = 0.001;
const double ms = 0.001;
const size_t _numnK_TC = 20;
const int64_t nca_TC = 4;
const int64_t nexp_TC = 1;
const size_t _numnot_refractory = 20;
const double pc_TC = 0.007;
const double phiH_TC = 3.73;
const double tauR_TC = 0.005;
const double tau_TC = 0.002;
const size_t _numv = 20;
const int64_t vShiftK_TC = - 25;
const int64_t vShiftNa_TC = - 40;
const int64_t vShiftT_TC = 2;
    ///// POINTERS ////////////
        
    double* __restrict  _ptr_array_TC_group_Ca_TC = _array_TC_group_Ca_TC;
    double* __restrict  _ptr_array_TC_group_Open = _array_TC_group_Open;
    double* __restrict  _ptr_array_TC_group_OpenLocked = _array_TC_group_OpenLocked;
    double* __restrict  _ptr_array_TC_group_Pone = _array_TC_group_Pone;
    double* __restrict  _ptr_array_TC_group_dCa_extra = _array_TC_group_dCa_extra;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_TC_group_gPoisson_TC = _array_TC_group_gPoisson_TC;
    double* __restrict  _ptr_array_TC_group_hNa_TC = _array_TC_group_hNa_TC;
    double* __restrict  _ptr_array_TC_group_hT_TC = _array_TC_group_hT_TC;
    double* __restrict  _ptr_array_TC_group_iAMPA_PYso_TC = _array_TC_group_iAMPA_PYso_TC;
    double* __restrict  _ptr_array_TC_group_iGABAA_RE_TC = _array_TC_group_iGABAA_RE_TC;
    double* __restrict  _ptr_array_TC_group_iGABAB_RE_TC = _array_TC_group_iGABAB_RE_TC;
    double* __restrict  _ptr_array_TC_group_mNa_TC = _array_TC_group_mNa_TC;
    double* __restrict  _ptr_array_TC_group_nK_TC = _array_TC_group_nK_TC;
    char* __restrict  _ptr_array_TC_group_not_refractory = _array_TC_group_not_refractory;
    double* __restrict  _ptr_array_TC_group_v = _array_TC_group_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = (- 25.0) * mV;
    const double _lio_2 = 1.0f*1.0/tauR_TC;
    const double _lio_3 = 1.0f*1.0/ms;
    const double _lio_4 = 1.0f*0.181818181818182/mV;
    const double _lio_5 = 1.0f*0.0704225352112676/mV;
    const double _lio_6 = 1.0f*(- 0.0862068965517241)/mV;
    const double _lio_7 = 1.0f*1.0/pc_TC;
    const double _lio_8 = 1.0f*1.0/Cac_TC;
    const double _lio_9 = 1.0f*(- dt)/tau_TC;
    const double _lio_10 = 1.0f*(- 4.0)/ms;
    const double _lio_11 = 8.0 + (0.2 * vShiftNa_TC);
    const double _lio_12 = 1.0f*0.2/mV;
    const double _lio_13 = (0.05555555555555555 * vShiftNa_TC) + 0.9444444444444444;
    const double _lio_14 = 1.0f*0.05555555555555555/mV;
    const double _lio_15 = 1.0f*(dt * phiH_TC)/ms;
    const double _lio_16 = (0.25 * vShiftT_TC) + 20.25;
    const double _lio_17 = 1.0f*0.25/mV;
    const double _lio_18 = 0.2 * vShiftT_TC;
    const double _lio_19 = 0.3125 * vShiftT_TC;
    const double _lio_20 = 1.0f*0.3125/mV;
    const double _lio_21 = (- 11.2) + ((- 0.28) * vShiftNa_TC);
    const double _lio_22 = 1.0f*0.28/mV;
    const double _lio_23 = (- 8.0) + (0.2 * (- vShiftNa_TC));
    const double _lio_24 = 4.16 + (0.32 * vShiftNa_TC);
    const double _lio_25 = 1.0f*0.32/mV;
    const double _lio_26 = (0.25 * vShiftNa_TC) + 3.25;
    const double _lio_27 = 1.0f*(- 0.5)/ms;
    const double _lio_28 = (0.025 * vShiftK_TC) + 0.25;
    const double _lio_29 = 1.0f*0.025/mV;
    const double _lio_30 = 0.48 + (0.032 * vShiftK_TC);
    const double _lio_31 = 1.0f*0.032/mV;
    const double _lio_32 = 3.0 + (0.2 * vShiftK_TC);
    const double _lio_33 = 1.0f*dt/Cm_TC;
    const double _lio_34 = - gH_TC;
    const double _lio_35 = - EH_TC;
    const double _lio_36 = - EKLeak_TC;
    const double _lio_37 = - EK_TC;
    const double _lio_38 = - ELeak_TC;
    const double _lio_39 = - ENa_TC;
    const double _lio_40 = - Eext_TC;
    const double _lio_41 = (- 13.2705524828078) * mV;
    const double _lio_42 = 2.0 * mM;
    const double _lio_43 = (- 0.161290322580645) * vShiftT_TC;
    const double _lio_44 = 1.0f*0.161290322580645/mV;


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double Ca_TC = _ptr_array_TC_group_Ca_TC[_idx];
        double Open = _ptr_array_TC_group_Open[_idx];
        double OpenLocked = _ptr_array_TC_group_OpenLocked[_idx];
        double Pone = _ptr_array_TC_group_Pone[_idx];
        const double dCa_extra = _ptr_array_TC_group_dCa_extra[_idx];
        double gPoisson_TC = _ptr_array_TC_group_gPoisson_TC[_idx];
        double hNa_TC = _ptr_array_TC_group_hNa_TC[_idx];
        double hT_TC = _ptr_array_TC_group_hT_TC[_idx];
        const double iAMPA_PYso_TC = _ptr_array_TC_group_iAMPA_PYso_TC[_idx];
        const double iGABAA_RE_TC = _ptr_array_TC_group_iGABAA_RE_TC[_idx];
        const double iGABAB_RE_TC = _ptr_array_TC_group_iGABAB_RE_TC[_idx];
        double mNa_TC = _ptr_array_TC_group_mNa_TC[_idx];
        double nK_TC = _ptr_array_TC_group_nK_TC[_idx];
        char not_refractory = _ptr_array_TC_group_not_refractory[_idx];
        double v = _ptr_array_TC_group_v[_idx];
        if(!not_refractory)
            not_refractory = false || (! (v > _lio_1));
        else 
            not_refractory = true || (! (v > _lio_1));
        const double _Ca_TC = Ca_TC + (dt * (dCa_extra + (_lio_2 * (Ca_inf_TC + (- Ca_TC)))));
        const double _Open = Open + (dt * ((1.0f*(_lio_3 * ((- Open) * (1.0 - (1.0f*1.0/(1.0 + (835983.066403625 * exp(_lio_4 * v)))))))/(20.0 + (1.0f*1000.0/((153.732067787116 * exp(_lio_5 * v)) + (0.000465492863075605 * exp(_lio_6 * v)))))) + (1.0f*(_lio_3 * ((1.0 + (- Open)) - OpenLocked))/((20.0 + (1.0f*1000.0/((153.732067787116 * exp(_lio_5 * v)) + (0.000465492863075605 * exp(_lio_6 * v))))) * (1.0 + (835983.066403625 * exp(_lio_4 * v)))))));
        const double _OpenLocked = OpenLocked + (dt * ((k4_TC * (Open * (_brian_pow(_lio_7 * Pone, nexp_TC)))) - (k4_TC * OpenLocked)));
        const double _Pone = Pone + (dt * ((k2_TC * (- Pone)) + (k2_TC * ((_brian_pow(_lio_8 * Ca_TC, nca_TC)) * (1.0 - Pone)))));
        const double _gPoisson_TC = (_lio_9 * gPoisson_TC) + gPoisson_TC;
        const double _hNa_TC = (dt * ((1.0f*(_lio_10 * hNa_TC)/(1.0 + exp(_lio_11 - (_lio_12 * v)))) + (_lio_3 * ((0.128 - (0.128 * hNa_TC)) * exp(_lio_13 - (_lio_14 * v)))))) + hNa_TC;
        const double _hT_TC = (1.0f*(_lio_15 * ((- hT_TC) + (1.0f*1.0/(1.0 + exp(_lio_16 + (_lio_17 * v))))))/(30.8 + (1.0f*(211.4 + (6798718666.66326 * exp(_lio_18 + (_lio_12 * v))))/(1.0 + (251321793304.994 * exp(_lio_19 + (_lio_20 * v))))))) + hT_TC;
        const double _mNa_TC = (dt * ((1.0f*(_lio_3 * ((- mNa_TC) * (_lio_21 + (_lio_22 * v))))/((- 1.0) + exp(_lio_23 + (_lio_12 * v)))) + (1.0f*(_lio_3 * ((1.0 - mNa_TC) * (_lio_24 - (_lio_25 * v))))/((- 1.0) + exp(_lio_26 - (_lio_17 * v)))))) + mNa_TC;
        const double _nK_TC = (dt * ((_lio_27 * (nK_TC * exp(_lio_28 - (_lio_29 * v)))) + (1.0f*(_lio_3 * ((1.0 - nK_TC) * (_lio_30 - (_lio_31 * v))))/((- 1.0) + exp(_lio_32 - (_lio_12 * v)))))) + nK_TC;
        const double _v = v + (_lio_33 * (((((_lio_34 * ((_lio_35 + v) * (Open + (gInc_TC * OpenLocked)))) + iAMPA_PYso_TC) + iGABAA_RE_TC) + iGABAB_RE_TC) - ((((((gKLeak_TC * (_lio_36 + v)) + (gK_TC * ((_brian_pow(nK_TC, 4)) * (_lio_37 + v)))) + (gLeak_TC * (_lio_38 + v))) + (gNa_TC * ((hNa_TC * (_brian_pow(mNa_TC, 3))) * (_lio_39 + v)))) + (gPoisson_TC * (_lio_40 + v))) + (1.0f*(gT_TC * (hT_TC * ((_lio_41 * log(1.0f*_lio_42/Ca_TC)) + v)))/(_brian_pow(1.0 + (0.000101693376272292 * exp(_lio_43 - (_lio_44 * v))), 2))))));
        Ca_TC = _Ca_TC;
        Open = _Open;
        OpenLocked = _OpenLocked;
        Pone = _Pone;
        gPoisson_TC = _gPoisson_TC;
        hNa_TC = _hNa_TC;
        hT_TC = _hT_TC;
        mNa_TC = _mNa_TC;
        nK_TC = _nK_TC;
        v = _v;
        _ptr_array_TC_group_Ca_TC[_idx] = Ca_TC;
        _ptr_array_TC_group_Open[_idx] = Open;
        _ptr_array_TC_group_OpenLocked[_idx] = OpenLocked;
        _ptr_array_TC_group_Pone[_idx] = Pone;
        _ptr_array_TC_group_gPoisson_TC[_idx] = gPoisson_TC;
        _ptr_array_TC_group_hNa_TC[_idx] = hNa_TC;
        _ptr_array_TC_group_hT_TC[_idx] = hT_TC;
        _ptr_array_TC_group_mNa_TC[_idx] = mNa_TC;
        _ptr_array_TC_group_nK_TC[_idx] = nK_TC;
        _ptr_array_TC_group_not_refractory[_idx] = not_refractory;
        _ptr_array_TC_group_v[_idx] = v;

    }

}


