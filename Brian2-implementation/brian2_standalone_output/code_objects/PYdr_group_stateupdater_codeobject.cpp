#include "code_objects/PYdr_group_stateupdater_codeobject.h"
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



void _run_PYdr_group_stateupdater_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numCaBuffer = 100;
const double Cm_PYdr = 0.009999999999999998;
const double EAR_PYdr = - 0.1;
const double EHVA_PYdr = 0.12;
const double EKCa_PYdr = - 0.1;
const double ELeak_PYdr = - 0.060950000000000004;
const double ENaP_PYdr = 0.055;
const double Eext_PYdr = 0.0;
const double KD_PYdr = 0.03;
const int64_t N = 100;
const double alphaCa_PYdr = 5000000.0;
const double areaDR = 3.5e-08;
const size_t _numdt = 1;
const double gAR_PYdr = 0.257;
const double gHVA_PYdr = 4.3;
const double gKCa_PYdr = 5.699999999999999;
const double gLeak_PYdr = 0.05;
const double gNaP_PYdr = 0.6859999999999999;
const size_t _numgPoisson_PYdr = 100;
const size_t _numiAMPA_PYso_PYdr = 100;
const size_t _numiAMPA_TC_PYdr = 100;
const size_t _numiCOM = 100;
const size_t _numiNMDA_PYso_PYdr = 100;
const double mV = 0.001;
const size_t _numnot_refractory = 100;
const double tauCa_PYdr = 0.15;
const double tau_PYdr = 0.002;
const size_t _numv = 100;
    ///// POINTERS ////////////
        
    double* __restrict  _ptr_array_PYdr_group_CaBuffer = _array_PYdr_group_CaBuffer;
    double*   _ptr_array_defaultclock_dt = _array_defaultclock_dt;
    double* __restrict  _ptr_array_PYdr_group_gPoisson_PYdr = _array_PYdr_group_gPoisson_PYdr;
    double* __restrict  _ptr_array_PYdr_group_iAMPA_PYso_PYdr = _array_PYdr_group_iAMPA_PYso_PYdr;
    double* __restrict  _ptr_array_PYdr_group_iAMPA_TC_PYdr = _array_PYdr_group_iAMPA_TC_PYdr;
    double* __restrict  _ptr_array_PYdr_group_iCOM = _array_PYdr_group_iCOM;
    double* __restrict  _ptr_array_PYdr_group_iNMDA_PYso_PYdr = _array_PYdr_group_iNMDA_PYso_PYdr;
    char* __restrict  _ptr_array_PYdr_group_not_refractory = _array_PYdr_group_not_refractory;
    double* __restrict  _ptr_array_PYdr_group_v = _array_PYdr_group_v;


    //// MAIN CODE ////////////
    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double dt = _ptr_array_defaultclock_dt[0];
    const double _lio_1 = (- 25.0) * mV;
    const double _lio_2 = 1.0f*1.0/tauCa_PYdr;
    const double _lio_3 = (alphaCa_PYdr * areaDR) * gHVA_PYdr;
    const double _lio_4 = - EHVA_PYdr;
    const double _lio_5 = 1.0f*0.1111111111111111/mV;
    const double _lio_6 = 1.0f*(- dt)/tau_PYdr;
    const double _lio_7 = 1.0f*dt/Cm_PYdr;
    const double _lio_8 = - EKCa_PYdr;
    const double _lio_9 = - EAR_PYdr;
    const double _lio_10 = 1.0f*0.25/mV;
    const double _lio_11 = - ELeak_PYdr;
    const double _lio_12 = - ENaP_PYdr;
    const double _lio_13 = 1.0f*(- 0.12987012987013)/mV;
    const double _lio_14 = - Eext_PYdr;


    const int _N = N;
    #pragma omp parallel for schedule(static)
    for(int _idx=0; _idx<_N; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        double CaBuffer = _ptr_array_PYdr_group_CaBuffer[_idx];
        double gPoisson_PYdr = _ptr_array_PYdr_group_gPoisson_PYdr[_idx];
        const double iAMPA_PYso_PYdr = _ptr_array_PYdr_group_iAMPA_PYso_PYdr[_idx];
        const double iAMPA_TC_PYdr = _ptr_array_PYdr_group_iAMPA_TC_PYdr[_idx];
        const double iCOM = _ptr_array_PYdr_group_iCOM[_idx];
        const double iNMDA_PYso_PYdr = _ptr_array_PYdr_group_iNMDA_PYso_PYdr[_idx];
        char not_refractory = _ptr_array_PYdr_group_not_refractory[_idx];
        double v = _ptr_array_PYdr_group_v[_idx];
        if(!not_refractory)
            not_refractory = false || (! (v > _lio_1));
        else 
            not_refractory = true || (! (v > _lio_1));
        const double _CaBuffer = CaBuffer + (dt * ((_lio_2 * (- CaBuffer)) - (1.0f*(_lio_3 * (_lio_4 + v))/(_brian_pow(1.0 + exp((- 2.2222222222222223) - (_lio_5 * v)), 2)))));
        const double _gPoisson_PYdr = (_lio_6 * gPoisson_PYdr) + gPoisson_PYdr;
        const double _v = v + (_lio_7 * ((((((1.0f*(gKCa_PYdr * ((- CaBuffer) * (_lio_8 + v)))/(KD_PYdr + CaBuffer)) + iAMPA_PYso_PYdr) + iAMPA_TC_PYdr) + iCOM) + iNMDA_PYso_PYdr) - (((((1.0f*(gAR_PYdr * (_lio_9 + v))/(1.0 + exp(18.75 + (_lio_10 * v)))) + (1.0f*(gHVA_PYdr * (_lio_4 + v))/(_brian_pow(1.0 + exp((- 2.2222222222222223) - (_lio_5 * v)), 2)))) + (gLeak_PYdr * (_lio_11 + v))) + (1.0f*(gNaP_PYdr * (_lio_12 + v))/(_brian_pow(1.0 + (0.000721797280255039 * exp(_lio_13 * v)), 3)))) + (gPoisson_PYdr * (_lio_14 + v)))));
        CaBuffer = _CaBuffer;
        gPoisson_PYdr = _gPoisson_PYdr;
        v = _v;
        _ptr_array_PYdr_group_CaBuffer[_idx] = CaBuffer;
        _ptr_array_PYdr_group_gPoisson_PYdr[_idx] = gPoisson_PYdr;
        _ptr_array_PYdr_group_not_refractory[_idx] = not_refractory;
        _ptr_array_PYdr_group_v[_idx] = v;

    }

}


