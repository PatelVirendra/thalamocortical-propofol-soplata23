#include "code_objects/IKNa_PYso_PYdr_summed_variable_iKNa_post_codeobject.h"
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



void _run_IKNa_PYso_PYdr_summed_variable_iKNa_post_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const double EKNa = - 0.1;
const size_t _numN = 1;
const int64_t N_post = 100;
int32_t* const _array_IKNa_PYso_PYdr__synaptic_post = _dynamic_array_IKNa_PYso_PYdr__synaptic_post.empty()? 0 : &_dynamic_array_IKNa_PYso_PYdr__synaptic_post[0];
const size_t _num_synaptic_post = _dynamic_array_IKNa_PYso_PYdr__synaptic_post.size();
double* const _array_IKNa_PYso_PYdr_concNa = _dynamic_array_IKNa_PYso_PYdr_concNa.empty()? 0 : &_dynamic_array_IKNa_PYso_PYdr_concNa[0];
const size_t _numconcNa = _dynamic_array_IKNa_PYso_PYdr_concNa.size();
const double gKNa_val = 13.299999999999999;
const size_t _numiKNa_post = 100;
const double mM = 1.0;
const size_t _numv_post = 100;
const size_t _num_postsynaptic_idx = _dynamic_array_IKNa_PYso_PYdr__synaptic_post.size();
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_IKNa_PYso_PYdr_N = _array_IKNa_PYso_PYdr_N;
    int32_t* __restrict  _ptr_array_IKNa_PYso_PYdr__synaptic_post = _array_IKNa_PYso_PYdr__synaptic_post;
    double* __restrict  _ptr_array_IKNa_PYso_PYdr_concNa = _array_IKNa_PYso_PYdr_concNa;
    double* __restrict  _ptr_array_PYso_group_iKNa = _array_PYso_group_iKNa;
    double* __restrict  _ptr_array_PYso_group_v = _array_PYso_group_v;


    //// MAIN CODE ////////////
    const int _target_size = N_post;

    // Set all the target variable values to zero
    #pragma omp parallel for schedule(static)
    for (int _target_idx=0; _target_idx<_target_size; _target_idx++)
    {
        _ptr_array_PYso_group_iKNa[_target_idx + 0] = 0;
    }

    // scalar code
    const size_t _vectorisation_idx = -1;
        
    const double _lio_1 = 0.37 * (- gKNa_val);
    const double _lio_2 = 0.0 - EKNa;
    const double _lio_3 = 38.7 * mM;


    for(int _idx=0; _idx<_ptr_array_IKNa_PYso_PYdr_N[0]; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _postsynaptic_idx = _ptr_array_IKNa_PYso_PYdr__synaptic_post[_idx];
        const double concNa = _ptr_array_IKNa_PYso_PYdr_concNa[_idx];
        const double v_post = _ptr_array_PYso_group_v[_postsynaptic_idx];
        const double _synaptic_var = 1.0f*(_lio_1 * (_lio_2 + v_post))/(1.0 + (_brian_pow(1.0f*_lio_3/concNa, 3.5)));

        _ptr_array_PYso_group_iKNa[_ptr_array_IKNa_PYso_PYdr__synaptic_post[_idx]] += _synaptic_var;
    }

}


