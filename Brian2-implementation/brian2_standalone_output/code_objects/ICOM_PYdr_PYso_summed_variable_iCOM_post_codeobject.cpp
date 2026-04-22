#include "code_objects/ICOM_PYdr_PYso_summed_variable_iCOM_post_codeobject.h"
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



void _run_ICOM_PYdr_PYso_summed_variable_iCOM_post_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    const size_t _numN = 1;
const int64_t N_post = 100;
int32_t* const _array_ICOM_PYdr_PYso__synaptic_post = _dynamic_array_ICOM_PYdr_PYso__synaptic_post.empty()? 0 : &_dynamic_array_ICOM_PYdr_PYso__synaptic_post[0];
const size_t _num_synaptic_post = _dynamic_array_ICOM_PYdr_PYso__synaptic_post.size();
double* const _array_ICOM_PYdr_PYso_gCOM = _dynamic_array_ICOM_PYdr_PYso_gCOM.empty()? 0 : &_dynamic_array_ICOM_PYdr_PYso_gCOM[0];
const size_t _numgCOM = _dynamic_array_ICOM_PYdr_PYso_gCOM.size();
const size_t _numiCOM_post = 100;
const size_t _numv_post = 100;
const size_t _numv_pre = 100;
const size_t _num_postsynaptic_idx = _dynamic_array_ICOM_PYdr_PYso__synaptic_post.size();
int32_t* const _array_ICOM_PYdr_PYso__synaptic_pre = _dynamic_array_ICOM_PYdr_PYso__synaptic_pre.empty()? 0 : &_dynamic_array_ICOM_PYdr_PYso__synaptic_pre[0];
const size_t _num_presynaptic_idx = _dynamic_array_ICOM_PYdr_PYso__synaptic_pre.size();
    ///// POINTERS ////////////
        
    int32_t*   _ptr_array_ICOM_PYdr_PYso_N = _array_ICOM_PYdr_PYso_N;
    int32_t* __restrict  _ptr_array_ICOM_PYdr_PYso__synaptic_post = _array_ICOM_PYdr_PYso__synaptic_post;
    double* __restrict  _ptr_array_ICOM_PYdr_PYso_gCOM = _array_ICOM_PYdr_PYso_gCOM;
    double* __restrict  _ptr_array_PYso_group_iCOM = _array_PYso_group_iCOM;
    double* __restrict  _ptr_array_PYso_group_v = _array_PYso_group_v;
    double* __restrict  _ptr_array_PYdr_group_v = _array_PYdr_group_v;
    int32_t* __restrict  _ptr_array_ICOM_PYdr_PYso__synaptic_pre = _array_ICOM_PYdr_PYso__synaptic_pre;


    //// MAIN CODE ////////////
    const int _target_size = N_post;

    // Set all the target variable values to zero
    #pragma omp parallel for schedule(static)
    for (int _target_idx=0; _target_idx<_target_size; _target_idx++)
    {
        _ptr_array_PYso_group_iCOM[_target_idx + 0] = 0;
    }

    // scalar code
    const size_t _vectorisation_idx = -1;
        


    for(int _idx=0; _idx<_ptr_array_ICOM_PYdr_PYso_N[0]; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _postsynaptic_idx = _ptr_array_ICOM_PYdr_PYso__synaptic_post[_idx];
        const int32_t _presynaptic_idx = _ptr_array_ICOM_PYdr_PYso__synaptic_pre[_idx];
        const double gCOM = _ptr_array_ICOM_PYdr_PYso_gCOM[_idx];
        const double v_post = _ptr_array_PYso_group_v[_postsynaptic_idx];
        const double v_pre = _ptr_array_PYdr_group_v[_presynaptic_idx];
        const double _synaptic_var = (- gCOM) * (v_post - v_pre);

        _ptr_array_PYso_group_iCOM[_ptr_array_ICOM_PYdr_PYso__synaptic_post[_idx]] += _synaptic_var;
    }

}


