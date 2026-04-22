#include "code_objects/IGABAA_RE_TC_summed_variable_iGABAA_RE_TC_post_codeobject.h"
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



void _run_IGABAA_RE_TC_summed_variable_iGABAA_RE_TC_post_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    double* const _array_IGABAA_RE_TC_EGABAA = _dynamic_array_IGABAA_RE_TC_EGABAA.empty()? 0 : &_dynamic_array_IGABAA_RE_TC_EGABAA[0];
const size_t _numEGABAA = _dynamic_array_IGABAA_RE_TC_EGABAA.size();
const size_t _numN = 1;
const int64_t N_post = 20;
double* const _array_IGABAA_RE_TC_Npre = _dynamic_array_IGABAA_RE_TC_Npre.empty()? 0 : &_dynamic_array_IGABAA_RE_TC_Npre[0];
const size_t _numNpre = _dynamic_array_IGABAA_RE_TC_Npre.size();
int32_t* const _array_IGABAA_RE_TC__synaptic_post = _dynamic_array_IGABAA_RE_TC__synaptic_post.empty()? 0 : &_dynamic_array_IGABAA_RE_TC__synaptic_post[0];
const size_t _num_synaptic_post = _dynamic_array_IGABAA_RE_TC__synaptic_post.size();
double* const _array_IGABAA_RE_TC_gGABAA = _dynamic_array_IGABAA_RE_TC_gGABAA.empty()? 0 : &_dynamic_array_IGABAA_RE_TC_gGABAA[0];
const size_t _numgGABAA = _dynamic_array_IGABAA_RE_TC_gGABAA.size();
const size_t _numiGABAA_RE_TC_post = 20;
double* const _array_IGABAA_RE_TC_propoCondMult = _dynamic_array_IGABAA_RE_TC_propoCondMult.empty()? 0 : &_dynamic_array_IGABAA_RE_TC_propoCondMult[0];
const size_t _numpropoCondMult = _dynamic_array_IGABAA_RE_TC_propoCondMult.size();
double* const _array_IGABAA_RE_TC_sGABAA = _dynamic_array_IGABAA_RE_TC_sGABAA.empty()? 0 : &_dynamic_array_IGABAA_RE_TC_sGABAA[0];
const size_t _numsGABAA = _dynamic_array_IGABAA_RE_TC_sGABAA.size();
const size_t _numv_post = 20;
const size_t _num_postsynaptic_idx = _dynamic_array_IGABAA_RE_TC__synaptic_post.size();
    ///// POINTERS ////////////
        
    double* __restrict  _ptr_array_IGABAA_RE_TC_EGABAA = _array_IGABAA_RE_TC_EGABAA;
    int32_t*   _ptr_array_IGABAA_RE_TC_N = _array_IGABAA_RE_TC_N;
    double* __restrict  _ptr_array_IGABAA_RE_TC_Npre = _array_IGABAA_RE_TC_Npre;
    int32_t* __restrict  _ptr_array_IGABAA_RE_TC__synaptic_post = _array_IGABAA_RE_TC__synaptic_post;
    double* __restrict  _ptr_array_IGABAA_RE_TC_gGABAA = _array_IGABAA_RE_TC_gGABAA;
    double* __restrict  _ptr_array_TC_group_iGABAA_RE_TC = _array_TC_group_iGABAA_RE_TC;
    double* __restrict  _ptr_array_IGABAA_RE_TC_propoCondMult = _array_IGABAA_RE_TC_propoCondMult;
    double* __restrict  _ptr_array_IGABAA_RE_TC_sGABAA = _array_IGABAA_RE_TC_sGABAA;
    double* __restrict  _ptr_array_TC_group_v = _array_TC_group_v;


    //// MAIN CODE ////////////
    const int _target_size = N_post;

    // Set all the target variable values to zero
    #pragma omp parallel for schedule(static)
    for (int _target_idx=0; _target_idx<_target_size; _target_idx++)
    {
        _ptr_array_TC_group_iGABAA_RE_TC[_target_idx + 0] = 0;
    }

    // scalar code
    const size_t _vectorisation_idx = -1;
        


    for(int _idx=0; _idx<_ptr_array_IGABAA_RE_TC_N[0]; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _postsynaptic_idx = _ptr_array_IGABAA_RE_TC__synaptic_post[_idx];
        const double EGABAA = _ptr_array_IGABAA_RE_TC_EGABAA[_idx];
        const double Npre = _ptr_array_IGABAA_RE_TC_Npre[_idx];
        const double gGABAA = _ptr_array_IGABAA_RE_TC_gGABAA[_idx];
        const double propoCondMult = _ptr_array_IGABAA_RE_TC_propoCondMult[_idx];
        const double sGABAA = _ptr_array_IGABAA_RE_TC_sGABAA[_idx];
        const double v_post = _ptr_array_TC_group_v[_postsynaptic_idx];
        const double _synaptic_var = 1.0f*((((- propoCondMult) * gGABAA) * sGABAA) * (v_post - EGABAA))/Npre;

        _ptr_array_TC_group_iGABAA_RE_TC[_ptr_array_IGABAA_RE_TC__synaptic_post[_idx]] += _synaptic_var;
    }

}


