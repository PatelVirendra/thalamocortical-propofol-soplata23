#include "code_objects/IAMPA_PYso_PYdr_summed_variable_iAMPA_PYso_PYdr_post_codeobject.h"
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
        
    template <typename T>
    static inline T _clip(const T value, const double a_min, const double a_max)
    {
        if (value < a_min)
            return a_min;
        if (value > a_max)
            return a_max;
        return value;
    }
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



void _run_IAMPA_PYso_PYdr_summed_variable_iAMPA_PYso_PYdr_post_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    double* const _array_IAMPA_PYso_PYdr_EAMPA = _dynamic_array_IAMPA_PYso_PYdr_EAMPA.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_EAMPA[0];
const size_t _numEAMPA = _dynamic_array_IAMPA_PYso_PYdr_EAMPA.size();
const size_t _numN = 1;
const int64_t N_post = 100;
double* const _array_IAMPA_PYso_PYdr_Npost = _dynamic_array_IAMPA_PYso_PYdr_Npost.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_Npost[0];
const size_t _numNpost = _dynamic_array_IAMPA_PYso_PYdr_Npost.size();
double* const _array_IAMPA_PYso_PYdr_Npre = _dynamic_array_IAMPA_PYso_PYdr_Npre.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_Npre[0];
const size_t _numNpre = _dynamic_array_IAMPA_PYso_PYdr_Npre.size();
int32_t* const _array_IAMPA_PYso_PYdr__synaptic_post = _dynamic_array_IAMPA_PYso_PYdr__synaptic_post.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr__synaptic_post[0];
const size_t _num_synaptic_post = _dynamic_array_IAMPA_PYso_PYdr__synaptic_post.size();
double* const _array_IAMPA_PYso_PYdr_gAMPA = _dynamic_array_IAMPA_PYso_PYdr_gAMPA.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_gAMPA[0];
const size_t _numgAMPA = _dynamic_array_IAMPA_PYso_PYdr_gAMPA.size();
const size_t _numiAMPA_PYso_PYdr_post = 100;
double* const _array_IAMPA_PYso_PYdr_radius = _dynamic_array_IAMPA_PYso_PYdr_radius.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_radius[0];
const size_t _numradius = _dynamic_array_IAMPA_PYso_PYdr_radius.size();
char* const _array_IAMPA_PYso_PYdr_remove_recurrent_bool = _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool[0];
const size_t _numremove_recurrent_bool = _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.size();
double* const _array_IAMPA_PYso_PYdr_res_AMPA = _dynamic_array_IAMPA_PYso_PYdr_res_AMPA.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_res_AMPA[0];
const size_t _numres_AMPA = _dynamic_array_IAMPA_PYso_PYdr_res_AMPA.size();
double* const _array_IAMPA_PYso_PYdr_sAMPA = _dynamic_array_IAMPA_PYso_PYdr_sAMPA.empty()? 0 : &_dynamic_array_IAMPA_PYso_PYdr_sAMPA[0];
const size_t _numsAMPA = _dynamic_array_IAMPA_PYso_PYdr_sAMPA.size();
const size_t _numv_post = 100;
const size_t _num_postsynaptic_idx = _dynamic_array_IAMPA_PYso_PYdr__synaptic_post.size();
    ///// POINTERS ////////////
        
    double* __restrict  _ptr_array_IAMPA_PYso_PYdr_EAMPA = _array_IAMPA_PYso_PYdr_EAMPA;
    int32_t*   _ptr_array_IAMPA_PYso_PYdr_N = _array_IAMPA_PYso_PYdr_N;
    double* __restrict  _ptr_array_IAMPA_PYso_PYdr_Npost = _array_IAMPA_PYso_PYdr_Npost;
    double* __restrict  _ptr_array_IAMPA_PYso_PYdr_Npre = _array_IAMPA_PYso_PYdr_Npre;
    int32_t* __restrict  _ptr_array_IAMPA_PYso_PYdr__synaptic_post = _array_IAMPA_PYso_PYdr__synaptic_post;
    double* __restrict  _ptr_array_IAMPA_PYso_PYdr_gAMPA = _array_IAMPA_PYso_PYdr_gAMPA;
    double* __restrict  _ptr_array_PYdr_group_iAMPA_PYso_PYdr = _array_PYdr_group_iAMPA_PYso_PYdr;
    double* __restrict  _ptr_array_IAMPA_PYso_PYdr_radius = _array_IAMPA_PYso_PYdr_radius;
    char* __restrict  _ptr_array_IAMPA_PYso_PYdr_remove_recurrent_bool = _array_IAMPA_PYso_PYdr_remove_recurrent_bool;
    double* __restrict  _ptr_array_IAMPA_PYso_PYdr_res_AMPA = _array_IAMPA_PYso_PYdr_res_AMPA;
    double* __restrict  _ptr_array_IAMPA_PYso_PYdr_sAMPA = _array_IAMPA_PYso_PYdr_sAMPA;
    double* __restrict  _ptr_array_PYdr_group_v = _array_PYdr_group_v;


    //// MAIN CODE ////////////
    const int _target_size = N_post;

    // Set all the target variable values to zero
    #pragma omp parallel for schedule(static)
    for (int _target_idx=0; _target_idx<_target_size; _target_idx++)
    {
        _ptr_array_PYdr_group_iAMPA_PYso_PYdr[_target_idx + 0] = 0;
    }

    // scalar code
    const size_t _vectorisation_idx = -1;
        


    for(int _idx=0; _idx<_ptr_array_IAMPA_PYso_PYdr_N[0]; _idx++)
    {
        // vector code
        const size_t _vectorisation_idx = _idx;
                
        const int32_t _postsynaptic_idx = _ptr_array_IAMPA_PYso_PYdr__synaptic_post[_idx];
        const double EAMPA = _ptr_array_IAMPA_PYso_PYdr_EAMPA[_idx];
        const double Npost = _ptr_array_IAMPA_PYso_PYdr_Npost[_idx];
        const double Npre = _ptr_array_IAMPA_PYso_PYdr_Npre[_idx];
        const double gAMPA = _ptr_array_IAMPA_PYso_PYdr_gAMPA[_idx];
        const double radius = _ptr_array_IAMPA_PYso_PYdr_radius[_idx];
        const char remove_recurrent_bool = _ptr_array_IAMPA_PYso_PYdr_remove_recurrent_bool[_idx];
        const double res_AMPA = _ptr_array_IAMPA_PYso_PYdr_res_AMPA[_idx];
        const double sAMPA = _ptr_array_IAMPA_PYso_PYdr_sAMPA[_idx];
        const double v_post = _ptr_array_PYdr_group_v[_postsynaptic_idx];
        double normalizing_factor;
        if(!remove_recurrent_bool)
            normalizing_factor = _clip(1.0f*((1.0 + (2.0 * radius)) * Npre)/Npost, 0, Npre);
        else 
            normalizing_factor = _clip(1.0f*((2.0 * radius) * Npre)/Npost, 0, Npre);
        const double _synaptic_var = (((- (1.0f*gAMPA/normalizing_factor)) * res_AMPA) * sAMPA) * (v_post - EAMPA);

        _ptr_array_PYdr_group_iAMPA_PYso_PYdr[_ptr_array_IAMPA_PYso_PYdr__synaptic_post[_idx]] += _synaptic_var;
    }

}


