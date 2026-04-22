#include "code_objects/INMDA_PYso_PYdr_up_codeobject.h"
#include "objects.h"
#include "brianlib/common_math.h"
#include "brianlib/stdint_compat.h"
#include<chrono>
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<climits>
#include "brianlib/stdint_compat.h"
#include "synapses_classes.h"

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



void _run_INMDA_PYso_PYdr_up_codeobject()
{
    using namespace brian;


    ///// CONSTANTS ///////////
    int32_t* const _array_INMDA_PYso_PYdr__synaptic_pre = _dynamic_array_INMDA_PYso_PYdr__synaptic_pre.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr__synaptic_pre[0];
const size_t _num_synaptic_pre = _dynamic_array_INMDA_PYso_PYdr__synaptic_pre.size();
double* const _array_INMDA_PYso_PYdr_spike_activity = _dynamic_array_INMDA_PYso_PYdr_spike_activity.empty()? 0 : &_dynamic_array_INMDA_PYso_PYdr_spike_activity[0];
const size_t _numspike_activity = _dynamic_array_INMDA_PYso_PYdr_spike_activity.size();
    ///// POINTERS ////////////
        
    int32_t* __restrict  _ptr_array_INMDA_PYso_PYdr__synaptic_pre = _array_INMDA_PYso_PYdr__synaptic_pre;
    double* __restrict  _ptr_array_INMDA_PYso_PYdr_spike_activity = _array_INMDA_PYso_PYdr_spike_activity;



    // This is only needed for the _debugmsg function below

    // scalar code
    const size_t _vectorisation_idx = -1;
        


    #pragma omp parallel
    {
    std::vector<int> *_spiking_synapses = INMDA_PYso_PYdr_up.peek();
    const int _num_spiking_synapses = _spiking_synapses->size();

    #pragma omp for schedule(static)
    for(int _spiking_synapse_idx=0;
        _spiking_synapse_idx<_num_spiking_synapses;
        _spiking_synapse_idx++)
    {
        const size_t _idx = (*_spiking_synapses)[_spiking_synapse_idx];
        const size_t _vectorisation_idx = _idx;
                
        double spike_activity;
        spike_activity = 1;
        _ptr_array_INMDA_PYso_PYdr_spike_activity[_idx] = spike_activity;

    }

    }

}

void _debugmsg_INMDA_PYso_PYdr_up_codeobject()
{
    using namespace brian;
    std::cout << "Number of synapses: " << _dynamic_array_INMDA_PYso_PYdr__synaptic_pre.size() << endl;
}

