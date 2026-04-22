#include <stdlib.h>
#include "objects.h"
#include <csignal>
#include <ctime>
#include <time.h>
#include <omp.h>
#include "run.h"
#include "brianlib/common_math.h"

#include "code_objects/IAMPA_PYso_IN_down_codeobject.h"
#include "code_objects/IAMPA_PYso_IN_down_push_spikes.h"
#include "code_objects/before_run_IAMPA_PYso_IN_down_push_spikes.h"
#include "code_objects/IAMPA_PYso_IN_stateupdater_codeobject.h"
#include "code_objects/IAMPA_PYso_IN_summed_variable_iAMPA_PYso_IN_post_codeobject.h"
#include "code_objects/IAMPA_PYso_IN_synapses_create_array_codeobject.h"
#include "code_objects/IAMPA_PYso_IN_up_codeobject.h"
#include "code_objects/IAMPA_PYso_IN_up_push_spikes.h"
#include "code_objects/before_run_IAMPA_PYso_IN_up_push_spikes.h"
#include "code_objects/IAMPA_PYso_PYdr_down_codeobject.h"
#include "code_objects/IAMPA_PYso_PYdr_down_push_spikes.h"
#include "code_objects/before_run_IAMPA_PYso_PYdr_down_push_spikes.h"
#include "code_objects/IAMPA_PYso_PYdr_stateupdater_codeobject.h"
#include "code_objects/IAMPA_PYso_PYdr_summed_variable_iAMPA_PYso_PYdr_post_codeobject.h"
#include "code_objects/IAMPA_PYso_PYdr_synapses_create_array_codeobject.h"
#include "code_objects/IAMPA_PYso_PYdr_up_codeobject.h"
#include "code_objects/IAMPA_PYso_PYdr_up_push_spikes.h"
#include "code_objects/before_run_IAMPA_PYso_PYdr_up_push_spikes.h"
#include "code_objects/IAMPA_PYso_RE_stateupdater_codeobject.h"
#include "code_objects/IAMPA_PYso_RE_summed_variable_iAMPA_PYso_RE_post_codeobject.h"
#include "code_objects/IAMPA_PYso_RE_synapses_create_array_codeobject.h"
#include "code_objects/IAMPA_PYso_TC_stateupdater_codeobject.h"
#include "code_objects/IAMPA_PYso_TC_summed_variable_iAMPA_PYso_TC_post_codeobject.h"
#include "code_objects/IAMPA_PYso_TC_synapses_create_array_codeobject.h"
#include "code_objects/IAMPA_TC_IN_stateupdater_codeobject.h"
#include "code_objects/IAMPA_TC_IN_summed_variable_iAMPA_TC_IN_post_codeobject.h"
#include "code_objects/IAMPA_TC_IN_synapses_create_array_codeobject.h"
#include "code_objects/IAMPA_TC_PYdr_stateupdater_codeobject.h"
#include "code_objects/IAMPA_TC_PYdr_summed_variable_iAMPA_TC_PYdr_post_codeobject.h"
#include "code_objects/IAMPA_TC_PYdr_synapses_create_array_codeobject.h"
#include "code_objects/IAMPA_TC_RE_stateupdater_codeobject.h"
#include "code_objects/IAMPA_TC_RE_summed_variable_iAMPA_TC_RE_post_codeobject.h"
#include "code_objects/IAMPA_TC_RE_synapses_create_generator_codeobject.h"
#include "code_objects/ICOM_PYdr_PYso_summed_variable_iCOM_post_codeobject.h"
#include "code_objects/ICOM_PYdr_PYso_synapses_create_generator_codeobject.h"
#include "code_objects/ICOM_PYso_PYdr_summed_variable_iCOM_post_codeobject.h"
#include "code_objects/ICOM_PYso_PYdr_synapses_create_generator_codeobject.h"
#include "code_objects/IGABAA_IN_IN_down_codeobject.h"
#include "code_objects/IGABAA_IN_IN_down_push_spikes.h"
#include "code_objects/before_run_IGABAA_IN_IN_down_push_spikes.h"
#include "code_objects/IGABAA_IN_IN_stateupdater_codeobject.h"
#include "code_objects/IGABAA_IN_IN_summed_variable_iGABAA_IN_IN_post_codeobject.h"
#include "code_objects/IGABAA_IN_IN_synapses_create_array_codeobject.h"
#include "code_objects/IGABAA_IN_IN_up_codeobject.h"
#include "code_objects/IGABAA_IN_IN_up_push_spikes.h"
#include "code_objects/before_run_IGABAA_IN_IN_up_push_spikes.h"
#include "code_objects/IGABAA_IN_PYso_down_codeobject.h"
#include "code_objects/IGABAA_IN_PYso_down_push_spikes.h"
#include "code_objects/before_run_IGABAA_IN_PYso_down_push_spikes.h"
#include "code_objects/IGABAA_IN_PYso_stateupdater_codeobject.h"
#include "code_objects/IGABAA_IN_PYso_summed_variable_iGABAA_IN_PYso_post_codeobject.h"
#include "code_objects/IGABAA_IN_PYso_synapses_create_array_codeobject.h"
#include "code_objects/IGABAA_IN_PYso_up_codeobject.h"
#include "code_objects/IGABAA_IN_PYso_up_push_spikes.h"
#include "code_objects/before_run_IGABAA_IN_PYso_up_push_spikes.h"
#include "code_objects/IGABAA_RE_RE_stateupdater_codeobject.h"
#include "code_objects/IGABAA_RE_RE_summed_variable_iGABAA_RE_RE_post_codeobject.h"
#include "code_objects/IGABAA_RE_RE_synapses_create_generator_codeobject.h"
#include "code_objects/IGABAA_RE_TC_stateupdater_codeobject.h"
#include "code_objects/IGABAA_RE_TC_summed_variable_iGABAA_RE_TC_post_codeobject.h"
#include "code_objects/IGABAA_RE_TC_synapses_create_generator_codeobject.h"
#include "code_objects/IGABAB_RE_TC_stateupdater_codeobject.h"
#include "code_objects/IGABAB_RE_TC_summed_variable_iGABAB_RE_TC_post_codeobject.h"
#include "code_objects/IGABAB_RE_TC_synapses_create_generator_codeobject.h"
#include "code_objects/IKNa_PYso_PYdr_stateupdater_codeobject.h"
#include "code_objects/IKNa_PYso_PYdr_summed_variable_iKNa_post_codeobject.h"
#include "code_objects/IKNa_PYso_PYdr_synapses_create_generator_codeobject.h"
#include "code_objects/IN_group_spike_resetter_codeobject.h"
#include "code_objects/IN_group_spike_thresholder_codeobject.h"
#include "code_objects/after_run_IN_group_spike_thresholder_codeobject.h"
#include "code_objects/IN_group_stateupdater_codeobject.h"
#include "code_objects/IN_spikemon_codeobject.h"
#include "code_objects/INMDA_PYso_IN_down_codeobject.h"
#include "code_objects/INMDA_PYso_IN_down_push_spikes.h"
#include "code_objects/before_run_INMDA_PYso_IN_down_push_spikes.h"
#include "code_objects/INMDA_PYso_IN_stateupdater_codeobject.h"
#include "code_objects/INMDA_PYso_IN_summed_variable_iNMDA_PYso_IN_post_codeobject.h"
#include "code_objects/INMDA_PYso_IN_synapses_create_array_codeobject.h"
#include "code_objects/INMDA_PYso_IN_up_codeobject.h"
#include "code_objects/INMDA_PYso_IN_up_push_spikes.h"
#include "code_objects/before_run_INMDA_PYso_IN_up_push_spikes.h"
#include "code_objects/INMDA_PYso_PYdr_down_codeobject.h"
#include "code_objects/INMDA_PYso_PYdr_down_push_spikes.h"
#include "code_objects/before_run_INMDA_PYso_PYdr_down_push_spikes.h"
#include "code_objects/INMDA_PYso_PYdr_stateupdater_codeobject.h"
#include "code_objects/INMDA_PYso_PYdr_summed_variable_iNMDA_PYso_PYdr_post_codeobject.h"
#include "code_objects/INMDA_PYso_PYdr_synapses_create_array_codeobject.h"
#include "code_objects/INMDA_PYso_PYdr_up_codeobject.h"
#include "code_objects/INMDA_PYso_PYdr_up_push_spikes.h"
#include "code_objects/before_run_INMDA_PYso_PYdr_up_push_spikes.h"
#include "code_objects/poissongroup_1_spike_thresholder_codeobject.h"
#include "code_objects/after_run_poissongroup_1_spike_thresholder_codeobject.h"
#include "code_objects/poissongroup_2_spike_thresholder_codeobject.h"
#include "code_objects/after_run_poissongroup_2_spike_thresholder_codeobject.h"
#include "code_objects/poissongroup_spike_thresholder_codeobject.h"
#include "code_objects/after_run_poissongroup_spike_thresholder_codeobject.h"
#include "code_objects/PYdr_group_spike_resetter_codeobject.h"
#include "code_objects/PYdr_group_spike_thresholder_codeobject.h"
#include "code_objects/after_run_PYdr_group_spike_thresholder_codeobject.h"
#include "code_objects/PYdr_group_stateupdater_codeobject.h"
#include "code_objects/PYdr_spikemon_codeobject.h"
#include "code_objects/PYso_group_spike_resetter_codeobject.h"
#include "code_objects/PYso_group_spike_thresholder_codeobject.h"
#include "code_objects/after_run_PYso_group_spike_thresholder_codeobject.h"
#include "code_objects/PYso_group_stateupdater_codeobject.h"
#include "code_objects/PYso_spikemon_codeobject.h"
#include "code_objects/RE_group_spike_resetter_codeobject.h"
#include "code_objects/RE_group_spike_thresholder_codeobject.h"
#include "code_objects/after_run_RE_group_spike_thresholder_codeobject.h"
#include "code_objects/RE_group_stateupdater_codeobject.h"
#include "code_objects/RE_spikemon_codeobject.h"
#include "code_objects/synapses_1_pre_codeobject.h"
#include "code_objects/synapses_1_pre_push_spikes.h"
#include "code_objects/before_run_synapses_1_pre_push_spikes.h"
#include "code_objects/synapses_1_synapses_create_generator_codeobject.h"
#include "code_objects/synapses_2_pre_codeobject.h"
#include "code_objects/synapses_2_pre_push_spikes.h"
#include "code_objects/before_run_synapses_2_pre_push_spikes.h"
#include "code_objects/synapses_2_synapses_create_generator_codeobject.h"
#include "code_objects/synapses_pre_codeobject.h"
#include "code_objects/synapses_pre_push_spikes.h"
#include "code_objects/before_run_synapses_pre_push_spikes.h"
#include "code_objects/synapses_synapses_create_generator_codeobject.h"
#include "code_objects/TC_group_spike_resetter_codeobject.h"
#include "code_objects/TC_group_spike_thresholder_codeobject.h"
#include "code_objects/after_run_TC_group_spike_thresholder_codeobject.h"
#include "code_objects/TC_group_stateupdater_codeobject.h"
#include "code_objects/TC_group_subexpression_update_codeobject.h"
#include "code_objects/TC_spikemon_codeobject.h"


#include <iostream>
#include <fstream>
#include <string>




void set_from_command_line(const std::vector<std::string> args)
{
    for (const auto& arg : args) {
		// Split into two parts
		size_t equal_sign = arg.find("=");
		auto name = arg.substr(0, equal_sign);
		auto value = arg.substr(equal_sign + 1, arg.length());
		brian::set_variable_by_name(name, value);
	}
}

void _int_handler(int signal_num) {
	if (Network::_globally_running && !Network::_globally_stopped) {
		Network::_globally_stopped = true;
	} else {
		std::signal(signal_num, SIG_DFL);
		std::raise(signal_num);
	}
}

int main(int argc, char **argv)
{
	std::signal(SIGINT, _int_handler);
	std::random_device _rd;
	std::vector<std::string> args(argv + 1, argv + argc);
	if (args.size() >=2 && args[0] == "--results_dir")
	{
		brian::results_dir = args[1];
		#ifdef DEBUG
		std::cout << "Setting results dir to '" << brian::results_dir << "'" << std::endl;
		#endif
		args.erase(args.begin(), args.begin()+2);
	}
        

	brian_start();
        

	{
		using namespace brian;

		omp_set_dynamic(0);
omp_set_num_threads(2);
                
        _array_defaultclock_timestep[0] = 0;
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 0.0001;
        _array_defaultclock_dt[0] = 1e-05;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYdr_group_lastspike; i++)
                        {
                            _array_PYdr_group_lastspike[i] = - 10000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYdr_group_not_refractory; i++)
                        {
                            _array_PYdr_group_not_refractory[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYdr_group_CaBuffer; i++)
                        {
                            _array_PYdr_group_CaBuffer[i] = _static_array__array_PYdr_group_CaBuffer[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYdr_group_v; i++)
                        {
                            _array_PYdr_group_v[i] = _static_array__array_PYdr_group_v[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYdr_group_gPoisson_PYdr; i++)
                        {
                            _array_PYdr_group_gPoisson_PYdr[i] = 0.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_poissongroup_rates; i++)
                        {
                            _array_poissongroup_rates[i] = 40.0;
                        }
                        
        _run_synapses_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYso_group_lastspike; i++)
                        {
                            _array_PYso_group_lastspike[i] = - 10000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYso_group_not_refractory; i++)
                        {
                            _array_PYso_group_not_refractory[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYso_group_hNa; i++)
                        {
                            _array_PYso_group_hNa[i] = _static_array__array_PYso_group_hNa[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYso_group_nK; i++)
                        {
                            _array_PYso_group_nK[i] = _static_array__array_PYso_group_nK[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYso_group_hA; i++)
                        {
                            _array_PYso_group_hA[i] = _static_array__array_PYso_group_hA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYso_group_mKS; i++)
                        {
                            _array_PYso_group_mKS[i] = _static_array__array_PYso_group_mKS[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_PYso_group_v; i++)
                        {
                            _array_PYso_group_v[i] = _static_array__array_PYso_group_v[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IN_group_lastspike; i++)
                        {
                            _array_IN_group_lastspike[i] = - 10000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IN_group_not_refractory; i++)
                        {
                            _array_IN_group_not_refractory[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IN_group_hNa_IN; i++)
                        {
                            _array_IN_group_hNa_IN[i] = _static_array__array_IN_group_hNa_IN[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IN_group_nK_IN; i++)
                        {
                            _array_IN_group_nK_IN[i] = _static_array__array_IN_group_nK_IN[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IN_group_v; i++)
                        {
                            _array_IN_group_v[i] = _static_array__array_IN_group_v[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_lastspike; i++)
                        {
                            _array_TC_group_lastspike[i] = - 10000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_not_refractory; i++)
                        {
                            _array_TC_group_not_refractory[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_Ca_TC; i++)
                        {
                            _array_TC_group_Ca_TC[i] = _static_array__array_TC_group_Ca_TC[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_mNa_TC; i++)
                        {
                            _array_TC_group_mNa_TC[i] = _static_array__array_TC_group_mNa_TC[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_hNa_TC; i++)
                        {
                            _array_TC_group_hNa_TC[i] = _static_array__array_TC_group_hNa_TC[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_nK_TC; i++)
                        {
                            _array_TC_group_nK_TC[i] = _static_array__array_TC_group_nK_TC[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_hT_TC; i++)
                        {
                            _array_TC_group_hT_TC[i] = _static_array__array_TC_group_hT_TC[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_Open; i++)
                        {
                            _array_TC_group_Open[i] = _static_array__array_TC_group_Open[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_Pone; i++)
                        {
                            _array_TC_group_Pone[i] = _static_array__array_TC_group_Pone[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_OpenLocked; i++)
                        {
                            _array_TC_group_OpenLocked[i] = _static_array__array_TC_group_OpenLocked[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_v; i++)
                        {
                            _array_TC_group_v[i] = _static_array__array_TC_group_v[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_TC_group_gPoisson_TC; i++)
                        {
                            _array_TC_group_gPoisson_TC[i] = 0.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_poissongroup_1_rates; i++)
                        {
                            _array_poissongroup_1_rates[i] = 40.0;
                        }
                        
        _run_synapses_1_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_lastspike; i++)
                        {
                            _array_RE_group_lastspike[i] = - 10000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_not_refractory; i++)
                        {
                            _array_RE_group_not_refractory[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_mNa_RE; i++)
                        {
                            _array_RE_group_mNa_RE[i] = _static_array__array_RE_group_mNa_RE[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_hNa_RE; i++)
                        {
                            _array_RE_group_hNa_RE[i] = _static_array__array_RE_group_hNa_RE[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_nK_RE; i++)
                        {
                            _array_RE_group_nK_RE[i] = _static_array__array_RE_group_nK_RE[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_mT_RE; i++)
                        {
                            _array_RE_group_mT_RE[i] = _static_array__array_RE_group_mT_RE[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_hT_RE; i++)
                        {
                            _array_RE_group_hT_RE[i] = _static_array__array_RE_group_hT_RE[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_v; i++)
                        {
                            _array_RE_group_v[i] = _static_array__array_RE_group_v[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_RE_group_gPoisson_RE; i++)
                        {
                            _array_RE_group_gPoisson_RE[i] = 0.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_poissongroup_2_rates; i++)
                        {
                            _array_poissongroup_2_rates[i] = 40.0;
                        }
                        
        _run_synapses_2_synapses_create_generator_codeobject();
        _run_ICOM_PYdr_PYso_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_ICOM_PYdr_PYso_gCOM.size(); i++)
                        {
                            _dynamic_array_ICOM_PYdr_PYso_gCOM[i] = _static_array__dynamic_array_ICOM_PYdr_PYso_gCOM[i];
                        }
                        
        _run_ICOM_PYso_PYdr_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_ICOM_PYso_PYdr_gCOM.size(); i++)
                        {
                            _dynamic_array_ICOM_PYso_PYdr_gCOM[i] = _static_array__dynamic_array_ICOM_PYso_PYdr_gCOM[i];
                        }
                        
        _run_IKNa_PYso_PYdr_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IKNa_PYso_PYdr_hNa_local.size(); i++)
                        {
                            _dynamic_array_IKNa_PYso_PYdr_hNa_local[i] = _static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IKNa_PYso_PYdr_concNa.size(); i++)
                        {
                            _dynamic_array_IKNa_PYso_PYdr_concNa[i] = _static_array__dynamic_array_IKNa_PYso_PYdr_concNa[i];
                        }
                        
        _dynamic_array_IAMPA_PYso_PYdr_delay.resize(1);
        _dynamic_array_IAMPA_PYso_PYdr_delay.resize(1);
        _dynamic_array_IAMPA_PYso_PYdr_delay[0] = 1e-05;
        _dynamic_array_IAMPA_PYso_PYdr_delay_1.resize(1);
        _dynamic_array_IAMPA_PYso_PYdr_delay_1.resize(1);
        _dynamic_array_IAMPA_PYso_PYdr_delay_1[0] = 3.0000000000000004e-05;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_PYdr_sources; i++)
                        {
                            _array_IAMPA_PYso_PYdr_sources[i] = _static_array__array_IAMPA_PYso_PYdr_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_PYdr_targets; i++)
                        {
                            _array_IAMPA_PYso_PYdr_targets[i] = _static_array__array_IAMPA_PYso_PYdr_targets[i];
                        }
                        
        _run_IAMPA_PYso_PYdr_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_Npre.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_Npre[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_Npost.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_Npost[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_gAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_gAMPA[i] = 0.05;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_EAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_EAMPA[i] = 0.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA[i] = 3480.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_tauAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_tauAMPA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_deprFactor.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_deprFactor[i] = 0.9;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA[i] = 0.4;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_radius.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_sAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_sAMPA[i] = _static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_res_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_res_AMPA[i] = _static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_PYdr_spike_activity.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_PYdr_spike_activity[i] = 0;
                        }
                        
        _dynamic_array_INMDA_PYso_PYdr_delay.resize(1);
        _dynamic_array_INMDA_PYso_PYdr_delay.resize(1);
        _dynamic_array_INMDA_PYso_PYdr_delay[0] = 1e-05;
        _dynamic_array_INMDA_PYso_PYdr_delay_1.resize(1);
        _dynamic_array_INMDA_PYso_PYdr_delay_1.resize(1);
        _dynamic_array_INMDA_PYso_PYdr_delay_1[0] = 3.0000000000000004e-05;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_INMDA_PYso_PYdr_sources; i++)
                        {
                            _array_INMDA_PYso_PYdr_sources[i] = _static_array__array_INMDA_PYso_PYdr_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_INMDA_PYso_PYdr_targets; i++)
                        {
                            _array_INMDA_PYso_PYdr_targets[i] = _static_array__array_INMDA_PYso_PYdr_targets[i];
                        }
                        
        _run_INMDA_PYso_PYdr_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_Npre.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_Npre[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_Npost.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_Npost[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_gNMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_gNMDA[i] = 0.025699999999999997;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_ENMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_ENMDA[i] = 0.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA[i] = 500.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA[i] = 0.1;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA[i] = 3480.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_deprFactor.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_deprFactor[i] = 0.9;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA[i] = 0.4;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_radius.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_sNMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_sNMDA[i] = _static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_xNMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_xNMDA[i] = _static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_res_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_res_NMDA[i] = _static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_PYdr_spike_activity.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_PYdr_spike_activity[i] = 0;
                        }
                        
        _dynamic_array_IAMPA_PYso_IN_delay.resize(1);
        _dynamic_array_IAMPA_PYso_IN_delay.resize(1);
        _dynamic_array_IAMPA_PYso_IN_delay[0] = 1e-05;
        _dynamic_array_IAMPA_PYso_IN_delay_1.resize(1);
        _dynamic_array_IAMPA_PYso_IN_delay_1.resize(1);
        _dynamic_array_IAMPA_PYso_IN_delay_1[0] = 3.0000000000000004e-05;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_IN_sources; i++)
                        {
                            _array_IAMPA_PYso_IN_sources[i] = _static_array__array_IAMPA_PYso_IN_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_IN_targets; i++)
                        {
                            _array_IAMPA_PYso_IN_targets[i] = _static_array__array_IAMPA_PYso_IN_targets[i];
                        }
                        
        _run_IAMPA_PYso_IN_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_Npre.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_Npre[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_Npost.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_gAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_gAMPA[i] = 10.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_EAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_EAMPA[i] = 0.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_alpha_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_alpha_AMPA[i] = 3480.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_tauAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_tauAMPA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_deprFactor.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_deprFactor[i] = 0.9;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_tauRes_AMPA[i] = 0.4;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_radius.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool[i] = false;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_sAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_sAMPA[i] = _static_array__dynamic_array_IAMPA_PYso_IN_sAMPA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_res_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_res_AMPA[i] = _static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_IN_spike_activity.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_IN_spike_activity[i] = 0;
                        }
                        
        _dynamic_array_INMDA_PYso_IN_delay.resize(1);
        _dynamic_array_INMDA_PYso_IN_delay.resize(1);
        _dynamic_array_INMDA_PYso_IN_delay[0] = 1e-05;
        _dynamic_array_INMDA_PYso_IN_delay_1.resize(1);
        _dynamic_array_INMDA_PYso_IN_delay_1.resize(1);
        _dynamic_array_INMDA_PYso_IN_delay_1[0] = 3.0000000000000004e-05;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_INMDA_PYso_IN_sources; i++)
                        {
                            _array_INMDA_PYso_IN_sources[i] = _static_array__array_INMDA_PYso_IN_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_INMDA_PYso_IN_targets; i++)
                        {
                            _array_INMDA_PYso_IN_targets[i] = _static_array__array_INMDA_PYso_IN_targets[i];
                        }
                        
        _run_INMDA_PYso_IN_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_Npre.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_Npre[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_Npost.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_gNMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_gNMDA[i] = 0.025;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_ENMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_ENMDA[i] = 0.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_alphaS_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_alphaS_NMDA[i] = 500.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_tauS_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_tauS_NMDA[i] = 0.1;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_alphaX_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_alphaX_NMDA[i] = 3480.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_tauX_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_tauX_NMDA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_deprFactor.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_deprFactor[i] = 0.9;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_tauRes_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_tauRes_NMDA[i] = 0.4;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_radius.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_remove_recurrent_bool[i] = false;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_sNMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_sNMDA[i] = _static_array__dynamic_array_INMDA_PYso_IN_sNMDA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_xNMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_xNMDA[i] = _static_array__dynamic_array_INMDA_PYso_IN_xNMDA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_res_NMDA.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_res_NMDA[i] = _static_array__dynamic_array_INMDA_PYso_IN_res_NMDA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_INMDA_PYso_IN_spike_activity.size(); i++)
                        {
                            _dynamic_array_INMDA_PYso_IN_spike_activity[i] = 0;
                        }
                        
        _dynamic_array_IGABAA_IN_PYso_delay.resize(1);
        _dynamic_array_IGABAA_IN_PYso_delay.resize(1);
        _dynamic_array_IGABAA_IN_PYso_delay[0] = 1e-05;
        _dynamic_array_IGABAA_IN_PYso_delay_1.resize(1);
        _dynamic_array_IGABAA_IN_PYso_delay_1.resize(1);
        _dynamic_array_IGABAA_IN_PYso_delay_1[0] = 3.0000000000000004e-05;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IGABAA_IN_PYso_sources; i++)
                        {
                            _array_IGABAA_IN_PYso_sources[i] = _static_array__array_IGABAA_IN_PYso_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IGABAA_IN_PYso_targets; i++)
                        {
                            _array_IGABAA_IN_PYso_targets[i] = _static_array__array_IGABAA_IN_PYso_targets[i];
                        }
                        
        _run_IGABAA_IN_PYso_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_Npre.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_Npost.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_Npost[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_gGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_gGABAA[i] = 1.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_EGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_EGABAA[i] = - 0.07;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_alpha_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_alpha_GABAA[i] = 1000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_tauGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_tauGABAA[i] = 0.005;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_deprFactor.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_deprFactor[i] = 0.9;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA[i] = 0.4;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_radius.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool[i] = false;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_propoCondMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_propoCondMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_propoTauMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_propoTauMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_sGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_sGABAA[i] = _static_array__dynamic_array_IGABAA_IN_PYso_sGABAA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_res_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_res_GABAA[i] = _static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_PYso_spike_activity.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_PYso_spike_activity[i] = 0;
                        }
                        
        _dynamic_array_IGABAA_IN_IN_delay.resize(1);
        _dynamic_array_IGABAA_IN_IN_delay.resize(1);
        _dynamic_array_IGABAA_IN_IN_delay[0] = 1e-05;
        _dynamic_array_IGABAA_IN_IN_delay_1.resize(1);
        _dynamic_array_IGABAA_IN_IN_delay_1.resize(1);
        _dynamic_array_IGABAA_IN_IN_delay_1[0] = 3.0000000000000004e-05;
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IGABAA_IN_IN_sources; i++)
                        {
                            _array_IGABAA_IN_IN_sources[i] = _static_array__array_IGABAA_IN_IN_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IGABAA_IN_IN_targets; i++)
                        {
                            _array_IGABAA_IN_IN_targets[i] = _static_array__array_IGABAA_IN_IN_targets[i];
                        }
                        
        _run_IGABAA_IN_IN_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_Npre.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_Npost.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_gGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_gGABAA[i] = 0.00825;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_EGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_EGABAA[i] = - 0.07;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_alpha_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_alpha_GABAA[i] = 1000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_tauGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_tauGABAA[i] = 0.005;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_deprFactor.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_deprFactor[i] = 0.9;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_tauRes_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_tauRes_GABAA[i] = 0.4;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_radius.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_remove_recurrent_bool[i] = true;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_propoCondMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_propoCondMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_propoTauMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_propoTauMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_sGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_sGABAA[i] = _static_array__dynamic_array_IGABAA_IN_IN_sGABAA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_res_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_res_GABAA[i] = _static_array__dynamic_array_IGABAA_IN_IN_res_GABAA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_IN_IN_spike_activity.size(); i++)
                        {
                            _dynamic_array_IGABAA_IN_IN_spike_activity[i] = 0;
                        }
                        
        _run_IAMPA_TC_RE_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_RE_Npre.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_RE_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_RE_Npost.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_RE_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_RE_gAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_RE_gAMPA[i] = 4.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_RE_EAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_RE_EAMPA[i] = 0.001;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_RE_tauAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_RE_tauAMPA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_RE_P_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_RE_P_AMPA[i] = 5000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_RE_sAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_RE_sAMPA[i] = _static_array__dynamic_array_IAMPA_TC_RE_sAMPA[i];
                        }
                        
        _run_IGABAA_RE_RE_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_Npre.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_Npost.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_gGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_gGABAA[i] = 1.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_EGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_EGABAA[i] = - 0.08;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_tauGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_tauGABAA[i] = 0.005;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_P_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_P_GABAA[i] = 2000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_propoCondMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_propoCondMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_propoTauMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_propoTauMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_RE_sGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_RE_sGABAA[i] = _static_array__dynamic_array_IGABAA_RE_RE_sGABAA[i];
                        }
                        
        _run_IGABAA_RE_TC_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_Npre.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_Npost.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_gGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_gGABAA[i] = 1.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_EGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_EGABAA[i] = - 0.08;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_tauGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_tauGABAA[i] = 0.005;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_P_GABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_P_GABAA[i] = 2000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_propoCondMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_propoCondMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_propoTauMult.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_propoTauMult[i] = 3;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAA_RE_TC_sGABAA.size(); i++)
                        {
                            _dynamic_array_IGABAA_RE_TC_sGABAA[i] = _static_array__dynamic_array_IGABAA_RE_TC_sGABAA[i];
                        }
                        
        _run_IGABAB_RE_TC_synapses_create_generator_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_Npre.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_Npost.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_gGABAB.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_gGABAB[i] = 0.009999999999999998;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_EGABAB.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_EGABAB[i] = - 0.095;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_K1.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_K1[i] = 500.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_K2.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_K2[i] = 1.2;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_K3.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_K3[i] = 180.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_K4.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_K4[i] = 34.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_rGABAB.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_rGABAB[i] = _static_array__dynamic_array_IGABAB_RE_TC_rGABAB[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IGABAB_RE_TC_sGABAB.size(); i++)
                        {
                            _dynamic_array_IGABAB_RE_TC_sGABAB[i] = _static_array__dynamic_array_IGABAB_RE_TC_sGABAB[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_TC_PYdr_sources; i++)
                        {
                            _array_IAMPA_TC_PYdr_sources[i] = _static_array__array_IAMPA_TC_PYdr_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_TC_PYdr_targets; i++)
                        {
                            _array_IAMPA_TC_PYdr_targets[i] = _static_array__array_IAMPA_TC_PYdr_targets[i];
                        }
                        
        _run_IAMPA_TC_PYdr_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_Npre.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_Npost.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_Npost[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_gAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_gAMPA[i] = 0.05;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_EAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_EAMPA[i] = 0.001;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_tauAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_tauAMPA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_P_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_P_AMPA[i] = 5000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_radius.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool[i] = false;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_PYdr_sAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_PYdr_sAMPA[i] = _static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_TC_IN_sources; i++)
                        {
                            _array_IAMPA_TC_IN_sources[i] = _static_array__array_IAMPA_TC_IN_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_TC_IN_targets; i++)
                        {
                            _array_IAMPA_TC_IN_targets[i] = _static_array__array_IAMPA_TC_IN_targets[i];
                        }
                        
        _run_IAMPA_TC_IN_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_Npre.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_Npre[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_Npost.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_gAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_gAMPA[i] = 1.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_EAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_EAMPA[i] = 0.001;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_tauAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_tauAMPA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_P_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_P_AMPA[i] = 5000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_radius.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_remove_recurrent_bool[i] = false;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_TC_IN_sAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_TC_IN_sAMPA[i] = _static_array__dynamic_array_IAMPA_TC_IN_sAMPA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_TC_sources; i++)
                        {
                            _array_IAMPA_PYso_TC_sources[i] = _static_array__array_IAMPA_PYso_TC_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_TC_targets; i++)
                        {
                            _array_IAMPA_PYso_TC_targets[i] = _static_array__array_IAMPA_PYso_TC_targets[i];
                        }
                        
        _run_IAMPA_PYso_TC_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_Npre.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_Npre[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_Npost.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_gAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_gAMPA[i] = 4.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_EAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_EAMPA[i] = 0.001;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_tauAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_tauAMPA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_P_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_P_AMPA[i] = 5000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_radius.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool[i] = false;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_TC_sAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_TC_sAMPA[i] = _static_array__dynamic_array_IAMPA_PYso_TC_sAMPA[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_RE_sources; i++)
                        {
                            _array_IAMPA_PYso_RE_sources[i] = _static_array__array_IAMPA_PYso_RE_sources[i];
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_num__array_IAMPA_PYso_RE_targets; i++)
                        {
                            _array_IAMPA_PYso_RE_targets[i] = _static_array__array_IAMPA_PYso_RE_targets[i];
                        }
                        
        _run_IAMPA_PYso_RE_synapses_create_array_codeobject();
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_Npre.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_Npre[i] = 100;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_Npost.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_Npost[i] = 20;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE[i] = 2.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE[i] = 0.001;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_tauAMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_tauAMPA[i] = 0.002;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_P_AMPA.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_P_AMPA[i] = 5000.0;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_radius.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_radius[i] = 10;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool[i] = false;
                        }
                        
        
                        #pragma omp for schedule(static)
                        for(int i=0; i<_dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.size(); i++)
                        {
                            _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE[i] = _static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE[i];
                        }
                        
        _array_clock_timestep[0] = 0;
        _array_clock_dt[0] = 0.0001;
        _array_clock_dt[0] = 0.0001;
        for (int _i=0; _i<2; _i++)
            brian::_random_generators[_i].seed(0L + _i);
        _array_defaultclock_timestep[0] = 0;
        _array_defaultclock_t[0] = 0.0;
        _before_run_IAMPA_PYso_IN_down_push_spikes();
        _before_run_IAMPA_PYso_IN_up_push_spikes();
        _before_run_IAMPA_PYso_PYdr_down_push_spikes();
        _before_run_IAMPA_PYso_PYdr_up_push_spikes();
        _before_run_IGABAA_IN_IN_down_push_spikes();
        _before_run_IGABAA_IN_IN_up_push_spikes();
        _before_run_IGABAA_IN_PYso_down_push_spikes();
        _before_run_IGABAA_IN_PYso_up_push_spikes();
        _before_run_INMDA_PYso_IN_down_push_spikes();
        _before_run_INMDA_PYso_IN_up_push_spikes();
        _before_run_INMDA_PYso_PYdr_down_push_spikes();
        _before_run_INMDA_PYso_PYdr_up_push_spikes();
        _before_run_synapses_1_pre_push_spikes();
        _before_run_synapses_2_pre_push_spikes();
        _before_run_synapses_pre_push_spikes();
        network.clear();
        network.add(&defaultclock, _run_TC_group_subexpression_update_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_IN_summed_variable_iAMPA_PYso_IN_post_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_PYdr_summed_variable_iAMPA_PYso_PYdr_post_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_RE_summed_variable_iAMPA_PYso_RE_post_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_TC_summed_variable_iAMPA_PYso_TC_post_codeobject);
        network.add(&defaultclock, _run_IAMPA_TC_IN_summed_variable_iAMPA_TC_IN_post_codeobject);
        network.add(&defaultclock, _run_IAMPA_TC_PYdr_summed_variable_iAMPA_TC_PYdr_post_codeobject);
        network.add(&defaultclock, _run_IAMPA_TC_RE_summed_variable_iAMPA_TC_RE_post_codeobject);
        network.add(&defaultclock, _run_ICOM_PYdr_PYso_summed_variable_iCOM_post_codeobject);
        network.add(&defaultclock, _run_ICOM_PYso_PYdr_summed_variable_iCOM_post_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_IN_summed_variable_iGABAA_IN_IN_post_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_PYso_summed_variable_iGABAA_IN_PYso_post_codeobject);
        network.add(&defaultclock, _run_IGABAA_RE_RE_summed_variable_iGABAA_RE_RE_post_codeobject);
        network.add(&defaultclock, _run_IGABAA_RE_TC_summed_variable_iGABAA_RE_TC_post_codeobject);
        network.add(&defaultclock, _run_IGABAB_RE_TC_summed_variable_iGABAB_RE_TC_post_codeobject);
        network.add(&defaultclock, _run_IKNa_PYso_PYdr_summed_variable_iKNa_post_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_IN_summed_variable_iNMDA_PYso_IN_post_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_PYdr_summed_variable_iNMDA_PYso_PYdr_post_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_IN_stateupdater_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_PYdr_stateupdater_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_RE_stateupdater_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_TC_stateupdater_codeobject);
        network.add(&defaultclock, _run_IAMPA_TC_IN_stateupdater_codeobject);
        network.add(&defaultclock, _run_IAMPA_TC_PYdr_stateupdater_codeobject);
        network.add(&defaultclock, _run_IAMPA_TC_RE_stateupdater_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_IN_stateupdater_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_PYso_stateupdater_codeobject);
        network.add(&defaultclock, _run_IGABAA_RE_RE_stateupdater_codeobject);
        network.add(&defaultclock, _run_IGABAA_RE_TC_stateupdater_codeobject);
        network.add(&defaultclock, _run_IGABAB_RE_TC_stateupdater_codeobject);
        network.add(&defaultclock, _run_IKNa_PYso_PYdr_stateupdater_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_IN_stateupdater_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_PYdr_stateupdater_codeobject);
        network.add(&defaultclock, _run_IN_group_stateupdater_codeobject);
        network.add(&defaultclock, _run_PYdr_group_stateupdater_codeobject);
        network.add(&defaultclock, _run_PYso_group_stateupdater_codeobject);
        network.add(&defaultclock, _run_RE_group_stateupdater_codeobject);
        network.add(&defaultclock, _run_TC_group_stateupdater_codeobject);
        network.add(&defaultclock, _run_IN_group_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_PYdr_group_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_PYso_group_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_RE_group_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_TC_group_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_poissongroup_1_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_poissongroup_2_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_poissongroup_spike_thresholder_codeobject);
        network.add(&defaultclock, _run_IN_spikemon_codeobject);
        network.add(&defaultclock, _run_PYdr_spikemon_codeobject);
        network.add(&defaultclock, _run_PYso_spikemon_codeobject);
        network.add(&defaultclock, _run_RE_spikemon_codeobject);
        network.add(&defaultclock, _run_TC_spikemon_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_IN_down_push_spikes);
        network.add(&defaultclock, _run_IAMPA_PYso_IN_down_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_IN_up_push_spikes);
        network.add(&defaultclock, _run_IAMPA_PYso_IN_up_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_PYdr_down_push_spikes);
        network.add(&defaultclock, _run_IAMPA_PYso_PYdr_down_codeobject);
        network.add(&defaultclock, _run_IAMPA_PYso_PYdr_up_push_spikes);
        network.add(&defaultclock, _run_IAMPA_PYso_PYdr_up_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_IN_down_push_spikes);
        network.add(&defaultclock, _run_IGABAA_IN_IN_down_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_IN_up_push_spikes);
        network.add(&defaultclock, _run_IGABAA_IN_IN_up_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_PYso_down_push_spikes);
        network.add(&defaultclock, _run_IGABAA_IN_PYso_down_codeobject);
        network.add(&defaultclock, _run_IGABAA_IN_PYso_up_push_spikes);
        network.add(&defaultclock, _run_IGABAA_IN_PYso_up_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_IN_down_push_spikes);
        network.add(&defaultclock, _run_INMDA_PYso_IN_down_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_IN_up_push_spikes);
        network.add(&defaultclock, _run_INMDA_PYso_IN_up_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_PYdr_down_push_spikes);
        network.add(&defaultclock, _run_INMDA_PYso_PYdr_down_codeobject);
        network.add(&defaultclock, _run_INMDA_PYso_PYdr_up_push_spikes);
        network.add(&defaultclock, _run_INMDA_PYso_PYdr_up_codeobject);
        network.add(&defaultclock, _run_synapses_1_pre_push_spikes);
        network.add(&defaultclock, _run_synapses_1_pre_codeobject);
        network.add(&defaultclock, _run_synapses_2_pre_push_spikes);
        network.add(&defaultclock, _run_synapses_2_pre_codeobject);
        network.add(&defaultclock, _run_synapses_pre_push_spikes);
        network.add(&defaultclock, _run_synapses_pre_codeobject);
        network.add(&defaultclock, _run_IN_group_spike_resetter_codeobject);
        network.add(&defaultclock, _run_PYdr_group_spike_resetter_codeobject);
        network.add(&defaultclock, _run_PYso_group_spike_resetter_codeobject);
        network.add(&defaultclock, _run_RE_group_spike_resetter_codeobject);
        network.add(&defaultclock, _run_TC_group_spike_resetter_codeobject);
        set_from_command_line(args);
        network.run(30.0, NULL, 10.0);
        _after_run_IN_group_spike_thresholder_codeobject();
        _after_run_PYdr_group_spike_thresholder_codeobject();
        _after_run_PYso_group_spike_thresholder_codeobject();
        _after_run_RE_group_spike_thresholder_codeobject();
        _after_run_TC_group_spike_thresholder_codeobject();
        _after_run_poissongroup_1_spike_thresholder_codeobject();
        _after_run_poissongroup_2_spike_thresholder_codeobject();
        _after_run_poissongroup_spike_thresholder_codeobject();
        #ifdef DEBUG
        _debugmsg_IN_spikemon_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_PYdr_spikemon_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_PYso_spikemon_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_RE_spikemon_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_TC_spikemon_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IAMPA_PYso_IN_down_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IAMPA_PYso_IN_up_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IAMPA_PYso_PYdr_down_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IAMPA_PYso_PYdr_up_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IGABAA_IN_IN_down_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IGABAA_IN_IN_up_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IGABAA_IN_PYso_down_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_IGABAA_IN_PYso_up_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_INMDA_PYso_IN_down_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_INMDA_PYso_IN_up_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_INMDA_PYso_PYdr_down_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_INMDA_PYso_PYdr_up_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_synapses_1_pre_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_synapses_2_pre_codeobject();
        #endif
        
        #ifdef DEBUG
        _debugmsg_synapses_pre_codeobject();
        #endif

	}
        

	brian_end();
        

	return 0;
}