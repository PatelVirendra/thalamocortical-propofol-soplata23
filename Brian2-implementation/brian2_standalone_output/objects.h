
#ifndef _BRIAN_OBJECTS_H
#define _BRIAN_OBJECTS_H

#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include<chrono>
#include<random>
#include<vector>
#include <omp.h>

namespace brian {

extern std::string results_dir;

class RandomGenerator {
    private:
        std::mt19937 gen;
        double stored_gauss;
        bool has_stored_gauss = false;
    public:
        RandomGenerator() {
            seed();
        }
        void seed() {
            std::random_device rd;
            gen.seed(rd());
            has_stored_gauss = false;
        }
        void seed(unsigned long seed) {
            gen.seed(seed);
            has_stored_gauss = false;
        }
        // Allow exporting/setting the internal state of the random generator
        friend std::ostream& operator<<(std::ostream& out, const RandomGenerator& rng);
        friend std::istream& operator>>(std::istream& in, RandomGenerator& rng);

        double rand() {
            /* shifts : 67108864 = 0x4000000, 9007199254740992 = 0x20000000000000 */
            const long a = gen() >> 5;
            const long b = gen() >> 6;
            return (a * 67108864.0 + b) / 9007199254740992.0;
        }

        double randn() {
            if (has_stored_gauss) {
                const double tmp = stored_gauss;
                has_stored_gauss = false;
                return tmp;
            }
            else {
                double f, x1, x2, r2;

                do {
                    x1 = 2.0*rand() - 1.0;
                    x2 = 2.0*rand() - 1.0;
                    r2 = x1*x1 + x2*x2;
                }
                while (r2 >= 1.0 || r2 == 0.0);

                /* Box-Muller transform */
                f = sqrt(-2.0*log(r2)/r2);
                /* Keep for next call */
                stored_gauss = f*x1;
                has_stored_gauss = true;
                return f*x2;
            }
        }
};

extern std::ostream& operator<<(std::ostream& out, const RandomGenerator& rng);
extern std::istream& operator>>(std::istream& in, RandomGenerator& rng);

// In OpenMP we need one state per thread
extern std::vector< RandomGenerator > _random_generators;

//////////////// clocks ///////////////////
extern Clock defaultclock;

//////////////// networks /////////////////
extern Network network;



void set_variable_by_name(std::string, std::string);

//////////////// dynamic arrays ///////////
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN__synaptic_pre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_alpha_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_delay;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_delay_1;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_deprFactor;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_EAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_gAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN_N_incoming;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN_N_outgoing;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_Npost;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_Npre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_radius;
extern std::vector<char> _dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_res_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_sAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_spike_activity;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_tauAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_IN_tauRes_AMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr__synaptic_pre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_delay;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_delay_1;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_deprFactor;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_EAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_gAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr_N_incoming;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr_N_outgoing;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_Npost;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_Npre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_radius;
extern std::vector<char> _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_sAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_spike_activity;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_tauAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE__synaptic_pre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE_N_incoming;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE_N_outgoing;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_Npost;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_Npre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_P_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_radius;
extern std::vector<char> _dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
extern std::vector<double> _dynamic_array_IAMPA_PYso_RE_tauAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC__synaptic_pre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_EAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_gAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC_N_incoming;
extern std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC_N_outgoing;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_Npost;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_Npre;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_P_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_radius;
extern std::vector<char> _dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_sAMPA;
extern std::vector<double> _dynamic_array_IAMPA_PYso_TC_tauAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_IN__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_IN__synaptic_pre;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_EAMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_gAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_IN_N_incoming;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_IN_N_outgoing;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_Npost;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_Npre;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_P_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_radius;
extern std::vector<char> _dynamic_array_IAMPA_TC_IN_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_sAMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_IN_tauAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr__synaptic_pre;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_EAMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_gAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr_N_incoming;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr_N_outgoing;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_Npost;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_Npre;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_P_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_radius;
extern std::vector<char> _dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_sAMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_PYdr_tauAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_RE__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_RE__synaptic_pre;
extern std::vector<double> _dynamic_array_IAMPA_TC_RE_EAMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_RE_gAMPA;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_RE_N_incoming;
extern std::vector<int32_t> _dynamic_array_IAMPA_TC_RE_N_outgoing;
extern std::vector<double> _dynamic_array_IAMPA_TC_RE_Npost;
extern std::vector<double> _dynamic_array_IAMPA_TC_RE_Npre;
extern std::vector<double> _dynamic_array_IAMPA_TC_RE_P_AMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_RE_sAMPA;
extern std::vector<double> _dynamic_array_IAMPA_TC_RE_tauAMPA;
extern std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso__synaptic_post;
extern std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso__synaptic_pre;
extern std::vector<double> _dynamic_array_ICOM_PYdr_PYso_gCOM;
extern std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso_N_incoming;
extern std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso_N_outgoing;
extern std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr__synaptic_post;
extern std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr__synaptic_pre;
extern std::vector<double> _dynamic_array_ICOM_PYso_PYdr_gCOM;
extern std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr_N_incoming;
extern std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr_N_outgoing;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_IN__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_IN__synaptic_pre;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_alpha_GABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_delay;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_delay_1;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_deprFactor;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_EGABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_gGABAA;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_IN_N_incoming;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_IN_N_outgoing;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_Npost;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_Npre;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_propoCondMult;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_propoTauMult;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_radius;
extern std::vector<char> _dynamic_array_IGABAA_IN_IN_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_res_GABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_sGABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_spike_activity;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_tauGABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_IN_tauRes_GABAA;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso__synaptic_pre;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_alpha_GABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_delay;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_delay_1;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_deprFactor;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_EGABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_gGABAA;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso_N_incoming;
extern std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso_N_outgoing;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_Npost;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_Npre;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_propoCondMult;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_propoTauMult;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_radius;
extern std::vector<char> _dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_res_GABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_sGABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_spike_activity;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_tauGABAA;
extern std::vector<double> _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_RE__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_RE__synaptic_pre;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_EGABAA;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_gGABAA;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_RE_N_incoming;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_RE_N_outgoing;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_Npost;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_Npre;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_P_GABAA;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_propoCondMult;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_propoTauMult;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_sGABAA;
extern std::vector<double> _dynamic_array_IGABAA_RE_RE_tauGABAA;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_TC__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_TC__synaptic_pre;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_EGABAA;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_gGABAA;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_TC_N_incoming;
extern std::vector<int32_t> _dynamic_array_IGABAA_RE_TC_N_outgoing;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_Npost;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_Npre;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_P_GABAA;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_propoCondMult;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_propoTauMult;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_sGABAA;
extern std::vector<double> _dynamic_array_IGABAA_RE_TC_tauGABAA;
extern std::vector<int32_t> _dynamic_array_IGABAB_RE_TC__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IGABAB_RE_TC__synaptic_pre;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_EGABAB;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_gGABAB;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_K1;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_K2;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_K3;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_K4;
extern std::vector<int32_t> _dynamic_array_IGABAB_RE_TC_N_incoming;
extern std::vector<int32_t> _dynamic_array_IGABAB_RE_TC_N_outgoing;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_Npost;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_Npre;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_rGABAB;
extern std::vector<double> _dynamic_array_IGABAB_RE_TC_sGABAB;
extern std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr__synaptic_post;
extern std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr__synaptic_pre;
extern std::vector<double> _dynamic_array_IKNa_PYso_PYdr_concNa;
extern std::vector<double> _dynamic_array_IKNa_PYso_PYdr_hNa_local;
extern std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr_N_incoming;
extern std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr_N_outgoing;
extern std::vector<int32_t> _dynamic_array_IN_spikemon_i;
extern std::vector<double> _dynamic_array_IN_spikemon_t;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_IN__synaptic_post;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_IN__synaptic_pre;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_alphaS_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_alphaX_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_delay;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_delay_1;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_deprFactor;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_ENMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_gNMDA;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_IN_N_incoming;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_IN_N_outgoing;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_Npost;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_Npre;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_radius;
extern std::vector<char> _dynamic_array_INMDA_PYso_IN_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_res_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_sNMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_spike_activity;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_tauRes_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_tauS_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_tauX_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_IN_xNMDA;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr__synaptic_post;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr__synaptic_pre;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_delay;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_delay_1;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_deprFactor;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_ENMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_gNMDA;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr_N_incoming;
extern std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr_N_outgoing;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_Npost;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_Npre;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_radius;
extern std::vector<char> _dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_res_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_sNMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_spike_activity;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA;
extern std::vector<double> _dynamic_array_INMDA_PYso_PYdr_xNMDA;
extern std::vector<int32_t> _dynamic_array_PYdr_spikemon_i;
extern std::vector<double> _dynamic_array_PYdr_spikemon_t;
extern std::vector<int32_t> _dynamic_array_PYso_spikemon_i;
extern std::vector<double> _dynamic_array_PYso_spikemon_t;
extern std::vector<int32_t> _dynamic_array_RE_spikemon_i;
extern std::vector<double> _dynamic_array_RE_spikemon_t;
extern std::vector<int32_t> _dynamic_array_synapses_1__synaptic_post;
extern std::vector<int32_t> _dynamic_array_synapses_1__synaptic_pre;
extern std::vector<double> _dynamic_array_synapses_1_delay;
extern std::vector<int32_t> _dynamic_array_synapses_1_N_incoming;
extern std::vector<int32_t> _dynamic_array_synapses_1_N_outgoing;
extern std::vector<int32_t> _dynamic_array_synapses_2__synaptic_post;
extern std::vector<int32_t> _dynamic_array_synapses_2__synaptic_pre;
extern std::vector<double> _dynamic_array_synapses_2_delay;
extern std::vector<int32_t> _dynamic_array_synapses_2_N_incoming;
extern std::vector<int32_t> _dynamic_array_synapses_2_N_outgoing;
extern std::vector<int32_t> _dynamic_array_synapses__synaptic_post;
extern std::vector<int32_t> _dynamic_array_synapses__synaptic_pre;
extern std::vector<double> _dynamic_array_synapses_delay;
extern std::vector<int32_t> _dynamic_array_synapses_N_incoming;
extern std::vector<int32_t> _dynamic_array_synapses_N_outgoing;
extern std::vector<int32_t> _dynamic_array_TC_spikemon_i;
extern std::vector<double> _dynamic_array_TC_spikemon_t;

//////////////// arrays ///////////////////
extern double *_array_clock_dt;
extern const int _num__array_clock_dt;
extern double *_array_clock_t;
extern const int _num__array_clock_t;
extern int64_t *_array_clock_timestep;
extern const int _num__array_clock_timestep;
extern double *_array_defaultclock_dt;
extern const int _num__array_defaultclock_dt;
extern double *_array_defaultclock_t;
extern const int _num__array_defaultclock_t;
extern int64_t *_array_defaultclock_timestep;
extern const int _num__array_defaultclock_timestep;
extern int32_t *_array_IAMPA_PYso_IN_N;
extern const int _num__array_IAMPA_PYso_IN_N;
extern int32_t *_array_IAMPA_PYso_IN_sources;
extern const int _num__array_IAMPA_PYso_IN_sources;
extern int32_t *_array_IAMPA_PYso_IN_targets;
extern const int _num__array_IAMPA_PYso_IN_targets;
extern int32_t *_array_IAMPA_PYso_PYdr_N;
extern const int _num__array_IAMPA_PYso_PYdr_N;
extern int32_t *_array_IAMPA_PYso_PYdr_sources;
extern const int _num__array_IAMPA_PYso_PYdr_sources;
extern int32_t *_array_IAMPA_PYso_PYdr_targets;
extern const int _num__array_IAMPA_PYso_PYdr_targets;
extern int32_t *_array_IAMPA_PYso_RE_N;
extern const int _num__array_IAMPA_PYso_RE_N;
extern int32_t *_array_IAMPA_PYso_RE_sources;
extern const int _num__array_IAMPA_PYso_RE_sources;
extern int32_t *_array_IAMPA_PYso_RE_targets;
extern const int _num__array_IAMPA_PYso_RE_targets;
extern int32_t *_array_IAMPA_PYso_TC_N;
extern const int _num__array_IAMPA_PYso_TC_N;
extern int32_t *_array_IAMPA_PYso_TC_sources;
extern const int _num__array_IAMPA_PYso_TC_sources;
extern int32_t *_array_IAMPA_PYso_TC_targets;
extern const int _num__array_IAMPA_PYso_TC_targets;
extern int32_t *_array_IAMPA_TC_IN_N;
extern const int _num__array_IAMPA_TC_IN_N;
extern int32_t *_array_IAMPA_TC_IN_sources;
extern const int _num__array_IAMPA_TC_IN_sources;
extern int32_t *_array_IAMPA_TC_IN_targets;
extern const int _num__array_IAMPA_TC_IN_targets;
extern int32_t *_array_IAMPA_TC_PYdr_N;
extern const int _num__array_IAMPA_TC_PYdr_N;
extern int32_t *_array_IAMPA_TC_PYdr_sources;
extern const int _num__array_IAMPA_TC_PYdr_sources;
extern int32_t *_array_IAMPA_TC_PYdr_targets;
extern const int _num__array_IAMPA_TC_PYdr_targets;
extern int32_t *_array_IAMPA_TC_RE_N;
extern const int _num__array_IAMPA_TC_RE_N;
extern int32_t *_array_ICOM_PYdr_PYso_N;
extern const int _num__array_ICOM_PYdr_PYso_N;
extern int32_t *_array_ICOM_PYso_PYdr_N;
extern const int _num__array_ICOM_PYso_PYdr_N;
extern int32_t *_array_IGABAA_IN_IN_N;
extern const int _num__array_IGABAA_IN_IN_N;
extern int32_t *_array_IGABAA_IN_IN_sources;
extern const int _num__array_IGABAA_IN_IN_sources;
extern int32_t *_array_IGABAA_IN_IN_targets;
extern const int _num__array_IGABAA_IN_IN_targets;
extern int32_t *_array_IGABAA_IN_PYso_N;
extern const int _num__array_IGABAA_IN_PYso_N;
extern int32_t *_array_IGABAA_IN_PYso_sources;
extern const int _num__array_IGABAA_IN_PYso_sources;
extern int32_t *_array_IGABAA_IN_PYso_targets;
extern const int _num__array_IGABAA_IN_PYso_targets;
extern int32_t *_array_IGABAA_RE_RE_N;
extern const int _num__array_IGABAA_RE_RE_N;
extern int32_t *_array_IGABAA_RE_TC_N;
extern const int _num__array_IGABAA_RE_TC_N;
extern int32_t *_array_IGABAB_RE_TC_N;
extern const int _num__array_IGABAB_RE_TC_N;
extern int32_t *_array_IKNa_PYso_PYdr_N;
extern const int _num__array_IKNa_PYso_PYdr_N;
extern int32_t *_array_IN_group__spikespace;
extern const int _num__array_IN_group__spikespace;
extern double *_array_IN_group_hNa_IN;
extern const int _num__array_IN_group_hNa_IN;
extern int32_t *_array_IN_group_i;
extern const int _num__array_IN_group_i;
extern double *_array_IN_group_iAMPA_PYso_IN;
extern const int _num__array_IN_group_iAMPA_PYso_IN;
extern double *_array_IN_group_iAMPA_TC_IN;
extern const int _num__array_IN_group_iAMPA_TC_IN;
extern double *_array_IN_group_iGABAA_IN_IN;
extern const int _num__array_IN_group_iGABAA_IN_IN;
extern double *_array_IN_group_iNMDA_PYso_IN;
extern const int _num__array_IN_group_iNMDA_PYso_IN;
extern double *_array_IN_group_lastspike;
extern const int _num__array_IN_group_lastspike;
extern double *_array_IN_group_nK_IN;
extern const int _num__array_IN_group_nK_IN;
extern char *_array_IN_group_not_refractory;
extern const int _num__array_IN_group_not_refractory;
extern double *_array_IN_group_v;
extern const int _num__array_IN_group_v;
extern int32_t *_array_IN_spikemon__source_idx;
extern const int _num__array_IN_spikemon__source_idx;
extern int32_t *_array_IN_spikemon_count;
extern const int _num__array_IN_spikemon_count;
extern int32_t *_array_IN_spikemon_N;
extern const int _num__array_IN_spikemon_N;
extern int32_t *_array_INMDA_PYso_IN_N;
extern const int _num__array_INMDA_PYso_IN_N;
extern int32_t *_array_INMDA_PYso_IN_sources;
extern const int _num__array_INMDA_PYso_IN_sources;
extern int32_t *_array_INMDA_PYso_IN_targets;
extern const int _num__array_INMDA_PYso_IN_targets;
extern int32_t *_array_INMDA_PYso_PYdr_N;
extern const int _num__array_INMDA_PYso_PYdr_N;
extern int32_t *_array_INMDA_PYso_PYdr_sources;
extern const int _num__array_INMDA_PYso_PYdr_sources;
extern int32_t *_array_INMDA_PYso_PYdr_targets;
extern const int _num__array_INMDA_PYso_PYdr_targets;
extern int32_t *_array_poissongroup_1__spikespace;
extern const int _num__array_poissongroup_1__spikespace;
extern int32_t *_array_poissongroup_1_i;
extern const int _num__array_poissongroup_1_i;
extern double *_array_poissongroup_1_rates;
extern const int _num__array_poissongroup_1_rates;
extern int32_t *_array_poissongroup_2__spikespace;
extern const int _num__array_poissongroup_2__spikespace;
extern int32_t *_array_poissongroup_2_i;
extern const int _num__array_poissongroup_2_i;
extern double *_array_poissongroup_2_rates;
extern const int _num__array_poissongroup_2_rates;
extern int32_t *_array_poissongroup__spikespace;
extern const int _num__array_poissongroup__spikespace;
extern int32_t *_array_poissongroup_i;
extern const int _num__array_poissongroup_i;
extern double *_array_poissongroup_rates;
extern const int _num__array_poissongroup_rates;
extern int32_t *_array_PYdr_group__spikespace;
extern const int _num__array_PYdr_group__spikespace;
extern double *_array_PYdr_group_CaBuffer;
extern const int _num__array_PYdr_group_CaBuffer;
extern double *_array_PYdr_group_gPoisson_PYdr;
extern const int _num__array_PYdr_group_gPoisson_PYdr;
extern int32_t *_array_PYdr_group_i;
extern const int _num__array_PYdr_group_i;
extern double *_array_PYdr_group_iAMPA_PYso_PYdr;
extern const int _num__array_PYdr_group_iAMPA_PYso_PYdr;
extern double *_array_PYdr_group_iAMPA_TC_PYdr;
extern const int _num__array_PYdr_group_iAMPA_TC_PYdr;
extern double *_array_PYdr_group_iCOM;
extern const int _num__array_PYdr_group_iCOM;
extern double *_array_PYdr_group_iNMDA_PYso_PYdr;
extern const int _num__array_PYdr_group_iNMDA_PYso_PYdr;
extern double *_array_PYdr_group_lastspike;
extern const int _num__array_PYdr_group_lastspike;
extern char *_array_PYdr_group_not_refractory;
extern const int _num__array_PYdr_group_not_refractory;
extern double *_array_PYdr_group_v;
extern const int _num__array_PYdr_group_v;
extern int32_t *_array_PYdr_spikemon__source_idx;
extern const int _num__array_PYdr_spikemon__source_idx;
extern int32_t *_array_PYdr_spikemon_count;
extern const int _num__array_PYdr_spikemon_count;
extern int32_t *_array_PYdr_spikemon_N;
extern const int _num__array_PYdr_spikemon_N;
extern int32_t *_array_PYso_group__spikespace;
extern const int _num__array_PYso_group__spikespace;
extern double *_array_PYso_group_hA;
extern const int _num__array_PYso_group_hA;
extern double *_array_PYso_group_hNa;
extern const int _num__array_PYso_group_hNa;
extern int32_t *_array_PYso_group_i;
extern const int _num__array_PYso_group_i;
extern double *_array_PYso_group_iCOM;
extern const int _num__array_PYso_group_iCOM;
extern double *_array_PYso_group_iGABAA_IN_PYso;
extern const int _num__array_PYso_group_iGABAA_IN_PYso;
extern double *_array_PYso_group_iKNa;
extern const int _num__array_PYso_group_iKNa;
extern double *_array_PYso_group_lastspike;
extern const int _num__array_PYso_group_lastspike;
extern double *_array_PYso_group_mKS;
extern const int _num__array_PYso_group_mKS;
extern double *_array_PYso_group_nK;
extern const int _num__array_PYso_group_nK;
extern char *_array_PYso_group_not_refractory;
extern const int _num__array_PYso_group_not_refractory;
extern double *_array_PYso_group_v;
extern const int _num__array_PYso_group_v;
extern int32_t *_array_PYso_spikemon__source_idx;
extern const int _num__array_PYso_spikemon__source_idx;
extern int32_t *_array_PYso_spikemon_count;
extern const int _num__array_PYso_spikemon_count;
extern int32_t *_array_PYso_spikemon_N;
extern const int _num__array_PYso_spikemon_N;
extern int32_t *_array_RE_group__spikespace;
extern const int _num__array_RE_group__spikespace;
extern double *_array_RE_group_gPoisson_RE;
extern const int _num__array_RE_group_gPoisson_RE;
extern double *_array_RE_group_hNa_RE;
extern const int _num__array_RE_group_hNa_RE;
extern double *_array_RE_group_hT_RE;
extern const int _num__array_RE_group_hT_RE;
extern int32_t *_array_RE_group_i;
extern const int _num__array_RE_group_i;
extern double *_array_RE_group_iAMPA_PYso_RE;
extern const int _num__array_RE_group_iAMPA_PYso_RE;
extern double *_array_RE_group_iAMPA_TC_RE;
extern const int _num__array_RE_group_iAMPA_TC_RE;
extern double *_array_RE_group_iGABAA_RE_RE;
extern const int _num__array_RE_group_iGABAA_RE_RE;
extern double *_array_RE_group_lastspike;
extern const int _num__array_RE_group_lastspike;
extern double *_array_RE_group_mNa_RE;
extern const int _num__array_RE_group_mNa_RE;
extern double *_array_RE_group_mT_RE;
extern const int _num__array_RE_group_mT_RE;
extern double *_array_RE_group_nK_RE;
extern const int _num__array_RE_group_nK_RE;
extern char *_array_RE_group_not_refractory;
extern const int _num__array_RE_group_not_refractory;
extern double *_array_RE_group_v;
extern const int _num__array_RE_group_v;
extern int32_t *_array_RE_spikemon__source_idx;
extern const int _num__array_RE_spikemon__source_idx;
extern int32_t *_array_RE_spikemon_count;
extern const int _num__array_RE_spikemon_count;
extern int32_t *_array_RE_spikemon_N;
extern const int _num__array_RE_spikemon_N;
extern int32_t *_array_synapses_1_N;
extern const int _num__array_synapses_1_N;
extern int32_t *_array_synapses_2_N;
extern const int _num__array_synapses_2_N;
extern int32_t *_array_synapses_N;
extern const int _num__array_synapses_N;
extern int32_t *_array_TC_group__spikespace;
extern const int _num__array_TC_group__spikespace;
extern double *_array_TC_group_Ca_TC;
extern const int _num__array_TC_group_Ca_TC;
extern double *_array_TC_group_dCa_extra;
extern const int _num__array_TC_group_dCa_extra;
extern double *_array_TC_group_gPoisson_TC;
extern const int _num__array_TC_group_gPoisson_TC;
extern double *_array_TC_group_hNa_TC;
extern const int _num__array_TC_group_hNa_TC;
extern double *_array_TC_group_hT_TC;
extern const int _num__array_TC_group_hT_TC;
extern int32_t *_array_TC_group_i;
extern const int _num__array_TC_group_i;
extern double *_array_TC_group_iAMPA_PYso_TC;
extern const int _num__array_TC_group_iAMPA_PYso_TC;
extern double *_array_TC_group_iGABAA_RE_TC;
extern const int _num__array_TC_group_iGABAA_RE_TC;
extern double *_array_TC_group_iGABAB_RE_TC;
extern const int _num__array_TC_group_iGABAB_RE_TC;
extern double *_array_TC_group_lastspike;
extern const int _num__array_TC_group_lastspike;
extern double *_array_TC_group_mNa_TC;
extern const int _num__array_TC_group_mNa_TC;
extern double *_array_TC_group_nK_TC;
extern const int _num__array_TC_group_nK_TC;
extern char *_array_TC_group_not_refractory;
extern const int _num__array_TC_group_not_refractory;
extern double *_array_TC_group_Open;
extern const int _num__array_TC_group_Open;
extern double *_array_TC_group_OpenLocked;
extern const int _num__array_TC_group_OpenLocked;
extern double *_array_TC_group_Pone;
extern const int _num__array_TC_group_Pone;
extern double *_array_TC_group_v;
extern const int _num__array_TC_group_v;
extern int32_t *_array_TC_spikemon__source_idx;
extern const int _num__array_TC_spikemon__source_idx;
extern int32_t *_array_TC_spikemon_count;
extern const int _num__array_TC_spikemon_count;
extern int32_t *_array_TC_spikemon_N;
extern const int _num__array_TC_spikemon_N;

//////////////// dynamic arrays 2d /////////

/////////////// static arrays /////////////
extern int32_t *_static_array__array_IAMPA_PYso_IN_sources;
extern const int _num__static_array__array_IAMPA_PYso_IN_sources;
extern int32_t *_static_array__array_IAMPA_PYso_IN_targets;
extern const int _num__static_array__array_IAMPA_PYso_IN_targets;
extern int32_t *_static_array__array_IAMPA_PYso_PYdr_sources;
extern const int _num__static_array__array_IAMPA_PYso_PYdr_sources;
extern int32_t *_static_array__array_IAMPA_PYso_PYdr_targets;
extern const int _num__static_array__array_IAMPA_PYso_PYdr_targets;
extern int32_t *_static_array__array_IAMPA_PYso_RE_sources;
extern const int _num__static_array__array_IAMPA_PYso_RE_sources;
extern int32_t *_static_array__array_IAMPA_PYso_RE_targets;
extern const int _num__static_array__array_IAMPA_PYso_RE_targets;
extern int32_t *_static_array__array_IAMPA_PYso_TC_sources;
extern const int _num__static_array__array_IAMPA_PYso_TC_sources;
extern int32_t *_static_array__array_IAMPA_PYso_TC_targets;
extern const int _num__static_array__array_IAMPA_PYso_TC_targets;
extern int32_t *_static_array__array_IAMPA_TC_IN_sources;
extern const int _num__static_array__array_IAMPA_TC_IN_sources;
extern int32_t *_static_array__array_IAMPA_TC_IN_targets;
extern const int _num__static_array__array_IAMPA_TC_IN_targets;
extern int32_t *_static_array__array_IAMPA_TC_PYdr_sources;
extern const int _num__static_array__array_IAMPA_TC_PYdr_sources;
extern int32_t *_static_array__array_IAMPA_TC_PYdr_targets;
extern const int _num__static_array__array_IAMPA_TC_PYdr_targets;
extern int32_t *_static_array__array_IGABAA_IN_IN_sources;
extern const int _num__static_array__array_IGABAA_IN_IN_sources;
extern int32_t *_static_array__array_IGABAA_IN_IN_targets;
extern const int _num__static_array__array_IGABAA_IN_IN_targets;
extern int32_t *_static_array__array_IGABAA_IN_PYso_sources;
extern const int _num__static_array__array_IGABAA_IN_PYso_sources;
extern int32_t *_static_array__array_IGABAA_IN_PYso_targets;
extern const int _num__static_array__array_IGABAA_IN_PYso_targets;
extern double *_static_array__array_IN_group_hNa_IN;
extern const int _num__static_array__array_IN_group_hNa_IN;
extern double *_static_array__array_IN_group_nK_IN;
extern const int _num__static_array__array_IN_group_nK_IN;
extern double *_static_array__array_IN_group_v;
extern const int _num__static_array__array_IN_group_v;
extern int32_t *_static_array__array_INMDA_PYso_IN_sources;
extern const int _num__static_array__array_INMDA_PYso_IN_sources;
extern int32_t *_static_array__array_INMDA_PYso_IN_targets;
extern const int _num__static_array__array_INMDA_PYso_IN_targets;
extern int32_t *_static_array__array_INMDA_PYso_PYdr_sources;
extern const int _num__static_array__array_INMDA_PYso_PYdr_sources;
extern int32_t *_static_array__array_INMDA_PYso_PYdr_targets;
extern const int _num__static_array__array_INMDA_PYso_PYdr_targets;
extern double *_static_array__array_PYdr_group_CaBuffer;
extern const int _num__static_array__array_PYdr_group_CaBuffer;
extern double *_static_array__array_PYdr_group_v;
extern const int _num__static_array__array_PYdr_group_v;
extern double *_static_array__array_PYso_group_hA;
extern const int _num__static_array__array_PYso_group_hA;
extern double *_static_array__array_PYso_group_hNa;
extern const int _num__static_array__array_PYso_group_hNa;
extern double *_static_array__array_PYso_group_mKS;
extern const int _num__static_array__array_PYso_group_mKS;
extern double *_static_array__array_PYso_group_nK;
extern const int _num__static_array__array_PYso_group_nK;
extern double *_static_array__array_PYso_group_v;
extern const int _num__static_array__array_PYso_group_v;
extern double *_static_array__array_RE_group_hNa_RE;
extern const int _num__static_array__array_RE_group_hNa_RE;
extern double *_static_array__array_RE_group_hT_RE;
extern const int _num__static_array__array_RE_group_hT_RE;
extern double *_static_array__array_RE_group_mNa_RE;
extern const int _num__static_array__array_RE_group_mNa_RE;
extern double *_static_array__array_RE_group_mT_RE;
extern const int _num__static_array__array_RE_group_mT_RE;
extern double *_static_array__array_RE_group_nK_RE;
extern const int _num__static_array__array_RE_group_nK_RE;
extern double *_static_array__array_RE_group_v;
extern const int _num__static_array__array_RE_group_v;
extern double *_static_array__array_TC_group_Ca_TC;
extern const int _num__static_array__array_TC_group_Ca_TC;
extern double *_static_array__array_TC_group_hNa_TC;
extern const int _num__static_array__array_TC_group_hNa_TC;
extern double *_static_array__array_TC_group_hT_TC;
extern const int _num__static_array__array_TC_group_hT_TC;
extern double *_static_array__array_TC_group_mNa_TC;
extern const int _num__static_array__array_TC_group_mNa_TC;
extern double *_static_array__array_TC_group_nK_TC;
extern const int _num__static_array__array_TC_group_nK_TC;
extern double *_static_array__array_TC_group_Open;
extern const int _num__static_array__array_TC_group_Open;
extern double *_static_array__array_TC_group_OpenLocked;
extern const int _num__static_array__array_TC_group_OpenLocked;
extern double *_static_array__array_TC_group_Pone;
extern const int _num__static_array__array_TC_group_Pone;
extern double *_static_array__array_TC_group_v;
extern const int _num__static_array__array_TC_group_v;
extern double *_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA;
extern const int _num__static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA;
extern double *_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA;
extern const int _num__static_array__dynamic_array_IAMPA_PYso_IN_sAMPA;
extern double *_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
extern const int _num__static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
extern double *_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA;
extern const int _num__static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA;
extern double *_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
extern const int _num__static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
extern double *_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA;
extern const int _num__static_array__dynamic_array_IAMPA_PYso_TC_sAMPA;
extern double *_static_array__dynamic_array_IAMPA_TC_IN_sAMPA;
extern const int _num__static_array__dynamic_array_IAMPA_TC_IN_sAMPA;
extern double *_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA;
extern const int _num__static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA;
extern double *_static_array__dynamic_array_IAMPA_TC_RE_sAMPA;
extern const int _num__static_array__dynamic_array_IAMPA_TC_RE_sAMPA;
extern double *_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM;
extern const int _num__static_array__dynamic_array_ICOM_PYdr_PYso_gCOM;
extern double *_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM;
extern const int _num__static_array__dynamic_array_ICOM_PYso_PYdr_gCOM;
extern double *_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA;
extern const int _num__static_array__dynamic_array_IGABAA_IN_IN_res_GABAA;
extern double *_static_array__dynamic_array_IGABAA_IN_IN_sGABAA;
extern const int _num__static_array__dynamic_array_IGABAA_IN_IN_sGABAA;
extern double *_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA;
extern const int _num__static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA;
extern double *_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA;
extern const int _num__static_array__dynamic_array_IGABAA_IN_PYso_sGABAA;
extern double *_static_array__dynamic_array_IGABAA_RE_RE_sGABAA;
extern const int _num__static_array__dynamic_array_IGABAA_RE_RE_sGABAA;
extern double *_static_array__dynamic_array_IGABAA_RE_TC_sGABAA;
extern const int _num__static_array__dynamic_array_IGABAA_RE_TC_sGABAA;
extern double *_static_array__dynamic_array_IGABAB_RE_TC_rGABAB;
extern const int _num__static_array__dynamic_array_IGABAB_RE_TC_rGABAB;
extern double *_static_array__dynamic_array_IGABAB_RE_TC_sGABAB;
extern const int _num__static_array__dynamic_array_IGABAB_RE_TC_sGABAB;
extern double *_static_array__dynamic_array_IKNa_PYso_PYdr_concNa;
extern const int _num__static_array__dynamic_array_IKNa_PYso_PYdr_concNa;
extern double *_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local;
extern const int _num__static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local;
extern double *_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA;
extern const int _num__static_array__dynamic_array_INMDA_PYso_IN_res_NMDA;
extern double *_static_array__dynamic_array_INMDA_PYso_IN_sNMDA;
extern const int _num__static_array__dynamic_array_INMDA_PYso_IN_sNMDA;
extern double *_static_array__dynamic_array_INMDA_PYso_IN_xNMDA;
extern const int _num__static_array__dynamic_array_INMDA_PYso_IN_xNMDA;
extern double *_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA;
extern const int _num__static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA;
extern double *_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA;
extern const int _num__static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA;
extern double *_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA;
extern const int _num__static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA;

//////////////// synapses /////////////////
// IAMPA_PYso_IN
extern SynapticPathway IAMPA_PYso_IN_down;
extern SynapticPathway IAMPA_PYso_IN_up;
// IAMPA_PYso_PYdr
extern SynapticPathway IAMPA_PYso_PYdr_down;
extern SynapticPathway IAMPA_PYso_PYdr_up;
// IAMPA_PYso_RE
// IAMPA_PYso_TC
// IAMPA_TC_IN
// IAMPA_TC_PYdr
// IAMPA_TC_RE
// ICOM_PYdr_PYso
// ICOM_PYso_PYdr
// IGABAA_IN_IN
extern SynapticPathway IGABAA_IN_IN_down;
extern SynapticPathway IGABAA_IN_IN_up;
// IGABAA_IN_PYso
extern SynapticPathway IGABAA_IN_PYso_down;
extern SynapticPathway IGABAA_IN_PYso_up;
// IGABAA_RE_RE
// IGABAA_RE_TC
// IGABAB_RE_TC
// IKNa_PYso_PYdr
// INMDA_PYso_IN
extern SynapticPathway INMDA_PYso_IN_down;
extern SynapticPathway INMDA_PYso_IN_up;
// INMDA_PYso_PYdr
extern SynapticPathway INMDA_PYso_PYdr_down;
extern SynapticPathway INMDA_PYso_PYdr_up;
// synapses
extern SynapticPathway synapses_pre;
// synapses_1
extern SynapticPathway synapses_1_pre;
// synapses_2
extern SynapticPathway synapses_2_pre;

// Profiling information for each code object
}

void _init_arrays();
void _load_arrays();
void _write_arrays();
void _dealloc_arrays();

#endif


