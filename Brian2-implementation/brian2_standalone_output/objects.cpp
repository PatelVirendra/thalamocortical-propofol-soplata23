

#include "objects.h"
#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include<chrono>
#include<random>
#include<vector>
#include<iostream>
#include<fstream>
#include<map>
#include<tuple>
#include<cstdlib>
#include<string>

namespace brian {

std::string results_dir = "results/";  // can be overwritten by --results_dir command line arg

// For multhreading, we need one generator for each thread.
std::vector< RandomGenerator > _random_generators;

std::ostream& operator<<(std::ostream& out, const RandomGenerator& rng)
{
    return out << rng.gen;
}

std::istream& operator>>(std::istream& in, RandomGenerator& rng)
{
    return in >> rng.gen;
}

//////////////// networks /////////////////
Network network;

void set_variable_from_value(std::string varname, char* var_pointer, size_t size, char value) {
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' to " << (value == 1 ? "True" : "False") << std::endl;
    #endif
    std::fill(var_pointer, var_pointer+size, value);
}

template<class T> void set_variable_from_value(std::string varname, T* var_pointer, size_t size, T value) {
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' to " << value << std::endl;
    #endif
    std::fill(var_pointer, var_pointer+size, value);
}

template<class T> void set_variable_from_file(std::string varname, T* var_pointer, size_t data_size, std::string filename) {
    ifstream f;
    streampos size;
    #ifdef DEBUG
    std::cout << "Setting '" << varname << "' from file '" << filename << "'" << std::endl;
    #endif
    f.open(filename, ios::in | ios::binary | ios::ate);
    size = f.tellg();
    if (size != data_size) {
        std::cerr << "Error reading '" << filename << "': file size " << size << " does not match expected size " << data_size << std::endl;
        return;
    }
    f.seekg(0, ios::beg);
    if (f.is_open())
        f.read(reinterpret_cast<char *>(var_pointer), data_size);
    else
        std::cerr << "Could not read '" << filename << "'" << std::endl;
    if (f.fail())
        std::cerr << "Error reading '" << filename << "'" << std::endl;
}

//////////////// set arrays by name ///////
void set_variable_by_name(std::string name, std::string s_value) {
    size_t var_size;
    size_t data_size;
    // C-style or Python-style capitalization is allowed for boolean values
    if (s_value == "true" || s_value == "True")
        s_value = "1";
    else if (s_value == "false" || s_value == "False")
        s_value = "0";
    // non-dynamic arrays
    if (name == "IN_group._spikespace") {
        var_size = 21;
        data_size = 21*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_IN_group__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.hNa_IN") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_hNa_IN, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_hNa_IN, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.iAMPA_PYso_IN") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_iAMPA_PYso_IN, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_iAMPA_PYso_IN, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.iAMPA_TC_IN") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_iAMPA_TC_IN, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_iAMPA_TC_IN, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.iGABAA_IN_IN") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_iGABAA_IN_IN, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_iGABAA_IN_IN, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.iNMDA_PYso_IN") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_iNMDA_PYso_IN, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_iNMDA_PYso_IN, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.lastspike") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_lastspike, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_lastspike, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.nK_IN") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_nK_IN, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_nK_IN, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.not_refractory") {
        var_size = 20;
        data_size = 20*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, _array_IN_group_not_refractory, var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_not_refractory, data_size, s_value);
        }
        return;
    }
    if (name == "IN_group.v") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_IN_group_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_IN_group_v, data_size, s_value);
        }
        return;
    }
    if (name == "poissongroup_1._spikespace") {
        var_size = 21;
        data_size = 21*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_poissongroup_1__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_poissongroup_1__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "poissongroup_1.rates") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_poissongroup_1_rates, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_poissongroup_1_rates, data_size, s_value);
        }
        return;
    }
    if (name == "poissongroup_2._spikespace") {
        var_size = 21;
        data_size = 21*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_poissongroup_2__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_poissongroup_2__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "poissongroup_2.rates") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_poissongroup_2_rates, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_poissongroup_2_rates, data_size, s_value);
        }
        return;
    }
    if (name == "poissongroup._spikespace") {
        var_size = 101;
        data_size = 101*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_poissongroup__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_poissongroup__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "poissongroup.rates") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_poissongroup_rates, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_poissongroup_rates, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group._spikespace") {
        var_size = 101;
        data_size = 101*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_PYdr_group__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.CaBuffer") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_CaBuffer, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_CaBuffer, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.gPoisson_PYdr") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_gPoisson_PYdr, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_gPoisson_PYdr, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.iAMPA_PYso_PYdr") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_iAMPA_PYso_PYdr, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_iAMPA_PYso_PYdr, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.iAMPA_TC_PYdr") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_iAMPA_TC_PYdr, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_iAMPA_TC_PYdr, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.iCOM") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_iCOM, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_iCOM, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.iNMDA_PYso_PYdr") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_iNMDA_PYso_PYdr, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_iNMDA_PYso_PYdr, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.lastspike") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_lastspike, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_lastspike, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.not_refractory") {
        var_size = 100;
        data_size = 100*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, _array_PYdr_group_not_refractory, var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_not_refractory, data_size, s_value);
        }
        return;
    }
    if (name == "PYdr_group.v") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYdr_group_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYdr_group_v, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group._spikespace") {
        var_size = 101;
        data_size = 101*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_PYso_group__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.hA") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_hA, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_hA, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.hNa") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_hNa, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_hNa, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.iCOM") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_iCOM, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_iCOM, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.iGABAA_IN_PYso") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_iGABAA_IN_PYso, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_iGABAA_IN_PYso, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.iKNa") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_iKNa, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_iKNa, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.lastspike") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_lastspike, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_lastspike, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.mKS") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_mKS, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_mKS, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.nK") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_nK, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_nK, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.not_refractory") {
        var_size = 100;
        data_size = 100*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, _array_PYso_group_not_refractory, var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_not_refractory, data_size, s_value);
        }
        return;
    }
    if (name == "PYso_group.v") {
        var_size = 100;
        data_size = 100*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_PYso_group_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_PYso_group_v, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group._spikespace") {
        var_size = 21;
        data_size = 21*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_RE_group__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.gPoisson_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_gPoisson_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_gPoisson_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.hNa_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_hNa_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_hNa_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.hT_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_hT_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_hT_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.iAMPA_PYso_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_iAMPA_PYso_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_iAMPA_PYso_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.iAMPA_TC_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_iAMPA_TC_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_iAMPA_TC_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.iGABAA_RE_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_iGABAA_RE_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_iGABAA_RE_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.lastspike") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_lastspike, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_lastspike, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.mNa_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_mNa_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_mNa_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.mT_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_mT_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_mT_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.nK_RE") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_nK_RE, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_nK_RE, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.not_refractory") {
        var_size = 20;
        data_size = 20*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, _array_RE_group_not_refractory, var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_not_refractory, data_size, s_value);
        }
        return;
    }
    if (name == "RE_group.v") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_RE_group_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_RE_group_v, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group._spikespace") {
        var_size = 21;
        data_size = 21*sizeof(int32_t);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<int32_t>(name, _array_TC_group__spikespace, var_size, (int32_t)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group__spikespace, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.Ca_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_Ca_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_Ca_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.dCa_extra") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_dCa_extra, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_dCa_extra, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.gPoisson_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_gPoisson_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_gPoisson_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.hNa_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_hNa_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_hNa_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.hT_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_hT_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_hT_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.iAMPA_PYso_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_iAMPA_PYso_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_iAMPA_PYso_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.iGABAA_RE_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_iGABAA_RE_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_iGABAA_RE_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.iGABAB_RE_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_iGABAB_RE_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_iGABAB_RE_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.lastspike") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_lastspike, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_lastspike, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.mNa_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_mNa_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_mNa_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.nK_TC") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_nK_TC, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_nK_TC, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.not_refractory") {
        var_size = 20;
        data_size = 20*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, _array_TC_group_not_refractory, var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_not_refractory, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.Open") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_Open, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_Open, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.OpenLocked") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_OpenLocked, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_OpenLocked, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.Pone") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_Pone, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_Pone, data_size, s_value);
        }
        return;
    }
    if (name == "TC_group.v") {
        var_size = 20;
        data_size = 20*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, _array_TC_group_v, var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, _array_TC_group_v, data_size, s_value);
        }
        return;
    }
    // dynamic arrays (1d)
    if (name == "IAMPA_PYso_IN.alpha_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_IN_alpha_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_alpha_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_alpha_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.delay") {
        var_size = _dynamic_array_IAMPA_PYso_IN_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.delay") {
        var_size = _dynamic_array_IAMPA_PYso_IN_delay_1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_delay_1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_delay_1[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.deprFactor") {
        var_size = _dynamic_array_IAMPA_PYso_IN_deprFactor.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_deprFactor[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_deprFactor[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.EAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_IN_EAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_EAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_EAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.gAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_IN_gAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_gAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_gAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.Npost") {
        var_size = _dynamic_array_IAMPA_PYso_IN_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.Npre") {
        var_size = _dynamic_array_IAMPA_PYso_IN_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.radius") {
        var_size = _dynamic_array_IAMPA_PYso_IN_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.remove_recurrent_bool") {
        var_size = _dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.res_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_IN_res_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_res_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_res_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.sAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_IN_sAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_sAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_sAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.spike_activity") {
        var_size = _dynamic_array_IAMPA_PYso_IN_spike_activity.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_spike_activity[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_spike_activity[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.tauAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_IN_tauAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_tauAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_tauAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_IN.tauRes_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_IN_tauRes_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_IN_tauRes_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.alpha_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.delay") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.delay") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_delay_1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_delay_1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_delay_1[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.deprFactor") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_deprFactor.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_deprFactor[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_deprFactor[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.EAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_EAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_EAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_EAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.gAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_gAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_gAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_gAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.Npost") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.Npre") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.radius") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.remove_recurrent_bool") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.res_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_res_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_res_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_res_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.sAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_sAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_sAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_sAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.spike_activity") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_spike_activity.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_spike_activity[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_spike_activity[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.tauAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_tauAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_tauAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_tauAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_PYdr.tauRes_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.EAMPA_PYso_RE") {
        var_size = _dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.gAMPA_PYso_RE") {
        var_size = _dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.Npost") {
        var_size = _dynamic_array_IAMPA_PYso_RE_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.Npre") {
        var_size = _dynamic_array_IAMPA_PYso_RE_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.P_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_RE_P_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_P_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_P_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.radius") {
        var_size = _dynamic_array_IAMPA_PYso_RE_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.remove_recurrent_bool") {
        var_size = _dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.sAMPA_PYso_RE") {
        var_size = _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_RE.tauAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_RE_tauAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_RE_tauAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_RE_tauAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.EAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_TC_EAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_EAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_EAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.gAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_TC_gAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_gAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_gAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.Npost") {
        var_size = _dynamic_array_IAMPA_PYso_TC_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.Npre") {
        var_size = _dynamic_array_IAMPA_PYso_TC_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.P_AMPA") {
        var_size = _dynamic_array_IAMPA_PYso_TC_P_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_P_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_P_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.radius") {
        var_size = _dynamic_array_IAMPA_PYso_TC_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.remove_recurrent_bool") {
        var_size = _dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.sAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_TC_sAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_sAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_sAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_PYso_TC.tauAMPA") {
        var_size = _dynamic_array_IAMPA_PYso_TC_tauAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_PYso_TC_tauAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_PYso_TC_tauAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.EAMPA") {
        var_size = _dynamic_array_IAMPA_TC_IN_EAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_EAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_EAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.gAMPA") {
        var_size = _dynamic_array_IAMPA_TC_IN_gAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_gAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_gAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.Npost") {
        var_size = _dynamic_array_IAMPA_TC_IN_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.Npre") {
        var_size = _dynamic_array_IAMPA_TC_IN_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.P_AMPA") {
        var_size = _dynamic_array_IAMPA_TC_IN_P_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_P_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_P_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.radius") {
        var_size = _dynamic_array_IAMPA_TC_IN_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.remove_recurrent_bool") {
        var_size = _dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IAMPA_TC_IN_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.sAMPA") {
        var_size = _dynamic_array_IAMPA_TC_IN_sAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_sAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_sAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_IN.tauAMPA") {
        var_size = _dynamic_array_IAMPA_TC_IN_tauAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_IN_tauAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_IN_tauAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.EAMPA") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_EAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_EAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_EAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.gAMPA") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_gAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_gAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_gAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.Npost") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.Npre") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.P_AMPA") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_P_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_P_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_P_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.radius") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.remove_recurrent_bool") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.sAMPA") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_sAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_sAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_sAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_PYdr.tauAMPA") {
        var_size = _dynamic_array_IAMPA_TC_PYdr_tauAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_PYdr_tauAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_PYdr_tauAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_RE.EAMPA") {
        var_size = _dynamic_array_IAMPA_TC_RE_EAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_RE_EAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_RE_EAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_RE.gAMPA") {
        var_size = _dynamic_array_IAMPA_TC_RE_gAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_RE_gAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_RE_gAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_RE.Npost") {
        var_size = _dynamic_array_IAMPA_TC_RE_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_RE_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_RE_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_RE.Npre") {
        var_size = _dynamic_array_IAMPA_TC_RE_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_RE_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_RE_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_RE.P_AMPA") {
        var_size = _dynamic_array_IAMPA_TC_RE_P_AMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_RE_P_AMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_RE_P_AMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_RE.sAMPA") {
        var_size = _dynamic_array_IAMPA_TC_RE_sAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_RE_sAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_RE_sAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IAMPA_TC_RE.tauAMPA") {
        var_size = _dynamic_array_IAMPA_TC_RE_tauAMPA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IAMPA_TC_RE_tauAMPA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IAMPA_TC_RE_tauAMPA[0], data_size, s_value);
        }
        return;
    }
    if (name == "ICOM_PYdr_PYso.gCOM") {
        var_size = _dynamic_array_ICOM_PYdr_PYso_gCOM.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_ICOM_PYdr_PYso_gCOM[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_ICOM_PYdr_PYso_gCOM[0], data_size, s_value);
        }
        return;
    }
    if (name == "ICOM_PYso_PYdr.gCOM") {
        var_size = _dynamic_array_ICOM_PYso_PYdr_gCOM.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_ICOM_PYso_PYdr_gCOM[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_ICOM_PYso_PYdr_gCOM[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.alpha_GABAA") {
        var_size = _dynamic_array_IGABAA_IN_IN_alpha_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_alpha_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_alpha_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.delay") {
        var_size = _dynamic_array_IGABAA_IN_IN_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.delay") {
        var_size = _dynamic_array_IGABAA_IN_IN_delay_1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_delay_1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_delay_1[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.deprFactor") {
        var_size = _dynamic_array_IGABAA_IN_IN_deprFactor.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_deprFactor[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_deprFactor[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.EGABAA") {
        var_size = _dynamic_array_IGABAA_IN_IN_EGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_EGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_EGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.gGABAA") {
        var_size = _dynamic_array_IGABAA_IN_IN_gGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_gGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_gGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.Npost") {
        var_size = _dynamic_array_IGABAA_IN_IN_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.Npre") {
        var_size = _dynamic_array_IGABAA_IN_IN_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.propoCondMult") {
        var_size = _dynamic_array_IGABAA_IN_IN_propoCondMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_propoCondMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_propoCondMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.propoTauMult") {
        var_size = _dynamic_array_IGABAA_IN_IN_propoTauMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_propoTauMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_propoTauMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.radius") {
        var_size = _dynamic_array_IGABAA_IN_IN_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.remove_recurrent_bool") {
        var_size = _dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IGABAA_IN_IN_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.res_GABAA") {
        var_size = _dynamic_array_IGABAA_IN_IN_res_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_res_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_res_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.sGABAA") {
        var_size = _dynamic_array_IGABAA_IN_IN_sGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_sGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_sGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.spike_activity") {
        var_size = _dynamic_array_IGABAA_IN_IN_spike_activity.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_spike_activity[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_spike_activity[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.tauGABAA") {
        var_size = _dynamic_array_IGABAA_IN_IN_tauGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_tauGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_tauGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_IN.tauRes_GABAA") {
        var_size = _dynamic_array_IGABAA_IN_IN_tauRes_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_IN_tauRes_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_IN_tauRes_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.alpha_GABAA") {
        var_size = _dynamic_array_IGABAA_IN_PYso_alpha_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_alpha_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_alpha_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.delay") {
        var_size = _dynamic_array_IGABAA_IN_PYso_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.delay") {
        var_size = _dynamic_array_IGABAA_IN_PYso_delay_1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_delay_1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_delay_1[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.deprFactor") {
        var_size = _dynamic_array_IGABAA_IN_PYso_deprFactor.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_deprFactor[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_deprFactor[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.EGABAA") {
        var_size = _dynamic_array_IGABAA_IN_PYso_EGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_EGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_EGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.gGABAA") {
        var_size = _dynamic_array_IGABAA_IN_PYso_gGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_gGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_gGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.Npost") {
        var_size = _dynamic_array_IGABAA_IN_PYso_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.Npre") {
        var_size = _dynamic_array_IGABAA_IN_PYso_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.propoCondMult") {
        var_size = _dynamic_array_IGABAA_IN_PYso_propoCondMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_propoCondMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_propoCondMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.propoTauMult") {
        var_size = _dynamic_array_IGABAA_IN_PYso_propoTauMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_propoTauMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_propoTauMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.radius") {
        var_size = _dynamic_array_IGABAA_IN_PYso_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.remove_recurrent_bool") {
        var_size = _dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.res_GABAA") {
        var_size = _dynamic_array_IGABAA_IN_PYso_res_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_res_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_res_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.sGABAA") {
        var_size = _dynamic_array_IGABAA_IN_PYso_sGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_sGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_sGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.spike_activity") {
        var_size = _dynamic_array_IGABAA_IN_PYso_spike_activity.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_spike_activity[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_spike_activity[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.tauGABAA") {
        var_size = _dynamic_array_IGABAA_IN_PYso_tauGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_tauGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_tauGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_IN_PYso.tauRes_GABAA") {
        var_size = _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_IN_PYso_tauRes_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_IN_PYso_tauRes_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.EGABAA") {
        var_size = _dynamic_array_IGABAA_RE_RE_EGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_EGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_EGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.gGABAA") {
        var_size = _dynamic_array_IGABAA_RE_RE_gGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_gGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_gGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.Npost") {
        var_size = _dynamic_array_IGABAA_RE_RE_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.Npre") {
        var_size = _dynamic_array_IGABAA_RE_RE_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.P_GABAA") {
        var_size = _dynamic_array_IGABAA_RE_RE_P_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_P_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_P_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.propoCondMult") {
        var_size = _dynamic_array_IGABAA_RE_RE_propoCondMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_propoCondMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_propoCondMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.propoTauMult") {
        var_size = _dynamic_array_IGABAA_RE_RE_propoTauMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_propoTauMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_propoTauMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.sGABAA") {
        var_size = _dynamic_array_IGABAA_RE_RE_sGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_sGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_sGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_RE.tauGABAA") {
        var_size = _dynamic_array_IGABAA_RE_RE_tauGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_RE_tauGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_RE_tauGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.EGABAA") {
        var_size = _dynamic_array_IGABAA_RE_TC_EGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_EGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_EGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.gGABAA") {
        var_size = _dynamic_array_IGABAA_RE_TC_gGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_gGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_gGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.Npost") {
        var_size = _dynamic_array_IGABAA_RE_TC_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.Npre") {
        var_size = _dynamic_array_IGABAA_RE_TC_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.P_GABAA") {
        var_size = _dynamic_array_IGABAA_RE_TC_P_GABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_P_GABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_P_GABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.propoCondMult") {
        var_size = _dynamic_array_IGABAA_RE_TC_propoCondMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_propoCondMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_propoCondMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.propoTauMult") {
        var_size = _dynamic_array_IGABAA_RE_TC_propoTauMult.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_propoTauMult[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_propoTauMult[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.sGABAA") {
        var_size = _dynamic_array_IGABAA_RE_TC_sGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_sGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_sGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAA_RE_TC.tauGABAA") {
        var_size = _dynamic_array_IGABAA_RE_TC_tauGABAA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAA_RE_TC_tauGABAA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAA_RE_TC_tauGABAA[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.EGABAB") {
        var_size = _dynamic_array_IGABAB_RE_TC_EGABAB.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_EGABAB[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_EGABAB[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.gGABAB") {
        var_size = _dynamic_array_IGABAB_RE_TC_gGABAB.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_gGABAB[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_gGABAB[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.K1") {
        var_size = _dynamic_array_IGABAB_RE_TC_K1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_K1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_K1[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.K2") {
        var_size = _dynamic_array_IGABAB_RE_TC_K2.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_K2[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_K2[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.K3") {
        var_size = _dynamic_array_IGABAB_RE_TC_K3.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_K3[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_K3[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.K4") {
        var_size = _dynamic_array_IGABAB_RE_TC_K4.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_K4[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_K4[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.Npost") {
        var_size = _dynamic_array_IGABAB_RE_TC_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.Npre") {
        var_size = _dynamic_array_IGABAB_RE_TC_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.rGABAB") {
        var_size = _dynamic_array_IGABAB_RE_TC_rGABAB.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_rGABAB[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_rGABAB[0], data_size, s_value);
        }
        return;
    }
    if (name == "IGABAB_RE_TC.sGABAB") {
        var_size = _dynamic_array_IGABAB_RE_TC_sGABAB.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IGABAB_RE_TC_sGABAB[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IGABAB_RE_TC_sGABAB[0], data_size, s_value);
        }
        return;
    }
    if (name == "IKNa_PYso_PYdr.concNa") {
        var_size = _dynamic_array_IKNa_PYso_PYdr_concNa.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IKNa_PYso_PYdr_concNa[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IKNa_PYso_PYdr_concNa[0], data_size, s_value);
        }
        return;
    }
    if (name == "IKNa_PYso_PYdr.hNa_local") {
        var_size = _dynamic_array_IKNa_PYso_PYdr_hNa_local.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_IKNa_PYso_PYdr_hNa_local[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_IKNa_PYso_PYdr_hNa_local[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.alphaS_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_alphaS_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_alphaS_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_alphaS_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.alphaX_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_alphaX_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_alphaX_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_alphaX_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.delay") {
        var_size = _dynamic_array_INMDA_PYso_IN_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.delay") {
        var_size = _dynamic_array_INMDA_PYso_IN_delay_1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_delay_1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_delay_1[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.deprFactor") {
        var_size = _dynamic_array_INMDA_PYso_IN_deprFactor.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_deprFactor[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_deprFactor[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.ENMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_ENMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_ENMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_ENMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.gNMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_gNMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_gNMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_gNMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.Npost") {
        var_size = _dynamic_array_INMDA_PYso_IN_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.Npre") {
        var_size = _dynamic_array_INMDA_PYso_IN_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.radius") {
        var_size = _dynamic_array_INMDA_PYso_IN_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.remove_recurrent_bool") {
        var_size = _dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_INMDA_PYso_IN_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.res_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_res_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_res_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_res_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.sNMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_sNMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_sNMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_sNMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.spike_activity") {
        var_size = _dynamic_array_INMDA_PYso_IN_spike_activity.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_spike_activity[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_spike_activity[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.tauRes_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_tauRes_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_tauRes_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_tauRes_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.tauS_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_tauS_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_tauS_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_tauS_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.tauX_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_tauX_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_tauX_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_tauX_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_IN.xNMDA") {
        var_size = _dynamic_array_INMDA_PYso_IN_xNMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_IN_xNMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_IN_xNMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.alphaS_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.alphaX_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.delay") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.delay") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_delay_1.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_delay_1[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_delay_1[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.deprFactor") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_deprFactor.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_deprFactor[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_deprFactor[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.ENMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_ENMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_ENMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_ENMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.gNMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_gNMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_gNMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_gNMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.Npost") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_Npost.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_Npost[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_Npost[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.Npre") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_Npre.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_Npre[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_Npre[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.radius") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_radius.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_radius[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_radius[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.remove_recurrent_bool") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.size();
        data_size = var_size*sizeof(char);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value(name, &_dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool[0], var_size, (char)atoi(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.res_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_res_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_res_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_res_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.sNMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_sNMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_sNMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_sNMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.spike_activity") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_spike_activity.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_spike_activity[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_spike_activity[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.tauRes_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.tauS_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_tauS_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_tauS_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.tauX_NMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_tauX_NMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_tauX_NMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "INMDA_PYso_PYdr.xNMDA") {
        var_size = _dynamic_array_INMDA_PYso_PYdr_xNMDA.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_INMDA_PYso_PYdr_xNMDA[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_INMDA_PYso_PYdr_xNMDA[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses_1.delay") {
        var_size = _dynamic_array_synapses_1_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_1_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_1_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses_2.delay") {
        var_size = _dynamic_array_synapses_2_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_2_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_2_delay[0], data_size, s_value);
        }
        return;
    }
    if (name == "synapses.delay") {
        var_size = _dynamic_array_synapses_delay.size();
        data_size = var_size*sizeof(double);
        if (s_value[0] == '-' || (s_value[0] >= '0' && s_value[0] <= '9')) {
            // set from single value
            set_variable_from_value<double>(name, &_dynamic_array_synapses_delay[0], var_size, (double)atof(s_value.c_str()));

        } else {
            // set from file
            set_variable_from_file(name, &_dynamic_array_synapses_delay[0], data_size, s_value);
        }
        return;
    }
    std::cerr << "Cannot set unknown variable '" << name << "'." << std::endl;
    exit(1);
}
//////////////// arrays ///////////////////
double * _array_clock_dt;
const int _num__array_clock_dt = 1;
double * _array_clock_t;
const int _num__array_clock_t = 1;
int64_t * _array_clock_timestep;
const int _num__array_clock_timestep = 1;
double * _array_defaultclock_dt;
const int _num__array_defaultclock_dt = 1;
double * _array_defaultclock_t;
const int _num__array_defaultclock_t = 1;
int64_t * _array_defaultclock_timestep;
const int _num__array_defaultclock_timestep = 1;
int32_t * _array_IAMPA_PYso_IN_N;
const int _num__array_IAMPA_PYso_IN_N = 1;
int32_t * _array_IAMPA_PYso_IN_sources;
const int _num__array_IAMPA_PYso_IN_sources = 2000;
int32_t * _array_IAMPA_PYso_IN_targets;
const int _num__array_IAMPA_PYso_IN_targets = 2000;
int32_t * _array_IAMPA_PYso_PYdr_N;
const int _num__array_IAMPA_PYso_PYdr_N = 1;
int32_t * _array_IAMPA_PYso_PYdr_sources;
const int _num__array_IAMPA_PYso_PYdr_sources = 2000;
int32_t * _array_IAMPA_PYso_PYdr_targets;
const int _num__array_IAMPA_PYso_PYdr_targets = 2000;
int32_t * _array_IAMPA_PYso_RE_N;
const int _num__array_IAMPA_PYso_RE_N = 1;
int32_t * _array_IAMPA_PYso_RE_sources;
const int _num__array_IAMPA_PYso_RE_sources = 2000;
int32_t * _array_IAMPA_PYso_RE_targets;
const int _num__array_IAMPA_PYso_RE_targets = 2000;
int32_t * _array_IAMPA_PYso_TC_N;
const int _num__array_IAMPA_PYso_TC_N = 1;
int32_t * _array_IAMPA_PYso_TC_sources;
const int _num__array_IAMPA_PYso_TC_sources = 2000;
int32_t * _array_IAMPA_PYso_TC_targets;
const int _num__array_IAMPA_PYso_TC_targets = 2000;
int32_t * _array_IAMPA_TC_IN_N;
const int _num__array_IAMPA_TC_IN_N = 1;
int32_t * _array_IAMPA_TC_IN_sources;
const int _num__array_IAMPA_TC_IN_sources = 400;
int32_t * _array_IAMPA_TC_IN_targets;
const int _num__array_IAMPA_TC_IN_targets = 400;
int32_t * _array_IAMPA_TC_PYdr_N;
const int _num__array_IAMPA_TC_PYdr_N = 1;
int32_t * _array_IAMPA_TC_PYdr_sources;
const int _num__array_IAMPA_TC_PYdr_sources = 420;
int32_t * _array_IAMPA_TC_PYdr_targets;
const int _num__array_IAMPA_TC_PYdr_targets = 420;
int32_t * _array_IAMPA_TC_RE_N;
const int _num__array_IAMPA_TC_RE_N = 1;
int32_t * _array_ICOM_PYdr_PYso_N;
const int _num__array_ICOM_PYdr_PYso_N = 1;
int32_t * _array_ICOM_PYso_PYdr_N;
const int _num__array_ICOM_PYso_PYdr_N = 1;
int32_t * _array_IGABAA_IN_IN_N;
const int _num__array_IGABAA_IN_IN_N = 1;
int32_t * _array_IGABAA_IN_IN_sources;
const int _num__array_IGABAA_IN_IN_sources = 380;
int32_t * _array_IGABAA_IN_IN_targets;
const int _num__array_IGABAA_IN_IN_targets = 380;
int32_t * _array_IGABAA_IN_PYso_N;
const int _num__array_IGABAA_IN_PYso_N = 1;
int32_t * _array_IGABAA_IN_PYso_sources;
const int _num__array_IGABAA_IN_PYso_sources = 420;
int32_t * _array_IGABAA_IN_PYso_targets;
const int _num__array_IGABAA_IN_PYso_targets = 420;
int32_t * _array_IGABAA_RE_RE_N;
const int _num__array_IGABAA_RE_RE_N = 1;
int32_t * _array_IGABAA_RE_TC_N;
const int _num__array_IGABAA_RE_TC_N = 1;
int32_t * _array_IGABAB_RE_TC_N;
const int _num__array_IGABAB_RE_TC_N = 1;
int32_t * _array_IKNa_PYso_PYdr_N;
const int _num__array_IKNa_PYso_PYdr_N = 1;
int32_t * _array_IN_group__spikespace;
const int _num__array_IN_group__spikespace = 21;
double * _array_IN_group_hNa_IN;
const int _num__array_IN_group_hNa_IN = 20;
int32_t * _array_IN_group_i;
const int _num__array_IN_group_i = 20;
double * _array_IN_group_iAMPA_PYso_IN;
const int _num__array_IN_group_iAMPA_PYso_IN = 20;
double * _array_IN_group_iAMPA_TC_IN;
const int _num__array_IN_group_iAMPA_TC_IN = 20;
double * _array_IN_group_iGABAA_IN_IN;
const int _num__array_IN_group_iGABAA_IN_IN = 20;
double * _array_IN_group_iNMDA_PYso_IN;
const int _num__array_IN_group_iNMDA_PYso_IN = 20;
double * _array_IN_group_lastspike;
const int _num__array_IN_group_lastspike = 20;
double * _array_IN_group_nK_IN;
const int _num__array_IN_group_nK_IN = 20;
char * _array_IN_group_not_refractory;
const int _num__array_IN_group_not_refractory = 20;
double * _array_IN_group_v;
const int _num__array_IN_group_v = 20;
int32_t * _array_IN_spikemon__source_idx;
const int _num__array_IN_spikemon__source_idx = 20;
int32_t * _array_IN_spikemon_count;
const int _num__array_IN_spikemon_count = 20;
int32_t * _array_IN_spikemon_N;
const int _num__array_IN_spikemon_N = 1;
int32_t * _array_INMDA_PYso_IN_N;
const int _num__array_INMDA_PYso_IN_N = 1;
int32_t * _array_INMDA_PYso_IN_sources;
const int _num__array_INMDA_PYso_IN_sources = 2000;
int32_t * _array_INMDA_PYso_IN_targets;
const int _num__array_INMDA_PYso_IN_targets = 2000;
int32_t * _array_INMDA_PYso_PYdr_N;
const int _num__array_INMDA_PYso_PYdr_N = 1;
int32_t * _array_INMDA_PYso_PYdr_sources;
const int _num__array_INMDA_PYso_PYdr_sources = 2000;
int32_t * _array_INMDA_PYso_PYdr_targets;
const int _num__array_INMDA_PYso_PYdr_targets = 2000;
int32_t * _array_poissongroup_1__spikespace;
const int _num__array_poissongroup_1__spikespace = 21;
int32_t * _array_poissongroup_1_i;
const int _num__array_poissongroup_1_i = 20;
double * _array_poissongroup_1_rates;
const int _num__array_poissongroup_1_rates = 20;
int32_t * _array_poissongroup_2__spikespace;
const int _num__array_poissongroup_2__spikespace = 21;
int32_t * _array_poissongroup_2_i;
const int _num__array_poissongroup_2_i = 20;
double * _array_poissongroup_2_rates;
const int _num__array_poissongroup_2_rates = 20;
int32_t * _array_poissongroup__spikespace;
const int _num__array_poissongroup__spikespace = 101;
int32_t * _array_poissongroup_i;
const int _num__array_poissongroup_i = 100;
double * _array_poissongroup_rates;
const int _num__array_poissongroup_rates = 100;
int32_t * _array_PYdr_group__spikespace;
const int _num__array_PYdr_group__spikespace = 101;
double * _array_PYdr_group_CaBuffer;
const int _num__array_PYdr_group_CaBuffer = 100;
double * _array_PYdr_group_gPoisson_PYdr;
const int _num__array_PYdr_group_gPoisson_PYdr = 100;
int32_t * _array_PYdr_group_i;
const int _num__array_PYdr_group_i = 100;
double * _array_PYdr_group_iAMPA_PYso_PYdr;
const int _num__array_PYdr_group_iAMPA_PYso_PYdr = 100;
double * _array_PYdr_group_iAMPA_TC_PYdr;
const int _num__array_PYdr_group_iAMPA_TC_PYdr = 100;
double * _array_PYdr_group_iCOM;
const int _num__array_PYdr_group_iCOM = 100;
double * _array_PYdr_group_iNMDA_PYso_PYdr;
const int _num__array_PYdr_group_iNMDA_PYso_PYdr = 100;
double * _array_PYdr_group_lastspike;
const int _num__array_PYdr_group_lastspike = 100;
char * _array_PYdr_group_not_refractory;
const int _num__array_PYdr_group_not_refractory = 100;
double * _array_PYdr_group_v;
const int _num__array_PYdr_group_v = 100;
int32_t * _array_PYdr_spikemon__source_idx;
const int _num__array_PYdr_spikemon__source_idx = 100;
int32_t * _array_PYdr_spikemon_count;
const int _num__array_PYdr_spikemon_count = 100;
int32_t * _array_PYdr_spikemon_N;
const int _num__array_PYdr_spikemon_N = 1;
int32_t * _array_PYso_group__spikespace;
const int _num__array_PYso_group__spikespace = 101;
double * _array_PYso_group_hA;
const int _num__array_PYso_group_hA = 100;
double * _array_PYso_group_hNa;
const int _num__array_PYso_group_hNa = 100;
int32_t * _array_PYso_group_i;
const int _num__array_PYso_group_i = 100;
double * _array_PYso_group_iCOM;
const int _num__array_PYso_group_iCOM = 100;
double * _array_PYso_group_iGABAA_IN_PYso;
const int _num__array_PYso_group_iGABAA_IN_PYso = 100;
double * _array_PYso_group_iKNa;
const int _num__array_PYso_group_iKNa = 100;
double * _array_PYso_group_lastspike;
const int _num__array_PYso_group_lastspike = 100;
double * _array_PYso_group_mKS;
const int _num__array_PYso_group_mKS = 100;
double * _array_PYso_group_nK;
const int _num__array_PYso_group_nK = 100;
char * _array_PYso_group_not_refractory;
const int _num__array_PYso_group_not_refractory = 100;
double * _array_PYso_group_v;
const int _num__array_PYso_group_v = 100;
int32_t * _array_PYso_spikemon__source_idx;
const int _num__array_PYso_spikemon__source_idx = 100;
int32_t * _array_PYso_spikemon_count;
const int _num__array_PYso_spikemon_count = 100;
int32_t * _array_PYso_spikemon_N;
const int _num__array_PYso_spikemon_N = 1;
int32_t * _array_RE_group__spikespace;
const int _num__array_RE_group__spikespace = 21;
double * _array_RE_group_gPoisson_RE;
const int _num__array_RE_group_gPoisson_RE = 20;
double * _array_RE_group_hNa_RE;
const int _num__array_RE_group_hNa_RE = 20;
double * _array_RE_group_hT_RE;
const int _num__array_RE_group_hT_RE = 20;
int32_t * _array_RE_group_i;
const int _num__array_RE_group_i = 20;
double * _array_RE_group_iAMPA_PYso_RE;
const int _num__array_RE_group_iAMPA_PYso_RE = 20;
double * _array_RE_group_iAMPA_TC_RE;
const int _num__array_RE_group_iAMPA_TC_RE = 20;
double * _array_RE_group_iGABAA_RE_RE;
const int _num__array_RE_group_iGABAA_RE_RE = 20;
double * _array_RE_group_lastspike;
const int _num__array_RE_group_lastspike = 20;
double * _array_RE_group_mNa_RE;
const int _num__array_RE_group_mNa_RE = 20;
double * _array_RE_group_mT_RE;
const int _num__array_RE_group_mT_RE = 20;
double * _array_RE_group_nK_RE;
const int _num__array_RE_group_nK_RE = 20;
char * _array_RE_group_not_refractory;
const int _num__array_RE_group_not_refractory = 20;
double * _array_RE_group_v;
const int _num__array_RE_group_v = 20;
int32_t * _array_RE_spikemon__source_idx;
const int _num__array_RE_spikemon__source_idx = 20;
int32_t * _array_RE_spikemon_count;
const int _num__array_RE_spikemon_count = 20;
int32_t * _array_RE_spikemon_N;
const int _num__array_RE_spikemon_N = 1;
int32_t * _array_synapses_1_N;
const int _num__array_synapses_1_N = 1;
int32_t * _array_synapses_2_N;
const int _num__array_synapses_2_N = 1;
int32_t * _array_synapses_N;
const int _num__array_synapses_N = 1;
int32_t * _array_TC_group__spikespace;
const int _num__array_TC_group__spikespace = 21;
double * _array_TC_group_Ca_TC;
const int _num__array_TC_group_Ca_TC = 20;
double * _array_TC_group_dCa_extra;
const int _num__array_TC_group_dCa_extra = 20;
double * _array_TC_group_gPoisson_TC;
const int _num__array_TC_group_gPoisson_TC = 20;
double * _array_TC_group_hNa_TC;
const int _num__array_TC_group_hNa_TC = 20;
double * _array_TC_group_hT_TC;
const int _num__array_TC_group_hT_TC = 20;
int32_t * _array_TC_group_i;
const int _num__array_TC_group_i = 20;
double * _array_TC_group_iAMPA_PYso_TC;
const int _num__array_TC_group_iAMPA_PYso_TC = 20;
double * _array_TC_group_iGABAA_RE_TC;
const int _num__array_TC_group_iGABAA_RE_TC = 20;
double * _array_TC_group_iGABAB_RE_TC;
const int _num__array_TC_group_iGABAB_RE_TC = 20;
double * _array_TC_group_lastspike;
const int _num__array_TC_group_lastspike = 20;
double * _array_TC_group_mNa_TC;
const int _num__array_TC_group_mNa_TC = 20;
double * _array_TC_group_nK_TC;
const int _num__array_TC_group_nK_TC = 20;
char * _array_TC_group_not_refractory;
const int _num__array_TC_group_not_refractory = 20;
double * _array_TC_group_Open;
const int _num__array_TC_group_Open = 20;
double * _array_TC_group_OpenLocked;
const int _num__array_TC_group_OpenLocked = 20;
double * _array_TC_group_Pone;
const int _num__array_TC_group_Pone = 20;
double * _array_TC_group_v;
const int _num__array_TC_group_v = 20;
int32_t * _array_TC_spikemon__source_idx;
const int _num__array_TC_spikemon__source_idx = 20;
int32_t * _array_TC_spikemon_count;
const int _num__array_TC_spikemon_count = 20;
int32_t * _array_TC_spikemon_N;
const int _num__array_TC_spikemon_N = 1;

//////////////// dynamic arrays 1d /////////
std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN__synaptic_post;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN__synaptic_pre;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_alpha_AMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_delay;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_delay_1;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_deprFactor;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_EAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_gAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN_N_incoming;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_IN_N_outgoing;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_Npost;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_Npre;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_radius;
std::vector<char> _dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_res_AMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_sAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_spike_activity;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_tauAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_IN_tauRes_AMPA;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr__synaptic_post;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr__synaptic_pre;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_delay;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_delay_1;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_deprFactor;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_EAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_gAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr_N_incoming;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_PYdr_N_outgoing;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_Npost;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_Npre;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_radius;
std::vector<char> _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_sAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_spike_activity;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_tauAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE__synaptic_post;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE__synaptic_pre;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE_N_incoming;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_RE_N_outgoing;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_Npost;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_Npre;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_P_AMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_radius;
std::vector<char> _dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
std::vector<double> _dynamic_array_IAMPA_PYso_RE_tauAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC__synaptic_post;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC__synaptic_pre;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_EAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_gAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC_N_incoming;
std::vector<int32_t> _dynamic_array_IAMPA_PYso_TC_N_outgoing;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_Npost;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_Npre;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_P_AMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_radius;
std::vector<char> _dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_sAMPA;
std::vector<double> _dynamic_array_IAMPA_PYso_TC_tauAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_TC_IN__synaptic_post;
std::vector<int32_t> _dynamic_array_IAMPA_TC_IN__synaptic_pre;
std::vector<double> _dynamic_array_IAMPA_TC_IN_EAMPA;
std::vector<double> _dynamic_array_IAMPA_TC_IN_gAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_TC_IN_N_incoming;
std::vector<int32_t> _dynamic_array_IAMPA_TC_IN_N_outgoing;
std::vector<double> _dynamic_array_IAMPA_TC_IN_Npost;
std::vector<double> _dynamic_array_IAMPA_TC_IN_Npre;
std::vector<double> _dynamic_array_IAMPA_TC_IN_P_AMPA;
std::vector<double> _dynamic_array_IAMPA_TC_IN_radius;
std::vector<char> _dynamic_array_IAMPA_TC_IN_remove_recurrent_bool;
std::vector<double> _dynamic_array_IAMPA_TC_IN_sAMPA;
std::vector<double> _dynamic_array_IAMPA_TC_IN_tauAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr__synaptic_post;
std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr__synaptic_pre;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_EAMPA;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_gAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr_N_incoming;
std::vector<int32_t> _dynamic_array_IAMPA_TC_PYdr_N_outgoing;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_Npost;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_Npre;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_P_AMPA;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_radius;
std::vector<char> _dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_sAMPA;
std::vector<double> _dynamic_array_IAMPA_TC_PYdr_tauAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_TC_RE__synaptic_post;
std::vector<int32_t> _dynamic_array_IAMPA_TC_RE__synaptic_pre;
std::vector<double> _dynamic_array_IAMPA_TC_RE_EAMPA;
std::vector<double> _dynamic_array_IAMPA_TC_RE_gAMPA;
std::vector<int32_t> _dynamic_array_IAMPA_TC_RE_N_incoming;
std::vector<int32_t> _dynamic_array_IAMPA_TC_RE_N_outgoing;
std::vector<double> _dynamic_array_IAMPA_TC_RE_Npost;
std::vector<double> _dynamic_array_IAMPA_TC_RE_Npre;
std::vector<double> _dynamic_array_IAMPA_TC_RE_P_AMPA;
std::vector<double> _dynamic_array_IAMPA_TC_RE_sAMPA;
std::vector<double> _dynamic_array_IAMPA_TC_RE_tauAMPA;
std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso__synaptic_post;
std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso__synaptic_pre;
std::vector<double> _dynamic_array_ICOM_PYdr_PYso_gCOM;
std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso_N_incoming;
std::vector<int32_t> _dynamic_array_ICOM_PYdr_PYso_N_outgoing;
std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr__synaptic_post;
std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr__synaptic_pre;
std::vector<double> _dynamic_array_ICOM_PYso_PYdr_gCOM;
std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr_N_incoming;
std::vector<int32_t> _dynamic_array_ICOM_PYso_PYdr_N_outgoing;
std::vector<int32_t> _dynamic_array_IGABAA_IN_IN__synaptic_post;
std::vector<int32_t> _dynamic_array_IGABAA_IN_IN__synaptic_pre;
std::vector<double> _dynamic_array_IGABAA_IN_IN_alpha_GABAA;
std::vector<double> _dynamic_array_IGABAA_IN_IN_delay;
std::vector<double> _dynamic_array_IGABAA_IN_IN_delay_1;
std::vector<double> _dynamic_array_IGABAA_IN_IN_deprFactor;
std::vector<double> _dynamic_array_IGABAA_IN_IN_EGABAA;
std::vector<double> _dynamic_array_IGABAA_IN_IN_gGABAA;
std::vector<int32_t> _dynamic_array_IGABAA_IN_IN_N_incoming;
std::vector<int32_t> _dynamic_array_IGABAA_IN_IN_N_outgoing;
std::vector<double> _dynamic_array_IGABAA_IN_IN_Npost;
std::vector<double> _dynamic_array_IGABAA_IN_IN_Npre;
std::vector<double> _dynamic_array_IGABAA_IN_IN_propoCondMult;
std::vector<double> _dynamic_array_IGABAA_IN_IN_propoTauMult;
std::vector<double> _dynamic_array_IGABAA_IN_IN_radius;
std::vector<char> _dynamic_array_IGABAA_IN_IN_remove_recurrent_bool;
std::vector<double> _dynamic_array_IGABAA_IN_IN_res_GABAA;
std::vector<double> _dynamic_array_IGABAA_IN_IN_sGABAA;
std::vector<double> _dynamic_array_IGABAA_IN_IN_spike_activity;
std::vector<double> _dynamic_array_IGABAA_IN_IN_tauGABAA;
std::vector<double> _dynamic_array_IGABAA_IN_IN_tauRes_GABAA;
std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso__synaptic_post;
std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso__synaptic_pre;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_alpha_GABAA;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_delay;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_delay_1;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_deprFactor;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_EGABAA;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_gGABAA;
std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso_N_incoming;
std::vector<int32_t> _dynamic_array_IGABAA_IN_PYso_N_outgoing;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_Npost;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_Npre;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_propoCondMult;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_propoTauMult;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_radius;
std::vector<char> _dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_res_GABAA;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_sGABAA;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_spike_activity;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_tauGABAA;
std::vector<double> _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA;
std::vector<int32_t> _dynamic_array_IGABAA_RE_RE__synaptic_post;
std::vector<int32_t> _dynamic_array_IGABAA_RE_RE__synaptic_pre;
std::vector<double> _dynamic_array_IGABAA_RE_RE_EGABAA;
std::vector<double> _dynamic_array_IGABAA_RE_RE_gGABAA;
std::vector<int32_t> _dynamic_array_IGABAA_RE_RE_N_incoming;
std::vector<int32_t> _dynamic_array_IGABAA_RE_RE_N_outgoing;
std::vector<double> _dynamic_array_IGABAA_RE_RE_Npost;
std::vector<double> _dynamic_array_IGABAA_RE_RE_Npre;
std::vector<double> _dynamic_array_IGABAA_RE_RE_P_GABAA;
std::vector<double> _dynamic_array_IGABAA_RE_RE_propoCondMult;
std::vector<double> _dynamic_array_IGABAA_RE_RE_propoTauMult;
std::vector<double> _dynamic_array_IGABAA_RE_RE_sGABAA;
std::vector<double> _dynamic_array_IGABAA_RE_RE_tauGABAA;
std::vector<int32_t> _dynamic_array_IGABAA_RE_TC__synaptic_post;
std::vector<int32_t> _dynamic_array_IGABAA_RE_TC__synaptic_pre;
std::vector<double> _dynamic_array_IGABAA_RE_TC_EGABAA;
std::vector<double> _dynamic_array_IGABAA_RE_TC_gGABAA;
std::vector<int32_t> _dynamic_array_IGABAA_RE_TC_N_incoming;
std::vector<int32_t> _dynamic_array_IGABAA_RE_TC_N_outgoing;
std::vector<double> _dynamic_array_IGABAA_RE_TC_Npost;
std::vector<double> _dynamic_array_IGABAA_RE_TC_Npre;
std::vector<double> _dynamic_array_IGABAA_RE_TC_P_GABAA;
std::vector<double> _dynamic_array_IGABAA_RE_TC_propoCondMult;
std::vector<double> _dynamic_array_IGABAA_RE_TC_propoTauMult;
std::vector<double> _dynamic_array_IGABAA_RE_TC_sGABAA;
std::vector<double> _dynamic_array_IGABAA_RE_TC_tauGABAA;
std::vector<int32_t> _dynamic_array_IGABAB_RE_TC__synaptic_post;
std::vector<int32_t> _dynamic_array_IGABAB_RE_TC__synaptic_pre;
std::vector<double> _dynamic_array_IGABAB_RE_TC_EGABAB;
std::vector<double> _dynamic_array_IGABAB_RE_TC_gGABAB;
std::vector<double> _dynamic_array_IGABAB_RE_TC_K1;
std::vector<double> _dynamic_array_IGABAB_RE_TC_K2;
std::vector<double> _dynamic_array_IGABAB_RE_TC_K3;
std::vector<double> _dynamic_array_IGABAB_RE_TC_K4;
std::vector<int32_t> _dynamic_array_IGABAB_RE_TC_N_incoming;
std::vector<int32_t> _dynamic_array_IGABAB_RE_TC_N_outgoing;
std::vector<double> _dynamic_array_IGABAB_RE_TC_Npost;
std::vector<double> _dynamic_array_IGABAB_RE_TC_Npre;
std::vector<double> _dynamic_array_IGABAB_RE_TC_rGABAB;
std::vector<double> _dynamic_array_IGABAB_RE_TC_sGABAB;
std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr__synaptic_post;
std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr__synaptic_pre;
std::vector<double> _dynamic_array_IKNa_PYso_PYdr_concNa;
std::vector<double> _dynamic_array_IKNa_PYso_PYdr_hNa_local;
std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr_N_incoming;
std::vector<int32_t> _dynamic_array_IKNa_PYso_PYdr_N_outgoing;
std::vector<int32_t> _dynamic_array_IN_spikemon_i;
std::vector<double> _dynamic_array_IN_spikemon_t;
std::vector<int32_t> _dynamic_array_INMDA_PYso_IN__synaptic_post;
std::vector<int32_t> _dynamic_array_INMDA_PYso_IN__synaptic_pre;
std::vector<double> _dynamic_array_INMDA_PYso_IN_alphaS_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_alphaX_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_delay;
std::vector<double> _dynamic_array_INMDA_PYso_IN_delay_1;
std::vector<double> _dynamic_array_INMDA_PYso_IN_deprFactor;
std::vector<double> _dynamic_array_INMDA_PYso_IN_ENMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_gNMDA;
std::vector<int32_t> _dynamic_array_INMDA_PYso_IN_N_incoming;
std::vector<int32_t> _dynamic_array_INMDA_PYso_IN_N_outgoing;
std::vector<double> _dynamic_array_INMDA_PYso_IN_Npost;
std::vector<double> _dynamic_array_INMDA_PYso_IN_Npre;
std::vector<double> _dynamic_array_INMDA_PYso_IN_radius;
std::vector<char> _dynamic_array_INMDA_PYso_IN_remove_recurrent_bool;
std::vector<double> _dynamic_array_INMDA_PYso_IN_res_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_sNMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_spike_activity;
std::vector<double> _dynamic_array_INMDA_PYso_IN_tauRes_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_tauS_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_tauX_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_IN_xNMDA;
std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr__synaptic_post;
std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr__synaptic_pre;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_delay;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_delay_1;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_deprFactor;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_ENMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_gNMDA;
std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr_N_incoming;
std::vector<int32_t> _dynamic_array_INMDA_PYso_PYdr_N_outgoing;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_Npost;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_Npre;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_radius;
std::vector<char> _dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_res_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_sNMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_spike_activity;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA;
std::vector<double> _dynamic_array_INMDA_PYso_PYdr_xNMDA;
std::vector<int32_t> _dynamic_array_PYdr_spikemon_i;
std::vector<double> _dynamic_array_PYdr_spikemon_t;
std::vector<int32_t> _dynamic_array_PYso_spikemon_i;
std::vector<double> _dynamic_array_PYso_spikemon_t;
std::vector<int32_t> _dynamic_array_RE_spikemon_i;
std::vector<double> _dynamic_array_RE_spikemon_t;
std::vector<int32_t> _dynamic_array_synapses_1__synaptic_post;
std::vector<int32_t> _dynamic_array_synapses_1__synaptic_pre;
std::vector<double> _dynamic_array_synapses_1_delay;
std::vector<int32_t> _dynamic_array_synapses_1_N_incoming;
std::vector<int32_t> _dynamic_array_synapses_1_N_outgoing;
std::vector<int32_t> _dynamic_array_synapses_2__synaptic_post;
std::vector<int32_t> _dynamic_array_synapses_2__synaptic_pre;
std::vector<double> _dynamic_array_synapses_2_delay;
std::vector<int32_t> _dynamic_array_synapses_2_N_incoming;
std::vector<int32_t> _dynamic_array_synapses_2_N_outgoing;
std::vector<int32_t> _dynamic_array_synapses__synaptic_post;
std::vector<int32_t> _dynamic_array_synapses__synaptic_pre;
std::vector<double> _dynamic_array_synapses_delay;
std::vector<int32_t> _dynamic_array_synapses_N_incoming;
std::vector<int32_t> _dynamic_array_synapses_N_outgoing;
std::vector<int32_t> _dynamic_array_TC_spikemon_i;
std::vector<double> _dynamic_array_TC_spikemon_t;

//////////////// dynamic arrays 2d /////////

/////////////// static arrays /////////////
int32_t * _static_array__array_IAMPA_PYso_IN_sources;
const int _num__static_array__array_IAMPA_PYso_IN_sources = 2000;
int32_t * _static_array__array_IAMPA_PYso_IN_targets;
const int _num__static_array__array_IAMPA_PYso_IN_targets = 2000;
int32_t * _static_array__array_IAMPA_PYso_PYdr_sources;
const int _num__static_array__array_IAMPA_PYso_PYdr_sources = 2000;
int32_t * _static_array__array_IAMPA_PYso_PYdr_targets;
const int _num__static_array__array_IAMPA_PYso_PYdr_targets = 2000;
int32_t * _static_array__array_IAMPA_PYso_RE_sources;
const int _num__static_array__array_IAMPA_PYso_RE_sources = 2000;
int32_t * _static_array__array_IAMPA_PYso_RE_targets;
const int _num__static_array__array_IAMPA_PYso_RE_targets = 2000;
int32_t * _static_array__array_IAMPA_PYso_TC_sources;
const int _num__static_array__array_IAMPA_PYso_TC_sources = 2000;
int32_t * _static_array__array_IAMPA_PYso_TC_targets;
const int _num__static_array__array_IAMPA_PYso_TC_targets = 2000;
int32_t * _static_array__array_IAMPA_TC_IN_sources;
const int _num__static_array__array_IAMPA_TC_IN_sources = 400;
int32_t * _static_array__array_IAMPA_TC_IN_targets;
const int _num__static_array__array_IAMPA_TC_IN_targets = 400;
int32_t * _static_array__array_IAMPA_TC_PYdr_sources;
const int _num__static_array__array_IAMPA_TC_PYdr_sources = 420;
int32_t * _static_array__array_IAMPA_TC_PYdr_targets;
const int _num__static_array__array_IAMPA_TC_PYdr_targets = 420;
int32_t * _static_array__array_IGABAA_IN_IN_sources;
const int _num__static_array__array_IGABAA_IN_IN_sources = 380;
int32_t * _static_array__array_IGABAA_IN_IN_targets;
const int _num__static_array__array_IGABAA_IN_IN_targets = 380;
int32_t * _static_array__array_IGABAA_IN_PYso_sources;
const int _num__static_array__array_IGABAA_IN_PYso_sources = 420;
int32_t * _static_array__array_IGABAA_IN_PYso_targets;
const int _num__static_array__array_IGABAA_IN_PYso_targets = 420;
int32_t * _static_array__array_INMDA_PYso_IN_sources;
const int _num__static_array__array_INMDA_PYso_IN_sources = 2000;
int32_t * _static_array__array_INMDA_PYso_IN_targets;
const int _num__static_array__array_INMDA_PYso_IN_targets = 2000;
int32_t * _static_array__array_INMDA_PYso_PYdr_sources;
const int _num__static_array__array_INMDA_PYso_PYdr_sources = 2000;
int32_t * _static_array__array_INMDA_PYso_PYdr_targets;
const int _num__static_array__array_INMDA_PYso_PYdr_targets = 2000;
double * _static_array__array_IN_group_hNa_IN;
const int _num__static_array__array_IN_group_hNa_IN = 20;
double * _static_array__array_IN_group_nK_IN;
const int _num__static_array__array_IN_group_nK_IN = 20;
double * _static_array__array_IN_group_v;
const int _num__static_array__array_IN_group_v = 20;
double * _static_array__array_PYdr_group_CaBuffer;
const int _num__static_array__array_PYdr_group_CaBuffer = 100;
double * _static_array__array_PYdr_group_v;
const int _num__static_array__array_PYdr_group_v = 100;
double * _static_array__array_PYso_group_hA;
const int _num__static_array__array_PYso_group_hA = 100;
double * _static_array__array_PYso_group_hNa;
const int _num__static_array__array_PYso_group_hNa = 100;
double * _static_array__array_PYso_group_mKS;
const int _num__static_array__array_PYso_group_mKS = 100;
double * _static_array__array_PYso_group_nK;
const int _num__static_array__array_PYso_group_nK = 100;
double * _static_array__array_PYso_group_v;
const int _num__static_array__array_PYso_group_v = 100;
double * _static_array__array_RE_group_hNa_RE;
const int _num__static_array__array_RE_group_hNa_RE = 20;
double * _static_array__array_RE_group_hT_RE;
const int _num__static_array__array_RE_group_hT_RE = 20;
double * _static_array__array_RE_group_mNa_RE;
const int _num__static_array__array_RE_group_mNa_RE = 20;
double * _static_array__array_RE_group_mT_RE;
const int _num__static_array__array_RE_group_mT_RE = 20;
double * _static_array__array_RE_group_nK_RE;
const int _num__static_array__array_RE_group_nK_RE = 20;
double * _static_array__array_RE_group_v;
const int _num__static_array__array_RE_group_v = 20;
double * _static_array__array_TC_group_Ca_TC;
const int _num__static_array__array_TC_group_Ca_TC = 20;
double * _static_array__array_TC_group_Open;
const int _num__static_array__array_TC_group_Open = 20;
double * _static_array__array_TC_group_OpenLocked;
const int _num__static_array__array_TC_group_OpenLocked = 20;
double * _static_array__array_TC_group_Pone;
const int _num__static_array__array_TC_group_Pone = 20;
double * _static_array__array_TC_group_hNa_TC;
const int _num__static_array__array_TC_group_hNa_TC = 20;
double * _static_array__array_TC_group_hT_TC;
const int _num__static_array__array_TC_group_hT_TC = 20;
double * _static_array__array_TC_group_mNa_TC;
const int _num__static_array__array_TC_group_mNa_TC = 20;
double * _static_array__array_TC_group_nK_TC;
const int _num__static_array__array_TC_group_nK_TC = 20;
double * _static_array__array_TC_group_v;
const int _num__static_array__array_TC_group_v = 20;
double * _static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA;
const int _num__static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA = 2000;
double * _static_array__dynamic_array_IAMPA_PYso_IN_sAMPA;
const int _num__static_array__dynamic_array_IAMPA_PYso_IN_sAMPA = 2000;
double * _static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
const int _num__static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA = 2000;
double * _static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA;
const int _num__static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA = 2000;
double * _static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
const int _num__static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE = 2000;
double * _static_array__dynamic_array_IAMPA_PYso_TC_sAMPA;
const int _num__static_array__dynamic_array_IAMPA_PYso_TC_sAMPA = 2000;
double * _static_array__dynamic_array_IAMPA_TC_IN_sAMPA;
const int _num__static_array__dynamic_array_IAMPA_TC_IN_sAMPA = 400;
double * _static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA;
const int _num__static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA = 420;
double * _static_array__dynamic_array_IAMPA_TC_RE_sAMPA;
const int _num__static_array__dynamic_array_IAMPA_TC_RE_sAMPA = 400;
double * _static_array__dynamic_array_ICOM_PYdr_PYso_gCOM;
const int _num__static_array__dynamic_array_ICOM_PYdr_PYso_gCOM = 100;
double * _static_array__dynamic_array_ICOM_PYso_PYdr_gCOM;
const int _num__static_array__dynamic_array_ICOM_PYso_PYdr_gCOM = 100;
double * _static_array__dynamic_array_IGABAA_IN_IN_res_GABAA;
const int _num__static_array__dynamic_array_IGABAA_IN_IN_res_GABAA = 380;
double * _static_array__dynamic_array_IGABAA_IN_IN_sGABAA;
const int _num__static_array__dynamic_array_IGABAA_IN_IN_sGABAA = 380;
double * _static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA;
const int _num__static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA = 420;
double * _static_array__dynamic_array_IGABAA_IN_PYso_sGABAA;
const int _num__static_array__dynamic_array_IGABAA_IN_PYso_sGABAA = 420;
double * _static_array__dynamic_array_IGABAA_RE_RE_sGABAA;
const int _num__static_array__dynamic_array_IGABAA_RE_RE_sGABAA = 400;
double * _static_array__dynamic_array_IGABAA_RE_TC_sGABAA;
const int _num__static_array__dynamic_array_IGABAA_RE_TC_sGABAA = 400;
double * _static_array__dynamic_array_IGABAB_RE_TC_rGABAB;
const int _num__static_array__dynamic_array_IGABAB_RE_TC_rGABAB = 400;
double * _static_array__dynamic_array_IGABAB_RE_TC_sGABAB;
const int _num__static_array__dynamic_array_IGABAB_RE_TC_sGABAB = 400;
double * _static_array__dynamic_array_IKNa_PYso_PYdr_concNa;
const int _num__static_array__dynamic_array_IKNa_PYso_PYdr_concNa = 100;
double * _static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local;
const int _num__static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local = 100;
double * _static_array__dynamic_array_INMDA_PYso_IN_res_NMDA;
const int _num__static_array__dynamic_array_INMDA_PYso_IN_res_NMDA = 2000;
double * _static_array__dynamic_array_INMDA_PYso_IN_sNMDA;
const int _num__static_array__dynamic_array_INMDA_PYso_IN_sNMDA = 2000;
double * _static_array__dynamic_array_INMDA_PYso_IN_xNMDA;
const int _num__static_array__dynamic_array_INMDA_PYso_IN_xNMDA = 2000;
double * _static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA;
const int _num__static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA = 2000;
double * _static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA;
const int _num__static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA = 2000;
double * _static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA;
const int _num__static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA = 2000;

//////////////// synapses /////////////////
// IAMPA_PYso_IN
SynapticPathway IAMPA_PYso_IN_down(
    _dynamic_array_IAMPA_PYso_IN__synaptic_pre,
    0, 100);
SynapticPathway IAMPA_PYso_IN_up(
    _dynamic_array_IAMPA_PYso_IN__synaptic_pre,
    0, 100);
// IAMPA_PYso_PYdr
SynapticPathway IAMPA_PYso_PYdr_down(
    _dynamic_array_IAMPA_PYso_PYdr__synaptic_pre,
    0, 100);
SynapticPathway IAMPA_PYso_PYdr_up(
    _dynamic_array_IAMPA_PYso_PYdr__synaptic_pre,
    0, 100);
// IAMPA_PYso_RE
// IAMPA_PYso_TC
// IAMPA_TC_IN
// IAMPA_TC_PYdr
// IAMPA_TC_RE
// ICOM_PYdr_PYso
// ICOM_PYso_PYdr
// IGABAA_IN_IN
SynapticPathway IGABAA_IN_IN_down(
    _dynamic_array_IGABAA_IN_IN__synaptic_pre,
    0, 20);
SynapticPathway IGABAA_IN_IN_up(
    _dynamic_array_IGABAA_IN_IN__synaptic_pre,
    0, 20);
// IGABAA_IN_PYso
SynapticPathway IGABAA_IN_PYso_down(
    _dynamic_array_IGABAA_IN_PYso__synaptic_pre,
    0, 20);
SynapticPathway IGABAA_IN_PYso_up(
    _dynamic_array_IGABAA_IN_PYso__synaptic_pre,
    0, 20);
// IGABAA_RE_RE
// IGABAA_RE_TC
// IGABAB_RE_TC
// IKNa_PYso_PYdr
// INMDA_PYso_IN
SynapticPathway INMDA_PYso_IN_down(
    _dynamic_array_INMDA_PYso_IN__synaptic_pre,
    0, 100);
SynapticPathway INMDA_PYso_IN_up(
    _dynamic_array_INMDA_PYso_IN__synaptic_pre,
    0, 100);
// INMDA_PYso_PYdr
SynapticPathway INMDA_PYso_PYdr_down(
    _dynamic_array_INMDA_PYso_PYdr__synaptic_pre,
    0, 100);
SynapticPathway INMDA_PYso_PYdr_up(
    _dynamic_array_INMDA_PYso_PYdr__synaptic_pre,
    0, 100);
// synapses
SynapticPathway synapses_pre(
    _dynamic_array_synapses__synaptic_pre,
    0, 100);
// synapses_1
SynapticPathway synapses_1_pre(
    _dynamic_array_synapses_1__synaptic_pre,
    0, 20);
// synapses_2
SynapticPathway synapses_2_pre(
    _dynamic_array_synapses_2__synaptic_pre,
    0, 20);

//////////////// clocks ///////////////////
// attributes will be set in run.cpp
Clock defaultclock;

// Profiling information for each code object
}

void _init_arrays()
{
    using namespace brian;

    // Arrays initialized to 0
    _array_clock_dt = new double[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_clock_dt[i] = 0;

    _array_clock_t = new double[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_clock_t[i] = 0;

    _array_clock_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_clock_timestep[i] = 0;

    _array_defaultclock_dt = new double[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_defaultclock_dt[i] = 0;

    _array_defaultclock_t = new double[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_defaultclock_t[i] = 0;

    _array_defaultclock_timestep = new int64_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_defaultclock_timestep[i] = 0;

    _array_IAMPA_PYso_IN_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IAMPA_PYso_IN_N[i] = 0;

    _array_IAMPA_PYso_IN_sources = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_IN_sources[i] = 0;

    _array_IAMPA_PYso_IN_targets = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_IN_targets[i] = 0;

    _array_IAMPA_PYso_PYdr_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IAMPA_PYso_PYdr_N[i] = 0;

    _array_IAMPA_PYso_PYdr_sources = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_PYdr_sources[i] = 0;

    _array_IAMPA_PYso_PYdr_targets = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_PYdr_targets[i] = 0;

    _array_IAMPA_PYso_RE_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IAMPA_PYso_RE_N[i] = 0;

    _array_IAMPA_PYso_RE_sources = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_RE_sources[i] = 0;

    _array_IAMPA_PYso_RE_targets = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_RE_targets[i] = 0;

    _array_IAMPA_PYso_TC_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IAMPA_PYso_TC_N[i] = 0;

    _array_IAMPA_PYso_TC_sources = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_TC_sources[i] = 0;

    _array_IAMPA_PYso_TC_targets = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_IAMPA_PYso_TC_targets[i] = 0;

    _array_IAMPA_TC_IN_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IAMPA_TC_IN_N[i] = 0;

    _array_IAMPA_TC_IN_sources = new int32_t[400];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<400; i++) _array_IAMPA_TC_IN_sources[i] = 0;

    _array_IAMPA_TC_IN_targets = new int32_t[400];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<400; i++) _array_IAMPA_TC_IN_targets[i] = 0;

    _array_IAMPA_TC_PYdr_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IAMPA_TC_PYdr_N[i] = 0;

    _array_IAMPA_TC_PYdr_sources = new int32_t[420];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<420; i++) _array_IAMPA_TC_PYdr_sources[i] = 0;

    _array_IAMPA_TC_PYdr_targets = new int32_t[420];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<420; i++) _array_IAMPA_TC_PYdr_targets[i] = 0;

    _array_IAMPA_TC_RE_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IAMPA_TC_RE_N[i] = 0;

    _array_ICOM_PYdr_PYso_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_ICOM_PYdr_PYso_N[i] = 0;

    _array_ICOM_PYso_PYdr_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_ICOM_PYso_PYdr_N[i] = 0;

    _array_IGABAA_IN_IN_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IGABAA_IN_IN_N[i] = 0;

    _array_IGABAA_IN_IN_sources = new int32_t[380];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<380; i++) _array_IGABAA_IN_IN_sources[i] = 0;

    _array_IGABAA_IN_IN_targets = new int32_t[380];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<380; i++) _array_IGABAA_IN_IN_targets[i] = 0;

    _array_IGABAA_IN_PYso_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IGABAA_IN_PYso_N[i] = 0;

    _array_IGABAA_IN_PYso_sources = new int32_t[420];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<420; i++) _array_IGABAA_IN_PYso_sources[i] = 0;

    _array_IGABAA_IN_PYso_targets = new int32_t[420];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<420; i++) _array_IGABAA_IN_PYso_targets[i] = 0;

    _array_IGABAA_RE_RE_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IGABAA_RE_RE_N[i] = 0;

    _array_IGABAA_RE_TC_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IGABAA_RE_TC_N[i] = 0;

    _array_IGABAB_RE_TC_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IGABAB_RE_TC_N[i] = 0;

    _array_IKNa_PYso_PYdr_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IKNa_PYso_PYdr_N[i] = 0;

    _array_IN_group__spikespace = new int32_t[21];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<21; i++) _array_IN_group__spikespace[i] = 0;

    _array_IN_group_hNa_IN = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_hNa_IN[i] = 0;

    _array_IN_group_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_i[i] = 0;

    _array_IN_group_iAMPA_PYso_IN = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_iAMPA_PYso_IN[i] = 0;

    _array_IN_group_iAMPA_TC_IN = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_iAMPA_TC_IN[i] = 0;

    _array_IN_group_iGABAA_IN_IN = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_iGABAA_IN_IN[i] = 0;

    _array_IN_group_iNMDA_PYso_IN = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_iNMDA_PYso_IN[i] = 0;

    _array_IN_group_lastspike = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_lastspike[i] = 0;

    _array_IN_group_nK_IN = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_nK_IN[i] = 0;

    _array_IN_group_not_refractory = new char[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_not_refractory[i] = 0;

    _array_IN_group_v = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_v[i] = 0;

    _array_IN_spikemon__source_idx = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_spikemon__source_idx[i] = 0;

    _array_IN_spikemon_count = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_spikemon_count[i] = 0;

    _array_IN_spikemon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_IN_spikemon_N[i] = 0;

    _array_INMDA_PYso_IN_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_INMDA_PYso_IN_N[i] = 0;

    _array_INMDA_PYso_IN_sources = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_INMDA_PYso_IN_sources[i] = 0;

    _array_INMDA_PYso_IN_targets = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_INMDA_PYso_IN_targets[i] = 0;

    _array_INMDA_PYso_PYdr_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_INMDA_PYso_PYdr_N[i] = 0;

    _array_INMDA_PYso_PYdr_sources = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_INMDA_PYso_PYdr_sources[i] = 0;

    _array_INMDA_PYso_PYdr_targets = new int32_t[2000];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<2000; i++) _array_INMDA_PYso_PYdr_targets[i] = 0;

    _array_poissongroup_1__spikespace = new int32_t[21];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<21; i++) _array_poissongroup_1__spikespace[i] = 0;

    _array_poissongroup_1_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_poissongroup_1_i[i] = 0;

    _array_poissongroup_1_rates = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_poissongroup_1_rates[i] = 0;

    _array_poissongroup_2__spikespace = new int32_t[21];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<21; i++) _array_poissongroup_2__spikespace[i] = 0;

    _array_poissongroup_2_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_poissongroup_2_i[i] = 0;

    _array_poissongroup_2_rates = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_poissongroup_2_rates[i] = 0;

    _array_poissongroup__spikespace = new int32_t[101];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<101; i++) _array_poissongroup__spikespace[i] = 0;

    _array_poissongroup_i = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_poissongroup_i[i] = 0;

    _array_poissongroup_rates = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_poissongroup_rates[i] = 0;

    _array_PYdr_group__spikespace = new int32_t[101];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<101; i++) _array_PYdr_group__spikespace[i] = 0;

    _array_PYdr_group_CaBuffer = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_CaBuffer[i] = 0;

    _array_PYdr_group_gPoisson_PYdr = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_gPoisson_PYdr[i] = 0;

    _array_PYdr_group_i = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_i[i] = 0;

    _array_PYdr_group_iAMPA_PYso_PYdr = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_iAMPA_PYso_PYdr[i] = 0;

    _array_PYdr_group_iAMPA_TC_PYdr = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_iAMPA_TC_PYdr[i] = 0;

    _array_PYdr_group_iCOM = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_iCOM[i] = 0;

    _array_PYdr_group_iNMDA_PYso_PYdr = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_iNMDA_PYso_PYdr[i] = 0;

    _array_PYdr_group_lastspike = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_lastspike[i] = 0;

    _array_PYdr_group_not_refractory = new char[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_not_refractory[i] = 0;

    _array_PYdr_group_v = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_v[i] = 0;

    _array_PYdr_spikemon__source_idx = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_spikemon__source_idx[i] = 0;

    _array_PYdr_spikemon_count = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_spikemon_count[i] = 0;

    _array_PYdr_spikemon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_PYdr_spikemon_N[i] = 0;

    _array_PYso_group__spikespace = new int32_t[101];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<101; i++) _array_PYso_group__spikespace[i] = 0;

    _array_PYso_group_hA = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_hA[i] = 0;

    _array_PYso_group_hNa = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_hNa[i] = 0;

    _array_PYso_group_i = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_i[i] = 0;

    _array_PYso_group_iCOM = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_iCOM[i] = 0;

    _array_PYso_group_iGABAA_IN_PYso = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_iGABAA_IN_PYso[i] = 0;

    _array_PYso_group_iKNa = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_iKNa[i] = 0;

    _array_PYso_group_lastspike = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_lastspike[i] = 0;

    _array_PYso_group_mKS = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_mKS[i] = 0;

    _array_PYso_group_nK = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_nK[i] = 0;

    _array_PYso_group_not_refractory = new char[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_not_refractory[i] = 0;

    _array_PYso_group_v = new double[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_v[i] = 0;

    _array_PYso_spikemon__source_idx = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_spikemon__source_idx[i] = 0;

    _array_PYso_spikemon_count = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_spikemon_count[i] = 0;

    _array_PYso_spikemon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_PYso_spikemon_N[i] = 0;

    _array_RE_group__spikespace = new int32_t[21];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<21; i++) _array_RE_group__spikespace[i] = 0;

    _array_RE_group_gPoisson_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_gPoisson_RE[i] = 0;

    _array_RE_group_hNa_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_hNa_RE[i] = 0;

    _array_RE_group_hT_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_hT_RE[i] = 0;

    _array_RE_group_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_i[i] = 0;

    _array_RE_group_iAMPA_PYso_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_iAMPA_PYso_RE[i] = 0;

    _array_RE_group_iAMPA_TC_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_iAMPA_TC_RE[i] = 0;

    _array_RE_group_iGABAA_RE_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_iGABAA_RE_RE[i] = 0;

    _array_RE_group_lastspike = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_lastspike[i] = 0;

    _array_RE_group_mNa_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_mNa_RE[i] = 0;

    _array_RE_group_mT_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_mT_RE[i] = 0;

    _array_RE_group_nK_RE = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_nK_RE[i] = 0;

    _array_RE_group_not_refractory = new char[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_not_refractory[i] = 0;

    _array_RE_group_v = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_v[i] = 0;

    _array_RE_spikemon__source_idx = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_spikemon__source_idx[i] = 0;

    _array_RE_spikemon_count = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_spikemon_count[i] = 0;

    _array_RE_spikemon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_RE_spikemon_N[i] = 0;

    _array_synapses_1_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_synapses_1_N[i] = 0;

    _array_synapses_2_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_synapses_2_N[i] = 0;

    _array_synapses_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_synapses_N[i] = 0;

    _array_TC_group__spikespace = new int32_t[21];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<21; i++) _array_TC_group__spikespace[i] = 0;

    _array_TC_group_Ca_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_Ca_TC[i] = 0;

    _array_TC_group_dCa_extra = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_dCa_extra[i] = 0;

    _array_TC_group_gPoisson_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_gPoisson_TC[i] = 0;

    _array_TC_group_hNa_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_hNa_TC[i] = 0;

    _array_TC_group_hT_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_hT_TC[i] = 0;

    _array_TC_group_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_i[i] = 0;

    _array_TC_group_iAMPA_PYso_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_iAMPA_PYso_TC[i] = 0;

    _array_TC_group_iGABAA_RE_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_iGABAA_RE_TC[i] = 0;

    _array_TC_group_iGABAB_RE_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_iGABAB_RE_TC[i] = 0;

    _array_TC_group_lastspike = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_lastspike[i] = 0;

    _array_TC_group_mNa_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_mNa_TC[i] = 0;

    _array_TC_group_nK_TC = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_nK_TC[i] = 0;

    _array_TC_group_not_refractory = new char[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_not_refractory[i] = 0;

    _array_TC_group_Open = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_Open[i] = 0;

    _array_TC_group_OpenLocked = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_OpenLocked[i] = 0;

    _array_TC_group_Pone = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_Pone[i] = 0;

    _array_TC_group_v = new double[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_v[i] = 0;

    _array_TC_spikemon__source_idx = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_spikemon__source_idx[i] = 0;

    _array_TC_spikemon_count = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_spikemon_count[i] = 0;

    _array_TC_spikemon_N = new int32_t[1];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _array_TC_spikemon_N[i] = 0;

    _dynamic_array_IAMPA_PYso_IN_delay.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IAMPA_PYso_IN_delay[i] = 0;

    _dynamic_array_IAMPA_PYso_IN_delay_1.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IAMPA_PYso_IN_delay_1[i] = 0;

    _dynamic_array_IAMPA_PYso_PYdr_delay.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IAMPA_PYso_PYdr_delay[i] = 0;

    _dynamic_array_IAMPA_PYso_PYdr_delay_1.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IAMPA_PYso_PYdr_delay_1[i] = 0;

    _dynamic_array_IGABAA_IN_IN_delay.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IGABAA_IN_IN_delay[i] = 0;

    _dynamic_array_IGABAA_IN_IN_delay_1.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IGABAA_IN_IN_delay_1[i] = 0;

    _dynamic_array_IGABAA_IN_PYso_delay.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IGABAA_IN_PYso_delay[i] = 0;

    _dynamic_array_IGABAA_IN_PYso_delay_1.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_IGABAA_IN_PYso_delay_1[i] = 0;

    _dynamic_array_INMDA_PYso_IN_delay.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_INMDA_PYso_IN_delay[i] = 0;

    _dynamic_array_INMDA_PYso_IN_delay_1.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_INMDA_PYso_IN_delay_1[i] = 0;

    _dynamic_array_INMDA_PYso_PYdr_delay.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_INMDA_PYso_PYdr_delay[i] = 0;

    _dynamic_array_INMDA_PYso_PYdr_delay_1.resize(1);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<1; i++) _dynamic_array_INMDA_PYso_PYdr_delay_1[i] = 0;


    // Arrays initialized to an "arange"
    _array_IN_group_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_group_i[i] = 0 + i;

    _array_IN_spikemon__source_idx = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_IN_spikemon__source_idx[i] = 0 + i;

    _array_poissongroup_1_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_poissongroup_1_i[i] = 0 + i;

    _array_poissongroup_2_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_poissongroup_2_i[i] = 0 + i;

    _array_poissongroup_i = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_poissongroup_i[i] = 0 + i;

    _array_PYdr_group_i = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_group_i[i] = 0 + i;

    _array_PYdr_spikemon__source_idx = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYdr_spikemon__source_idx[i] = 0 + i;

    _array_PYso_group_i = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_group_i[i] = 0 + i;

    _array_PYso_spikemon__source_idx = new int32_t[100];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<100; i++) _array_PYso_spikemon__source_idx[i] = 0 + i;

    _array_RE_group_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_group_i[i] = 0 + i;

    _array_RE_spikemon__source_idx = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_RE_spikemon__source_idx[i] = 0 + i;

    _array_TC_group_i = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_group_i[i] = 0 + i;

    _array_TC_spikemon__source_idx = new int32_t[20];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<20; i++) _array_TC_spikemon__source_idx[i] = 0 + i;


    // static arrays
    _static_array__array_IAMPA_PYso_IN_sources = new int32_t[2000];
    _static_array__array_IAMPA_PYso_IN_targets = new int32_t[2000];
    _static_array__array_IAMPA_PYso_PYdr_sources = new int32_t[2000];
    _static_array__array_IAMPA_PYso_PYdr_targets = new int32_t[2000];
    _static_array__array_IAMPA_PYso_RE_sources = new int32_t[2000];
    _static_array__array_IAMPA_PYso_RE_targets = new int32_t[2000];
    _static_array__array_IAMPA_PYso_TC_sources = new int32_t[2000];
    _static_array__array_IAMPA_PYso_TC_targets = new int32_t[2000];
    _static_array__array_IAMPA_TC_IN_sources = new int32_t[400];
    _static_array__array_IAMPA_TC_IN_targets = new int32_t[400];
    _static_array__array_IAMPA_TC_PYdr_sources = new int32_t[420];
    _static_array__array_IAMPA_TC_PYdr_targets = new int32_t[420];
    _static_array__array_IGABAA_IN_IN_sources = new int32_t[380];
    _static_array__array_IGABAA_IN_IN_targets = new int32_t[380];
    _static_array__array_IGABAA_IN_PYso_sources = new int32_t[420];
    _static_array__array_IGABAA_IN_PYso_targets = new int32_t[420];
    _static_array__array_INMDA_PYso_IN_sources = new int32_t[2000];
    _static_array__array_INMDA_PYso_IN_targets = new int32_t[2000];
    _static_array__array_INMDA_PYso_PYdr_sources = new int32_t[2000];
    _static_array__array_INMDA_PYso_PYdr_targets = new int32_t[2000];
    _static_array__array_IN_group_hNa_IN = new double[20];
    _static_array__array_IN_group_nK_IN = new double[20];
    _static_array__array_IN_group_v = new double[20];
    _static_array__array_PYdr_group_CaBuffer = new double[100];
    _static_array__array_PYdr_group_v = new double[100];
    _static_array__array_PYso_group_hA = new double[100];
    _static_array__array_PYso_group_hNa = new double[100];
    _static_array__array_PYso_group_mKS = new double[100];
    _static_array__array_PYso_group_nK = new double[100];
    _static_array__array_PYso_group_v = new double[100];
    _static_array__array_RE_group_hNa_RE = new double[20];
    _static_array__array_RE_group_hT_RE = new double[20];
    _static_array__array_RE_group_mNa_RE = new double[20];
    _static_array__array_RE_group_mT_RE = new double[20];
    _static_array__array_RE_group_nK_RE = new double[20];
    _static_array__array_RE_group_v = new double[20];
    _static_array__array_TC_group_Ca_TC = new double[20];
    _static_array__array_TC_group_Open = new double[20];
    _static_array__array_TC_group_OpenLocked = new double[20];
    _static_array__array_TC_group_Pone = new double[20];
    _static_array__array_TC_group_hNa_TC = new double[20];
    _static_array__array_TC_group_hT_TC = new double[20];
    _static_array__array_TC_group_mNa_TC = new double[20];
    _static_array__array_TC_group_nK_TC = new double[20];
    _static_array__array_TC_group_v = new double[20];
    _static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA = new double[2000];
    _static_array__dynamic_array_IAMPA_PYso_IN_sAMPA = new double[2000];
    _static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA = new double[2000];
    _static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA = new double[2000];
    _static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE = new double[2000];
    _static_array__dynamic_array_IAMPA_PYso_TC_sAMPA = new double[2000];
    _static_array__dynamic_array_IAMPA_TC_IN_sAMPA = new double[400];
    _static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA = new double[420];
    _static_array__dynamic_array_IAMPA_TC_RE_sAMPA = new double[400];
    _static_array__dynamic_array_ICOM_PYdr_PYso_gCOM = new double[100];
    _static_array__dynamic_array_ICOM_PYso_PYdr_gCOM = new double[100];
    _static_array__dynamic_array_IGABAA_IN_IN_res_GABAA = new double[380];
    _static_array__dynamic_array_IGABAA_IN_IN_sGABAA = new double[380];
    _static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA = new double[420];
    _static_array__dynamic_array_IGABAA_IN_PYso_sGABAA = new double[420];
    _static_array__dynamic_array_IGABAA_RE_RE_sGABAA = new double[400];
    _static_array__dynamic_array_IGABAA_RE_TC_sGABAA = new double[400];
    _static_array__dynamic_array_IGABAB_RE_TC_rGABAB = new double[400];
    _static_array__dynamic_array_IGABAB_RE_TC_sGABAB = new double[400];
    _static_array__dynamic_array_IKNa_PYso_PYdr_concNa = new double[100];
    _static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local = new double[100];
    _static_array__dynamic_array_INMDA_PYso_IN_res_NMDA = new double[2000];
    _static_array__dynamic_array_INMDA_PYso_IN_sNMDA = new double[2000];
    _static_array__dynamic_array_INMDA_PYso_IN_xNMDA = new double[2000];
    _static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA = new double[2000];
    _static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA = new double[2000];
    _static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA = new double[2000];

    // Random number generator states
    std::random_device rd;
    for (int i=0; i<2; i++)
        _random_generators.push_back(RandomGenerator());
}

void _load_arrays()
{
    using namespace brian;

    ifstream f_static_array__array_IAMPA_PYso_IN_sources;
    f_static_array__array_IAMPA_PYso_IN_sources.open("static_arrays/_static_array__array_IAMPA_PYso_IN_sources", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_IN_sources.is_open())
    {
        f_static_array__array_IAMPA_PYso_IN_sources.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_IN_sources), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_IN_sources." << endl;
    }
    ifstream f_static_array__array_IAMPA_PYso_IN_targets;
    f_static_array__array_IAMPA_PYso_IN_targets.open("static_arrays/_static_array__array_IAMPA_PYso_IN_targets", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_IN_targets.is_open())
    {
        f_static_array__array_IAMPA_PYso_IN_targets.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_IN_targets), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_IN_targets." << endl;
    }
    ifstream f_static_array__array_IAMPA_PYso_PYdr_sources;
    f_static_array__array_IAMPA_PYso_PYdr_sources.open("static_arrays/_static_array__array_IAMPA_PYso_PYdr_sources", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_PYdr_sources.is_open())
    {
        f_static_array__array_IAMPA_PYso_PYdr_sources.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_PYdr_sources), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_PYdr_sources." << endl;
    }
    ifstream f_static_array__array_IAMPA_PYso_PYdr_targets;
    f_static_array__array_IAMPA_PYso_PYdr_targets.open("static_arrays/_static_array__array_IAMPA_PYso_PYdr_targets", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_PYdr_targets.is_open())
    {
        f_static_array__array_IAMPA_PYso_PYdr_targets.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_PYdr_targets), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_PYdr_targets." << endl;
    }
    ifstream f_static_array__array_IAMPA_PYso_RE_sources;
    f_static_array__array_IAMPA_PYso_RE_sources.open("static_arrays/_static_array__array_IAMPA_PYso_RE_sources", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_RE_sources.is_open())
    {
        f_static_array__array_IAMPA_PYso_RE_sources.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_RE_sources), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_RE_sources." << endl;
    }
    ifstream f_static_array__array_IAMPA_PYso_RE_targets;
    f_static_array__array_IAMPA_PYso_RE_targets.open("static_arrays/_static_array__array_IAMPA_PYso_RE_targets", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_RE_targets.is_open())
    {
        f_static_array__array_IAMPA_PYso_RE_targets.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_RE_targets), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_RE_targets." << endl;
    }
    ifstream f_static_array__array_IAMPA_PYso_TC_sources;
    f_static_array__array_IAMPA_PYso_TC_sources.open("static_arrays/_static_array__array_IAMPA_PYso_TC_sources", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_TC_sources.is_open())
    {
        f_static_array__array_IAMPA_PYso_TC_sources.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_TC_sources), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_TC_sources." << endl;
    }
    ifstream f_static_array__array_IAMPA_PYso_TC_targets;
    f_static_array__array_IAMPA_PYso_TC_targets.open("static_arrays/_static_array__array_IAMPA_PYso_TC_targets", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_PYso_TC_targets.is_open())
    {
        f_static_array__array_IAMPA_PYso_TC_targets.read(reinterpret_cast<char*>(_static_array__array_IAMPA_PYso_TC_targets), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_PYso_TC_targets." << endl;
    }
    ifstream f_static_array__array_IAMPA_TC_IN_sources;
    f_static_array__array_IAMPA_TC_IN_sources.open("static_arrays/_static_array__array_IAMPA_TC_IN_sources", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_TC_IN_sources.is_open())
    {
        f_static_array__array_IAMPA_TC_IN_sources.read(reinterpret_cast<char*>(_static_array__array_IAMPA_TC_IN_sources), 400*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_TC_IN_sources." << endl;
    }
    ifstream f_static_array__array_IAMPA_TC_IN_targets;
    f_static_array__array_IAMPA_TC_IN_targets.open("static_arrays/_static_array__array_IAMPA_TC_IN_targets", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_TC_IN_targets.is_open())
    {
        f_static_array__array_IAMPA_TC_IN_targets.read(reinterpret_cast<char*>(_static_array__array_IAMPA_TC_IN_targets), 400*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_TC_IN_targets." << endl;
    }
    ifstream f_static_array__array_IAMPA_TC_PYdr_sources;
    f_static_array__array_IAMPA_TC_PYdr_sources.open("static_arrays/_static_array__array_IAMPA_TC_PYdr_sources", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_TC_PYdr_sources.is_open())
    {
        f_static_array__array_IAMPA_TC_PYdr_sources.read(reinterpret_cast<char*>(_static_array__array_IAMPA_TC_PYdr_sources), 420*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_TC_PYdr_sources." << endl;
    }
    ifstream f_static_array__array_IAMPA_TC_PYdr_targets;
    f_static_array__array_IAMPA_TC_PYdr_targets.open("static_arrays/_static_array__array_IAMPA_TC_PYdr_targets", ios::in | ios::binary);
    if(f_static_array__array_IAMPA_TC_PYdr_targets.is_open())
    {
        f_static_array__array_IAMPA_TC_PYdr_targets.read(reinterpret_cast<char*>(_static_array__array_IAMPA_TC_PYdr_targets), 420*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IAMPA_TC_PYdr_targets." << endl;
    }
    ifstream f_static_array__array_IGABAA_IN_IN_sources;
    f_static_array__array_IGABAA_IN_IN_sources.open("static_arrays/_static_array__array_IGABAA_IN_IN_sources", ios::in | ios::binary);
    if(f_static_array__array_IGABAA_IN_IN_sources.is_open())
    {
        f_static_array__array_IGABAA_IN_IN_sources.read(reinterpret_cast<char*>(_static_array__array_IGABAA_IN_IN_sources), 380*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IGABAA_IN_IN_sources." << endl;
    }
    ifstream f_static_array__array_IGABAA_IN_IN_targets;
    f_static_array__array_IGABAA_IN_IN_targets.open("static_arrays/_static_array__array_IGABAA_IN_IN_targets", ios::in | ios::binary);
    if(f_static_array__array_IGABAA_IN_IN_targets.is_open())
    {
        f_static_array__array_IGABAA_IN_IN_targets.read(reinterpret_cast<char*>(_static_array__array_IGABAA_IN_IN_targets), 380*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IGABAA_IN_IN_targets." << endl;
    }
    ifstream f_static_array__array_IGABAA_IN_PYso_sources;
    f_static_array__array_IGABAA_IN_PYso_sources.open("static_arrays/_static_array__array_IGABAA_IN_PYso_sources", ios::in | ios::binary);
    if(f_static_array__array_IGABAA_IN_PYso_sources.is_open())
    {
        f_static_array__array_IGABAA_IN_PYso_sources.read(reinterpret_cast<char*>(_static_array__array_IGABAA_IN_PYso_sources), 420*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IGABAA_IN_PYso_sources." << endl;
    }
    ifstream f_static_array__array_IGABAA_IN_PYso_targets;
    f_static_array__array_IGABAA_IN_PYso_targets.open("static_arrays/_static_array__array_IGABAA_IN_PYso_targets", ios::in | ios::binary);
    if(f_static_array__array_IGABAA_IN_PYso_targets.is_open())
    {
        f_static_array__array_IGABAA_IN_PYso_targets.read(reinterpret_cast<char*>(_static_array__array_IGABAA_IN_PYso_targets), 420*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IGABAA_IN_PYso_targets." << endl;
    }
    ifstream f_static_array__array_INMDA_PYso_IN_sources;
    f_static_array__array_INMDA_PYso_IN_sources.open("static_arrays/_static_array__array_INMDA_PYso_IN_sources", ios::in | ios::binary);
    if(f_static_array__array_INMDA_PYso_IN_sources.is_open())
    {
        f_static_array__array_INMDA_PYso_IN_sources.read(reinterpret_cast<char*>(_static_array__array_INMDA_PYso_IN_sources), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_INMDA_PYso_IN_sources." << endl;
    }
    ifstream f_static_array__array_INMDA_PYso_IN_targets;
    f_static_array__array_INMDA_PYso_IN_targets.open("static_arrays/_static_array__array_INMDA_PYso_IN_targets", ios::in | ios::binary);
    if(f_static_array__array_INMDA_PYso_IN_targets.is_open())
    {
        f_static_array__array_INMDA_PYso_IN_targets.read(reinterpret_cast<char*>(_static_array__array_INMDA_PYso_IN_targets), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_INMDA_PYso_IN_targets." << endl;
    }
    ifstream f_static_array__array_INMDA_PYso_PYdr_sources;
    f_static_array__array_INMDA_PYso_PYdr_sources.open("static_arrays/_static_array__array_INMDA_PYso_PYdr_sources", ios::in | ios::binary);
    if(f_static_array__array_INMDA_PYso_PYdr_sources.is_open())
    {
        f_static_array__array_INMDA_PYso_PYdr_sources.read(reinterpret_cast<char*>(_static_array__array_INMDA_PYso_PYdr_sources), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_INMDA_PYso_PYdr_sources." << endl;
    }
    ifstream f_static_array__array_INMDA_PYso_PYdr_targets;
    f_static_array__array_INMDA_PYso_PYdr_targets.open("static_arrays/_static_array__array_INMDA_PYso_PYdr_targets", ios::in | ios::binary);
    if(f_static_array__array_INMDA_PYso_PYdr_targets.is_open())
    {
        f_static_array__array_INMDA_PYso_PYdr_targets.read(reinterpret_cast<char*>(_static_array__array_INMDA_PYso_PYdr_targets), 2000*sizeof(int32_t));
    } else
    {
        std::cout << "Error opening static array _static_array__array_INMDA_PYso_PYdr_targets." << endl;
    }
    ifstream f_static_array__array_IN_group_hNa_IN;
    f_static_array__array_IN_group_hNa_IN.open("static_arrays/_static_array__array_IN_group_hNa_IN", ios::in | ios::binary);
    if(f_static_array__array_IN_group_hNa_IN.is_open())
    {
        f_static_array__array_IN_group_hNa_IN.read(reinterpret_cast<char*>(_static_array__array_IN_group_hNa_IN), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IN_group_hNa_IN." << endl;
    }
    ifstream f_static_array__array_IN_group_nK_IN;
    f_static_array__array_IN_group_nK_IN.open("static_arrays/_static_array__array_IN_group_nK_IN", ios::in | ios::binary);
    if(f_static_array__array_IN_group_nK_IN.is_open())
    {
        f_static_array__array_IN_group_nK_IN.read(reinterpret_cast<char*>(_static_array__array_IN_group_nK_IN), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IN_group_nK_IN." << endl;
    }
    ifstream f_static_array__array_IN_group_v;
    f_static_array__array_IN_group_v.open("static_arrays/_static_array__array_IN_group_v", ios::in | ios::binary);
    if(f_static_array__array_IN_group_v.is_open())
    {
        f_static_array__array_IN_group_v.read(reinterpret_cast<char*>(_static_array__array_IN_group_v), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_IN_group_v." << endl;
    }
    ifstream f_static_array__array_PYdr_group_CaBuffer;
    f_static_array__array_PYdr_group_CaBuffer.open("static_arrays/_static_array__array_PYdr_group_CaBuffer", ios::in | ios::binary);
    if(f_static_array__array_PYdr_group_CaBuffer.is_open())
    {
        f_static_array__array_PYdr_group_CaBuffer.read(reinterpret_cast<char*>(_static_array__array_PYdr_group_CaBuffer), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_PYdr_group_CaBuffer." << endl;
    }
    ifstream f_static_array__array_PYdr_group_v;
    f_static_array__array_PYdr_group_v.open("static_arrays/_static_array__array_PYdr_group_v", ios::in | ios::binary);
    if(f_static_array__array_PYdr_group_v.is_open())
    {
        f_static_array__array_PYdr_group_v.read(reinterpret_cast<char*>(_static_array__array_PYdr_group_v), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_PYdr_group_v." << endl;
    }
    ifstream f_static_array__array_PYso_group_hA;
    f_static_array__array_PYso_group_hA.open("static_arrays/_static_array__array_PYso_group_hA", ios::in | ios::binary);
    if(f_static_array__array_PYso_group_hA.is_open())
    {
        f_static_array__array_PYso_group_hA.read(reinterpret_cast<char*>(_static_array__array_PYso_group_hA), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_PYso_group_hA." << endl;
    }
    ifstream f_static_array__array_PYso_group_hNa;
    f_static_array__array_PYso_group_hNa.open("static_arrays/_static_array__array_PYso_group_hNa", ios::in | ios::binary);
    if(f_static_array__array_PYso_group_hNa.is_open())
    {
        f_static_array__array_PYso_group_hNa.read(reinterpret_cast<char*>(_static_array__array_PYso_group_hNa), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_PYso_group_hNa." << endl;
    }
    ifstream f_static_array__array_PYso_group_mKS;
    f_static_array__array_PYso_group_mKS.open("static_arrays/_static_array__array_PYso_group_mKS", ios::in | ios::binary);
    if(f_static_array__array_PYso_group_mKS.is_open())
    {
        f_static_array__array_PYso_group_mKS.read(reinterpret_cast<char*>(_static_array__array_PYso_group_mKS), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_PYso_group_mKS." << endl;
    }
    ifstream f_static_array__array_PYso_group_nK;
    f_static_array__array_PYso_group_nK.open("static_arrays/_static_array__array_PYso_group_nK", ios::in | ios::binary);
    if(f_static_array__array_PYso_group_nK.is_open())
    {
        f_static_array__array_PYso_group_nK.read(reinterpret_cast<char*>(_static_array__array_PYso_group_nK), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_PYso_group_nK." << endl;
    }
    ifstream f_static_array__array_PYso_group_v;
    f_static_array__array_PYso_group_v.open("static_arrays/_static_array__array_PYso_group_v", ios::in | ios::binary);
    if(f_static_array__array_PYso_group_v.is_open())
    {
        f_static_array__array_PYso_group_v.read(reinterpret_cast<char*>(_static_array__array_PYso_group_v), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_PYso_group_v." << endl;
    }
    ifstream f_static_array__array_RE_group_hNa_RE;
    f_static_array__array_RE_group_hNa_RE.open("static_arrays/_static_array__array_RE_group_hNa_RE", ios::in | ios::binary);
    if(f_static_array__array_RE_group_hNa_RE.is_open())
    {
        f_static_array__array_RE_group_hNa_RE.read(reinterpret_cast<char*>(_static_array__array_RE_group_hNa_RE), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_RE_group_hNa_RE." << endl;
    }
    ifstream f_static_array__array_RE_group_hT_RE;
    f_static_array__array_RE_group_hT_RE.open("static_arrays/_static_array__array_RE_group_hT_RE", ios::in | ios::binary);
    if(f_static_array__array_RE_group_hT_RE.is_open())
    {
        f_static_array__array_RE_group_hT_RE.read(reinterpret_cast<char*>(_static_array__array_RE_group_hT_RE), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_RE_group_hT_RE." << endl;
    }
    ifstream f_static_array__array_RE_group_mNa_RE;
    f_static_array__array_RE_group_mNa_RE.open("static_arrays/_static_array__array_RE_group_mNa_RE", ios::in | ios::binary);
    if(f_static_array__array_RE_group_mNa_RE.is_open())
    {
        f_static_array__array_RE_group_mNa_RE.read(reinterpret_cast<char*>(_static_array__array_RE_group_mNa_RE), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_RE_group_mNa_RE." << endl;
    }
    ifstream f_static_array__array_RE_group_mT_RE;
    f_static_array__array_RE_group_mT_RE.open("static_arrays/_static_array__array_RE_group_mT_RE", ios::in | ios::binary);
    if(f_static_array__array_RE_group_mT_RE.is_open())
    {
        f_static_array__array_RE_group_mT_RE.read(reinterpret_cast<char*>(_static_array__array_RE_group_mT_RE), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_RE_group_mT_RE." << endl;
    }
    ifstream f_static_array__array_RE_group_nK_RE;
    f_static_array__array_RE_group_nK_RE.open("static_arrays/_static_array__array_RE_group_nK_RE", ios::in | ios::binary);
    if(f_static_array__array_RE_group_nK_RE.is_open())
    {
        f_static_array__array_RE_group_nK_RE.read(reinterpret_cast<char*>(_static_array__array_RE_group_nK_RE), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_RE_group_nK_RE." << endl;
    }
    ifstream f_static_array__array_RE_group_v;
    f_static_array__array_RE_group_v.open("static_arrays/_static_array__array_RE_group_v", ios::in | ios::binary);
    if(f_static_array__array_RE_group_v.is_open())
    {
        f_static_array__array_RE_group_v.read(reinterpret_cast<char*>(_static_array__array_RE_group_v), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_RE_group_v." << endl;
    }
    ifstream f_static_array__array_TC_group_Ca_TC;
    f_static_array__array_TC_group_Ca_TC.open("static_arrays/_static_array__array_TC_group_Ca_TC", ios::in | ios::binary);
    if(f_static_array__array_TC_group_Ca_TC.is_open())
    {
        f_static_array__array_TC_group_Ca_TC.read(reinterpret_cast<char*>(_static_array__array_TC_group_Ca_TC), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_Ca_TC." << endl;
    }
    ifstream f_static_array__array_TC_group_Open;
    f_static_array__array_TC_group_Open.open("static_arrays/_static_array__array_TC_group_Open", ios::in | ios::binary);
    if(f_static_array__array_TC_group_Open.is_open())
    {
        f_static_array__array_TC_group_Open.read(reinterpret_cast<char*>(_static_array__array_TC_group_Open), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_Open." << endl;
    }
    ifstream f_static_array__array_TC_group_OpenLocked;
    f_static_array__array_TC_group_OpenLocked.open("static_arrays/_static_array__array_TC_group_OpenLocked", ios::in | ios::binary);
    if(f_static_array__array_TC_group_OpenLocked.is_open())
    {
        f_static_array__array_TC_group_OpenLocked.read(reinterpret_cast<char*>(_static_array__array_TC_group_OpenLocked), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_OpenLocked." << endl;
    }
    ifstream f_static_array__array_TC_group_Pone;
    f_static_array__array_TC_group_Pone.open("static_arrays/_static_array__array_TC_group_Pone", ios::in | ios::binary);
    if(f_static_array__array_TC_group_Pone.is_open())
    {
        f_static_array__array_TC_group_Pone.read(reinterpret_cast<char*>(_static_array__array_TC_group_Pone), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_Pone." << endl;
    }
    ifstream f_static_array__array_TC_group_hNa_TC;
    f_static_array__array_TC_group_hNa_TC.open("static_arrays/_static_array__array_TC_group_hNa_TC", ios::in | ios::binary);
    if(f_static_array__array_TC_group_hNa_TC.is_open())
    {
        f_static_array__array_TC_group_hNa_TC.read(reinterpret_cast<char*>(_static_array__array_TC_group_hNa_TC), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_hNa_TC." << endl;
    }
    ifstream f_static_array__array_TC_group_hT_TC;
    f_static_array__array_TC_group_hT_TC.open("static_arrays/_static_array__array_TC_group_hT_TC", ios::in | ios::binary);
    if(f_static_array__array_TC_group_hT_TC.is_open())
    {
        f_static_array__array_TC_group_hT_TC.read(reinterpret_cast<char*>(_static_array__array_TC_group_hT_TC), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_hT_TC." << endl;
    }
    ifstream f_static_array__array_TC_group_mNa_TC;
    f_static_array__array_TC_group_mNa_TC.open("static_arrays/_static_array__array_TC_group_mNa_TC", ios::in | ios::binary);
    if(f_static_array__array_TC_group_mNa_TC.is_open())
    {
        f_static_array__array_TC_group_mNa_TC.read(reinterpret_cast<char*>(_static_array__array_TC_group_mNa_TC), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_mNa_TC." << endl;
    }
    ifstream f_static_array__array_TC_group_nK_TC;
    f_static_array__array_TC_group_nK_TC.open("static_arrays/_static_array__array_TC_group_nK_TC", ios::in | ios::binary);
    if(f_static_array__array_TC_group_nK_TC.is_open())
    {
        f_static_array__array_TC_group_nK_TC.read(reinterpret_cast<char*>(_static_array__array_TC_group_nK_TC), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_nK_TC." << endl;
    }
    ifstream f_static_array__array_TC_group_v;
    f_static_array__array_TC_group_v.open("static_arrays/_static_array__array_TC_group_v", ios::in | ios::binary);
    if(f_static_array__array_TC_group_v.is_open())
    {
        f_static_array__array_TC_group_v.read(reinterpret_cast<char*>(_static_array__array_TC_group_v), 20*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__array_TC_group_v." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA;
    f_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA;
    f_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_PYso_IN_sAMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
    f_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA;
    f_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
    f_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.open("static_arrays/_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.is_open())
    {
        f_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA;
    f_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_PYso_TC_sAMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_TC_IN_sAMPA;
    f_static_array__dynamic_array_IAMPA_TC_IN_sAMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_TC_IN_sAMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_TC_IN_sAMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_TC_IN_sAMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_TC_IN_sAMPA), 400*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_TC_IN_sAMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA;
    f_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA), 420*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_IAMPA_TC_RE_sAMPA;
    f_static_array__dynamic_array_IAMPA_TC_RE_sAMPA.open("static_arrays/_static_array__dynamic_array_IAMPA_TC_RE_sAMPA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IAMPA_TC_RE_sAMPA.is_open())
    {
        f_static_array__dynamic_array_IAMPA_TC_RE_sAMPA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IAMPA_TC_RE_sAMPA), 400*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IAMPA_TC_RE_sAMPA." << endl;
    }
    ifstream f_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM;
    f_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM.open("static_arrays/_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM", ios::in | ios::binary);
    if(f_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM.is_open())
    {
        f_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM.read(reinterpret_cast<char*>(_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_ICOM_PYdr_PYso_gCOM." << endl;
    }
    ifstream f_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM;
    f_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM.open("static_arrays/_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM", ios::in | ios::binary);
    if(f_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM.is_open())
    {
        f_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM.read(reinterpret_cast<char*>(_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_ICOM_PYso_PYdr_gCOM." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA;
    f_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA.open("static_arrays/_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA.is_open())
    {
        f_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA), 380*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAA_IN_IN_res_GABAA." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAA_IN_IN_sGABAA;
    f_static_array__dynamic_array_IGABAA_IN_IN_sGABAA.open("static_arrays/_static_array__dynamic_array_IGABAA_IN_IN_sGABAA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAA_IN_IN_sGABAA.is_open())
    {
        f_static_array__dynamic_array_IGABAA_IN_IN_sGABAA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAA_IN_IN_sGABAA), 380*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAA_IN_IN_sGABAA." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA;
    f_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA.open("static_arrays/_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA.is_open())
    {
        f_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA), 420*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA;
    f_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA.open("static_arrays/_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA.is_open())
    {
        f_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA), 420*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAA_IN_PYso_sGABAA." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAA_RE_RE_sGABAA;
    f_static_array__dynamic_array_IGABAA_RE_RE_sGABAA.open("static_arrays/_static_array__dynamic_array_IGABAA_RE_RE_sGABAA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAA_RE_RE_sGABAA.is_open())
    {
        f_static_array__dynamic_array_IGABAA_RE_RE_sGABAA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAA_RE_RE_sGABAA), 400*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAA_RE_RE_sGABAA." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAA_RE_TC_sGABAA;
    f_static_array__dynamic_array_IGABAA_RE_TC_sGABAA.open("static_arrays/_static_array__dynamic_array_IGABAA_RE_TC_sGABAA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAA_RE_TC_sGABAA.is_open())
    {
        f_static_array__dynamic_array_IGABAA_RE_TC_sGABAA.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAA_RE_TC_sGABAA), 400*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAA_RE_TC_sGABAA." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAB_RE_TC_rGABAB;
    f_static_array__dynamic_array_IGABAB_RE_TC_rGABAB.open("static_arrays/_static_array__dynamic_array_IGABAB_RE_TC_rGABAB", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAB_RE_TC_rGABAB.is_open())
    {
        f_static_array__dynamic_array_IGABAB_RE_TC_rGABAB.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAB_RE_TC_rGABAB), 400*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAB_RE_TC_rGABAB." << endl;
    }
    ifstream f_static_array__dynamic_array_IGABAB_RE_TC_sGABAB;
    f_static_array__dynamic_array_IGABAB_RE_TC_sGABAB.open("static_arrays/_static_array__dynamic_array_IGABAB_RE_TC_sGABAB", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IGABAB_RE_TC_sGABAB.is_open())
    {
        f_static_array__dynamic_array_IGABAB_RE_TC_sGABAB.read(reinterpret_cast<char*>(_static_array__dynamic_array_IGABAB_RE_TC_sGABAB), 400*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IGABAB_RE_TC_sGABAB." << endl;
    }
    ifstream f_static_array__dynamic_array_IKNa_PYso_PYdr_concNa;
    f_static_array__dynamic_array_IKNa_PYso_PYdr_concNa.open("static_arrays/_static_array__dynamic_array_IKNa_PYso_PYdr_concNa", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IKNa_PYso_PYdr_concNa.is_open())
    {
        f_static_array__dynamic_array_IKNa_PYso_PYdr_concNa.read(reinterpret_cast<char*>(_static_array__dynamic_array_IKNa_PYso_PYdr_concNa), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IKNa_PYso_PYdr_concNa." << endl;
    }
    ifstream f_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local;
    f_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local.open("static_arrays/_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local", ios::in | ios::binary);
    if(f_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local.is_open())
    {
        f_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local.read(reinterpret_cast<char*>(_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local), 100*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local." << endl;
    }
    ifstream f_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA;
    f_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA.open("static_arrays/_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA.is_open())
    {
        f_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA.read(reinterpret_cast<char*>(_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_INMDA_PYso_IN_res_NMDA." << endl;
    }
    ifstream f_static_array__dynamic_array_INMDA_PYso_IN_sNMDA;
    f_static_array__dynamic_array_INMDA_PYso_IN_sNMDA.open("static_arrays/_static_array__dynamic_array_INMDA_PYso_IN_sNMDA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_INMDA_PYso_IN_sNMDA.is_open())
    {
        f_static_array__dynamic_array_INMDA_PYso_IN_sNMDA.read(reinterpret_cast<char*>(_static_array__dynamic_array_INMDA_PYso_IN_sNMDA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_INMDA_PYso_IN_sNMDA." << endl;
    }
    ifstream f_static_array__dynamic_array_INMDA_PYso_IN_xNMDA;
    f_static_array__dynamic_array_INMDA_PYso_IN_xNMDA.open("static_arrays/_static_array__dynamic_array_INMDA_PYso_IN_xNMDA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_INMDA_PYso_IN_xNMDA.is_open())
    {
        f_static_array__dynamic_array_INMDA_PYso_IN_xNMDA.read(reinterpret_cast<char*>(_static_array__dynamic_array_INMDA_PYso_IN_xNMDA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_INMDA_PYso_IN_xNMDA." << endl;
    }
    ifstream f_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA;
    f_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA.open("static_arrays/_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA.is_open())
    {
        f_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA.read(reinterpret_cast<char*>(_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA." << endl;
    }
    ifstream f_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA;
    f_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA.open("static_arrays/_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA.is_open())
    {
        f_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA.read(reinterpret_cast<char*>(_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA." << endl;
    }
    ifstream f_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA;
    f_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA.open("static_arrays/_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA", ios::in | ios::binary);
    if(f_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA.is_open())
    {
        f_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA.read(reinterpret_cast<char*>(_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA), 2000*sizeof(double));
    } else
    {
        std::cout << "Error opening static array _static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA." << endl;
    }
}

void _write_arrays()
{
    using namespace brian;

    ofstream outfile__array_clock_dt;
    outfile__array_clock_dt.open(results_dir + "_array_clock_dt_3495478892", ios::binary | ios::out);
    if(outfile__array_clock_dt.is_open())
    {
        outfile__array_clock_dt.write(reinterpret_cast<char*>(_array_clock_dt), 1*sizeof(_array_clock_dt[0]));
        outfile__array_clock_dt.close();
    } else
    {
        std::cout << "Error writing output file for _array_clock_dt." << endl;
    }
    ofstream outfile__array_clock_t;
    outfile__array_clock_t.open(results_dir + "_array_clock_t_705559507", ios::binary | ios::out);
    if(outfile__array_clock_t.is_open())
    {
        outfile__array_clock_t.write(reinterpret_cast<char*>(_array_clock_t), 1*sizeof(_array_clock_t[0]));
        outfile__array_clock_t.close();
    } else
    {
        std::cout << "Error writing output file for _array_clock_t." << endl;
    }
    ofstream outfile__array_clock_timestep;
    outfile__array_clock_timestep.open(results_dir + "_array_clock_timestep_3597035576", ios::binary | ios::out);
    if(outfile__array_clock_timestep.is_open())
    {
        outfile__array_clock_timestep.write(reinterpret_cast<char*>(_array_clock_timestep), 1*sizeof(_array_clock_timestep[0]));
        outfile__array_clock_timestep.close();
    } else
    {
        std::cout << "Error writing output file for _array_clock_timestep." << endl;
    }
    ofstream outfile__array_defaultclock_dt;
    outfile__array_defaultclock_dt.open(results_dir + "_array_defaultclock_dt_1978099143", ios::binary | ios::out);
    if(outfile__array_defaultclock_dt.is_open())
    {
        outfile__array_defaultclock_dt.write(reinterpret_cast<char*>(_array_defaultclock_dt), 1*sizeof(_array_defaultclock_dt[0]));
        outfile__array_defaultclock_dt.close();
    } else
    {
        std::cout << "Error writing output file for _array_defaultclock_dt." << endl;
    }
    ofstream outfile__array_defaultclock_t;
    outfile__array_defaultclock_t.open(results_dir + "_array_defaultclock_t_2669362164", ios::binary | ios::out);
    if(outfile__array_defaultclock_t.is_open())
    {
        outfile__array_defaultclock_t.write(reinterpret_cast<char*>(_array_defaultclock_t), 1*sizeof(_array_defaultclock_t[0]));
        outfile__array_defaultclock_t.close();
    } else
    {
        std::cout << "Error writing output file for _array_defaultclock_t." << endl;
    }
    ofstream outfile__array_defaultclock_timestep;
    outfile__array_defaultclock_timestep.open(results_dir + "_array_defaultclock_timestep_144223508", ios::binary | ios::out);
    if(outfile__array_defaultclock_timestep.is_open())
    {
        outfile__array_defaultclock_timestep.write(reinterpret_cast<char*>(_array_defaultclock_timestep), 1*sizeof(_array_defaultclock_timestep[0]));
        outfile__array_defaultclock_timestep.close();
    } else
    {
        std::cout << "Error writing output file for _array_defaultclock_timestep." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_IN_N;
    outfile__array_IAMPA_PYso_IN_N.open(results_dir + "_array_IAMPA_PYso_IN_N_3332304076", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_IN_N.is_open())
    {
        outfile__array_IAMPA_PYso_IN_N.write(reinterpret_cast<char*>(_array_IAMPA_PYso_IN_N), 1*sizeof(_array_IAMPA_PYso_IN_N[0]));
        outfile__array_IAMPA_PYso_IN_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_IN_N." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_IN_sources;
    outfile__array_IAMPA_PYso_IN_sources.open(results_dir + "_array_IAMPA_PYso_IN_sources_2842991999", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_IN_sources.is_open())
    {
        outfile__array_IAMPA_PYso_IN_sources.write(reinterpret_cast<char*>(_array_IAMPA_PYso_IN_sources), 2000*sizeof(_array_IAMPA_PYso_IN_sources[0]));
        outfile__array_IAMPA_PYso_IN_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_IN_sources." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_IN_targets;
    outfile__array_IAMPA_PYso_IN_targets.open(results_dir + "_array_IAMPA_PYso_IN_targets_3563775646", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_IN_targets.is_open())
    {
        outfile__array_IAMPA_PYso_IN_targets.write(reinterpret_cast<char*>(_array_IAMPA_PYso_IN_targets), 2000*sizeof(_array_IAMPA_PYso_IN_targets[0]));
        outfile__array_IAMPA_PYso_IN_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_IN_targets." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_PYdr_N;
    outfile__array_IAMPA_PYso_PYdr_N.open(results_dir + "_array_IAMPA_PYso_PYdr_N_3744445673", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_PYdr_N.is_open())
    {
        outfile__array_IAMPA_PYso_PYdr_N.write(reinterpret_cast<char*>(_array_IAMPA_PYso_PYdr_N), 1*sizeof(_array_IAMPA_PYso_PYdr_N[0]));
        outfile__array_IAMPA_PYso_PYdr_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_PYdr_N." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_PYdr_sources;
    outfile__array_IAMPA_PYso_PYdr_sources.open(results_dir + "_array_IAMPA_PYso_PYdr_sources_1078934768", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_PYdr_sources.is_open())
    {
        outfile__array_IAMPA_PYso_PYdr_sources.write(reinterpret_cast<char*>(_array_IAMPA_PYso_PYdr_sources), 2000*sizeof(_array_IAMPA_PYso_PYdr_sources[0]));
        outfile__array_IAMPA_PYso_PYdr_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_PYdr_sources." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_PYdr_targets;
    outfile__array_IAMPA_PYso_PYdr_targets.open(results_dir + "_array_IAMPA_PYso_PYdr_targets_1028736785", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_PYdr_targets.is_open())
    {
        outfile__array_IAMPA_PYso_PYdr_targets.write(reinterpret_cast<char*>(_array_IAMPA_PYso_PYdr_targets), 2000*sizeof(_array_IAMPA_PYso_PYdr_targets[0]));
        outfile__array_IAMPA_PYso_PYdr_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_PYdr_targets." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_RE_N;
    outfile__array_IAMPA_PYso_RE_N.open(results_dir + "_array_IAMPA_PYso_RE_N_1305727923", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_RE_N.is_open())
    {
        outfile__array_IAMPA_PYso_RE_N.write(reinterpret_cast<char*>(_array_IAMPA_PYso_RE_N), 1*sizeof(_array_IAMPA_PYso_RE_N[0]));
        outfile__array_IAMPA_PYso_RE_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_RE_N." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_RE_sources;
    outfile__array_IAMPA_PYso_RE_sources.open(results_dir + "_array_IAMPA_PYso_RE_sources_4100946365", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_RE_sources.is_open())
    {
        outfile__array_IAMPA_PYso_RE_sources.write(reinterpret_cast<char*>(_array_IAMPA_PYso_RE_sources), 2000*sizeof(_array_IAMPA_PYso_RE_sources[0]));
        outfile__array_IAMPA_PYso_RE_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_RE_sources." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_RE_targets;
    outfile__array_IAMPA_PYso_RE_targets.open(results_dir + "_array_IAMPA_PYso_RE_targets_2305884764", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_RE_targets.is_open())
    {
        outfile__array_IAMPA_PYso_RE_targets.write(reinterpret_cast<char*>(_array_IAMPA_PYso_RE_targets), 2000*sizeof(_array_IAMPA_PYso_RE_targets[0]));
        outfile__array_IAMPA_PYso_RE_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_RE_targets." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_TC_N;
    outfile__array_IAMPA_PYso_TC_N.open(results_dir + "_array_IAMPA_PYso_TC_N_1815474397", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_TC_N.is_open())
    {
        outfile__array_IAMPA_PYso_TC_N.write(reinterpret_cast<char*>(_array_IAMPA_PYso_TC_N), 1*sizeof(_array_IAMPA_PYso_TC_N[0]));
        outfile__array_IAMPA_PYso_TC_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_TC_N." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_TC_sources;
    outfile__array_IAMPA_PYso_TC_sources.open(results_dir + "_array_IAMPA_PYso_TC_sources_2322296944", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_TC_sources.is_open())
    {
        outfile__array_IAMPA_PYso_TC_sources.write(reinterpret_cast<char*>(_array_IAMPA_PYso_TC_sources), 2000*sizeof(_array_IAMPA_PYso_TC_sources[0]));
        outfile__array_IAMPA_PYso_TC_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_TC_sources." << endl;
    }
    ofstream outfile__array_IAMPA_PYso_TC_targets;
    outfile__array_IAMPA_PYso_TC_targets.open(results_dir + "_array_IAMPA_PYso_TC_targets_4151643025", ios::binary | ios::out);
    if(outfile__array_IAMPA_PYso_TC_targets.is_open())
    {
        outfile__array_IAMPA_PYso_TC_targets.write(reinterpret_cast<char*>(_array_IAMPA_PYso_TC_targets), 2000*sizeof(_array_IAMPA_PYso_TC_targets[0]));
        outfile__array_IAMPA_PYso_TC_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_PYso_TC_targets." << endl;
    }
    ofstream outfile__array_IAMPA_TC_IN_N;
    outfile__array_IAMPA_TC_IN_N.open(results_dir + "_array_IAMPA_TC_IN_N_1409664613", ios::binary | ios::out);
    if(outfile__array_IAMPA_TC_IN_N.is_open())
    {
        outfile__array_IAMPA_TC_IN_N.write(reinterpret_cast<char*>(_array_IAMPA_TC_IN_N), 1*sizeof(_array_IAMPA_TC_IN_N[0]));
        outfile__array_IAMPA_TC_IN_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_TC_IN_N." << endl;
    }
    ofstream outfile__array_IAMPA_TC_IN_sources;
    outfile__array_IAMPA_TC_IN_sources.open(results_dir + "_array_IAMPA_TC_IN_sources_2294341831", ios::binary | ios::out);
    if(outfile__array_IAMPA_TC_IN_sources.is_open())
    {
        outfile__array_IAMPA_TC_IN_sources.write(reinterpret_cast<char*>(_array_IAMPA_TC_IN_sources), 400*sizeof(_array_IAMPA_TC_IN_sources[0]));
        outfile__array_IAMPA_TC_IN_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_TC_IN_sources." << endl;
    }
    ofstream outfile__array_IAMPA_TC_IN_targets;
    outfile__array_IAMPA_TC_IN_targets.open(results_dir + "_array_IAMPA_TC_IN_targets_4125006630", ios::binary | ios::out);
    if(outfile__array_IAMPA_TC_IN_targets.is_open())
    {
        outfile__array_IAMPA_TC_IN_targets.write(reinterpret_cast<char*>(_array_IAMPA_TC_IN_targets), 400*sizeof(_array_IAMPA_TC_IN_targets[0]));
        outfile__array_IAMPA_TC_IN_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_TC_IN_targets." << endl;
    }
    ofstream outfile__array_IAMPA_TC_PYdr_N;
    outfile__array_IAMPA_TC_PYdr_N.open(results_dir + "_array_IAMPA_TC_PYdr_N_2068965380", ios::binary | ios::out);
    if(outfile__array_IAMPA_TC_PYdr_N.is_open())
    {
        outfile__array_IAMPA_TC_PYdr_N.write(reinterpret_cast<char*>(_array_IAMPA_TC_PYdr_N), 1*sizeof(_array_IAMPA_TC_PYdr_N[0]));
        outfile__array_IAMPA_TC_PYdr_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_TC_PYdr_N." << endl;
    }
    ofstream outfile__array_IAMPA_TC_PYdr_sources;
    outfile__array_IAMPA_TC_PYdr_sources.open(results_dir + "_array_IAMPA_TC_PYdr_sources_1122448501", ios::binary | ios::out);
    if(outfile__array_IAMPA_TC_PYdr_sources.is_open())
    {
        outfile__array_IAMPA_TC_PYdr_sources.write(reinterpret_cast<char*>(_array_IAMPA_TC_PYdr_sources), 420*sizeof(_array_IAMPA_TC_PYdr_sources[0]));
        outfile__array_IAMPA_TC_PYdr_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_TC_PYdr_sources." << endl;
    }
    ofstream outfile__array_IAMPA_TC_PYdr_targets;
    outfile__array_IAMPA_TC_PYdr_targets.open(results_dir + "_array_IAMPA_TC_PYdr_targets_1073303444", ios::binary | ios::out);
    if(outfile__array_IAMPA_TC_PYdr_targets.is_open())
    {
        outfile__array_IAMPA_TC_PYdr_targets.write(reinterpret_cast<char*>(_array_IAMPA_TC_PYdr_targets), 420*sizeof(_array_IAMPA_TC_PYdr_targets[0]));
        outfile__array_IAMPA_TC_PYdr_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_TC_PYdr_targets." << endl;
    }
    ofstream outfile__array_IAMPA_TC_RE_N;
    outfile__array_IAMPA_TC_RE_N.open(results_dir + "_array_IAMPA_TC_RE_N_3746101530", ios::binary | ios::out);
    if(outfile__array_IAMPA_TC_RE_N.is_open())
    {
        outfile__array_IAMPA_TC_RE_N.write(reinterpret_cast<char*>(_array_IAMPA_TC_RE_N), 1*sizeof(_array_IAMPA_TC_RE_N[0]));
        outfile__array_IAMPA_TC_RE_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IAMPA_TC_RE_N." << endl;
    }
    ofstream outfile__array_ICOM_PYdr_PYso_N;
    outfile__array_ICOM_PYdr_PYso_N.open(results_dir + "_array_ICOM_PYdr_PYso_N_3196551417", ios::binary | ios::out);
    if(outfile__array_ICOM_PYdr_PYso_N.is_open())
    {
        outfile__array_ICOM_PYdr_PYso_N.write(reinterpret_cast<char*>(_array_ICOM_PYdr_PYso_N), 1*sizeof(_array_ICOM_PYdr_PYso_N[0]));
        outfile__array_ICOM_PYdr_PYso_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_ICOM_PYdr_PYso_N." << endl;
    }
    ofstream outfile__array_ICOM_PYso_PYdr_N;
    outfile__array_ICOM_PYso_PYdr_N.open(results_dir + "_array_ICOM_PYso_PYdr_N_1219024062", ios::binary | ios::out);
    if(outfile__array_ICOM_PYso_PYdr_N.is_open())
    {
        outfile__array_ICOM_PYso_PYdr_N.write(reinterpret_cast<char*>(_array_ICOM_PYso_PYdr_N), 1*sizeof(_array_ICOM_PYso_PYdr_N[0]));
        outfile__array_ICOM_PYso_PYdr_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_ICOM_PYso_PYdr_N." << endl;
    }
    ofstream outfile__array_IGABAA_IN_IN_N;
    outfile__array_IGABAA_IN_IN_N.open(results_dir + "_array_IGABAA_IN_IN_N_1145522004", ios::binary | ios::out);
    if(outfile__array_IGABAA_IN_IN_N.is_open())
    {
        outfile__array_IGABAA_IN_IN_N.write(reinterpret_cast<char*>(_array_IGABAA_IN_IN_N), 1*sizeof(_array_IGABAA_IN_IN_N[0]));
        outfile__array_IGABAA_IN_IN_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_IN_IN_N." << endl;
    }
    ofstream outfile__array_IGABAA_IN_IN_sources;
    outfile__array_IGABAA_IN_IN_sources.open(results_dir + "_array_IGABAA_IN_IN_sources_3864763947", ios::binary | ios::out);
    if(outfile__array_IGABAA_IN_IN_sources.is_open())
    {
        outfile__array_IGABAA_IN_IN_sources.write(reinterpret_cast<char*>(_array_IGABAA_IN_IN_sources), 380*sizeof(_array_IGABAA_IN_IN_sources[0]));
        outfile__array_IGABAA_IN_IN_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_IN_IN_sources." << endl;
    }
    ofstream outfile__array_IGABAA_IN_IN_targets;
    outfile__array_IGABAA_IN_IN_targets.open(results_dir + "_array_IGABAA_IN_IN_targets_2605049290", ios::binary | ios::out);
    if(outfile__array_IGABAA_IN_IN_targets.is_open())
    {
        outfile__array_IGABAA_IN_IN_targets.write(reinterpret_cast<char*>(_array_IGABAA_IN_IN_targets), 380*sizeof(_array_IGABAA_IN_IN_targets[0]));
        outfile__array_IGABAA_IN_IN_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_IN_IN_targets." << endl;
    }
    ofstream outfile__array_IGABAA_IN_PYso_N;
    outfile__array_IGABAA_IN_PYso_N.open(results_dir + "_array_IGABAA_IN_PYso_N_4270003015", ios::binary | ios::out);
    if(outfile__array_IGABAA_IN_PYso_N.is_open())
    {
        outfile__array_IGABAA_IN_PYso_N.write(reinterpret_cast<char*>(_array_IGABAA_IN_PYso_N), 1*sizeof(_array_IGABAA_IN_PYso_N[0]));
        outfile__array_IGABAA_IN_PYso_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_IN_PYso_N." << endl;
    }
    ofstream outfile__array_IGABAA_IN_PYso_sources;
    outfile__array_IGABAA_IN_PYso_sources.open(results_dir + "_array_IGABAA_IN_PYso_sources_2674689410", ios::binary | ios::out);
    if(outfile__array_IGABAA_IN_PYso_sources.is_open())
    {
        outfile__array_IGABAA_IN_PYso_sources.write(reinterpret_cast<char*>(_array_IGABAA_IN_PYso_sources), 420*sizeof(_array_IGABAA_IN_PYso_sources[0]));
        outfile__array_IGABAA_IN_PYso_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_IN_PYso_sources." << endl;
    }
    ofstream outfile__array_IGABAA_IN_PYso_targets;
    outfile__array_IGABAA_IN_PYso_targets.open(results_dir + "_array_IGABAA_IN_PYso_targets_3799187043", ios::binary | ios::out);
    if(outfile__array_IGABAA_IN_PYso_targets.is_open())
    {
        outfile__array_IGABAA_IN_PYso_targets.write(reinterpret_cast<char*>(_array_IGABAA_IN_PYso_targets), 420*sizeof(_array_IGABAA_IN_PYso_targets[0]));
        outfile__array_IGABAA_IN_PYso_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_IN_PYso_targets." << endl;
    }
    ofstream outfile__array_IGABAA_RE_RE_N;
    outfile__array_IGABAA_RE_RE_N.open(results_dir + "_array_IGABAA_RE_RE_N_3222552796", ios::binary | ios::out);
    if(outfile__array_IGABAA_RE_RE_N.is_open())
    {
        outfile__array_IGABAA_RE_RE_N.write(reinterpret_cast<char*>(_array_IGABAA_RE_RE_N), 1*sizeof(_array_IGABAA_RE_RE_N[0]));
        outfile__array_IGABAA_RE_RE_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_RE_RE_N." << endl;
    }
    ofstream outfile__array_IGABAA_RE_TC_N;
    outfile__array_IGABAA_RE_TC_N.open(results_dir + "_array_IGABAA_RE_TC_N_3790758834", ios::binary | ios::out);
    if(outfile__array_IGABAA_RE_TC_N.is_open())
    {
        outfile__array_IGABAA_RE_TC_N.write(reinterpret_cast<char*>(_array_IGABAA_RE_TC_N), 1*sizeof(_array_IGABAA_RE_TC_N[0]));
        outfile__array_IGABAA_RE_TC_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAA_RE_TC_N." << endl;
    }
    ofstream outfile__array_IGABAB_RE_TC_N;
    outfile__array_IGABAB_RE_TC_N.open(results_dir + "_array_IGABAB_RE_TC_N_3632226167", ios::binary | ios::out);
    if(outfile__array_IGABAB_RE_TC_N.is_open())
    {
        outfile__array_IGABAB_RE_TC_N.write(reinterpret_cast<char*>(_array_IGABAB_RE_TC_N), 1*sizeof(_array_IGABAB_RE_TC_N[0]));
        outfile__array_IGABAB_RE_TC_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IGABAB_RE_TC_N." << endl;
    }
    ofstream outfile__array_IKNa_PYso_PYdr_N;
    outfile__array_IKNa_PYso_PYdr_N.open(results_dir + "_array_IKNa_PYso_PYdr_N_870976321", ios::binary | ios::out);
    if(outfile__array_IKNa_PYso_PYdr_N.is_open())
    {
        outfile__array_IKNa_PYso_PYdr_N.write(reinterpret_cast<char*>(_array_IKNa_PYso_PYdr_N), 1*sizeof(_array_IKNa_PYso_PYdr_N[0]));
        outfile__array_IKNa_PYso_PYdr_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IKNa_PYso_PYdr_N." << endl;
    }
    ofstream outfile__array_IN_group__spikespace;
    outfile__array_IN_group__spikespace.open(results_dir + "_array_IN_group__spikespace_3415669226", ios::binary | ios::out);
    if(outfile__array_IN_group__spikespace.is_open())
    {
        outfile__array_IN_group__spikespace.write(reinterpret_cast<char*>(_array_IN_group__spikespace), 21*sizeof(_array_IN_group__spikespace[0]));
        outfile__array_IN_group__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group__spikespace." << endl;
    }
    ofstream outfile__array_IN_group_hNa_IN;
    outfile__array_IN_group_hNa_IN.open(results_dir + "_array_IN_group_hNa_IN_4289110967", ios::binary | ios::out);
    if(outfile__array_IN_group_hNa_IN.is_open())
    {
        outfile__array_IN_group_hNa_IN.write(reinterpret_cast<char*>(_array_IN_group_hNa_IN), 20*sizeof(_array_IN_group_hNa_IN[0]));
        outfile__array_IN_group_hNa_IN.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_hNa_IN." << endl;
    }
    ofstream outfile__array_IN_group_i;
    outfile__array_IN_group_i.open(results_dir + "_array_IN_group_i_477439055", ios::binary | ios::out);
    if(outfile__array_IN_group_i.is_open())
    {
        outfile__array_IN_group_i.write(reinterpret_cast<char*>(_array_IN_group_i), 20*sizeof(_array_IN_group_i[0]));
        outfile__array_IN_group_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_i." << endl;
    }
    ofstream outfile__array_IN_group_iAMPA_PYso_IN;
    outfile__array_IN_group_iAMPA_PYso_IN.open(results_dir + "_array_IN_group_iAMPA_PYso_IN_491654791", ios::binary | ios::out);
    if(outfile__array_IN_group_iAMPA_PYso_IN.is_open())
    {
        outfile__array_IN_group_iAMPA_PYso_IN.write(reinterpret_cast<char*>(_array_IN_group_iAMPA_PYso_IN), 20*sizeof(_array_IN_group_iAMPA_PYso_IN[0]));
        outfile__array_IN_group_iAMPA_PYso_IN.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_iAMPA_PYso_IN." << endl;
    }
    ofstream outfile__array_IN_group_iAMPA_TC_IN;
    outfile__array_IN_group_iAMPA_TC_IN.open(results_dir + "_array_IN_group_iAMPA_TC_IN_2428657363", ios::binary | ios::out);
    if(outfile__array_IN_group_iAMPA_TC_IN.is_open())
    {
        outfile__array_IN_group_iAMPA_TC_IN.write(reinterpret_cast<char*>(_array_IN_group_iAMPA_TC_IN), 20*sizeof(_array_IN_group_iAMPA_TC_IN[0]));
        outfile__array_IN_group_iAMPA_TC_IN.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_iAMPA_TC_IN." << endl;
    }
    ofstream outfile__array_IN_group_iGABAA_IN_IN;
    outfile__array_IN_group_iGABAA_IN_IN.open(results_dir + "_array_IN_group_iGABAA_IN_IN_2219504670", ios::binary | ios::out);
    if(outfile__array_IN_group_iGABAA_IN_IN.is_open())
    {
        outfile__array_IN_group_iGABAA_IN_IN.write(reinterpret_cast<char*>(_array_IN_group_iGABAA_IN_IN), 20*sizeof(_array_IN_group_iGABAA_IN_IN[0]));
        outfile__array_IN_group_iGABAA_IN_IN.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_iGABAA_IN_IN." << endl;
    }
    ofstream outfile__array_IN_group_iNMDA_PYso_IN;
    outfile__array_IN_group_iNMDA_PYso_IN.open(results_dir + "_array_IN_group_iNMDA_PYso_IN_1100984420", ios::binary | ios::out);
    if(outfile__array_IN_group_iNMDA_PYso_IN.is_open())
    {
        outfile__array_IN_group_iNMDA_PYso_IN.write(reinterpret_cast<char*>(_array_IN_group_iNMDA_PYso_IN), 20*sizeof(_array_IN_group_iNMDA_PYso_IN[0]));
        outfile__array_IN_group_iNMDA_PYso_IN.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_iNMDA_PYso_IN." << endl;
    }
    ofstream outfile__array_IN_group_lastspike;
    outfile__array_IN_group_lastspike.open(results_dir + "_array_IN_group_lastspike_1659730378", ios::binary | ios::out);
    if(outfile__array_IN_group_lastspike.is_open())
    {
        outfile__array_IN_group_lastspike.write(reinterpret_cast<char*>(_array_IN_group_lastspike), 20*sizeof(_array_IN_group_lastspike[0]));
        outfile__array_IN_group_lastspike.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_lastspike." << endl;
    }
    ofstream outfile__array_IN_group_nK_IN;
    outfile__array_IN_group_nK_IN.open(results_dir + "_array_IN_group_nK_IN_2989775565", ios::binary | ios::out);
    if(outfile__array_IN_group_nK_IN.is_open())
    {
        outfile__array_IN_group_nK_IN.write(reinterpret_cast<char*>(_array_IN_group_nK_IN), 20*sizeof(_array_IN_group_nK_IN[0]));
        outfile__array_IN_group_nK_IN.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_nK_IN." << endl;
    }
    ofstream outfile__array_IN_group_not_refractory;
    outfile__array_IN_group_not_refractory.open(results_dir + "_array_IN_group_not_refractory_3139046054", ios::binary | ios::out);
    if(outfile__array_IN_group_not_refractory.is_open())
    {
        outfile__array_IN_group_not_refractory.write(reinterpret_cast<char*>(_array_IN_group_not_refractory), 20*sizeof(_array_IN_group_not_refractory[0]));
        outfile__array_IN_group_not_refractory.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_not_refractory." << endl;
    }
    ofstream outfile__array_IN_group_v;
    outfile__array_IN_group_v.open(results_dir + "_array_IN_group_v_2440899002", ios::binary | ios::out);
    if(outfile__array_IN_group_v.is_open())
    {
        outfile__array_IN_group_v.write(reinterpret_cast<char*>(_array_IN_group_v), 20*sizeof(_array_IN_group_v[0]));
        outfile__array_IN_group_v.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_group_v." << endl;
    }
    ofstream outfile__array_IN_spikemon__source_idx;
    outfile__array_IN_spikemon__source_idx.open(results_dir + "_array_IN_spikemon__source_idx_3130554247", ios::binary | ios::out);
    if(outfile__array_IN_spikemon__source_idx.is_open())
    {
        outfile__array_IN_spikemon__source_idx.write(reinterpret_cast<char*>(_array_IN_spikemon__source_idx), 20*sizeof(_array_IN_spikemon__source_idx[0]));
        outfile__array_IN_spikemon__source_idx.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_spikemon__source_idx." << endl;
    }
    ofstream outfile__array_IN_spikemon_count;
    outfile__array_IN_spikemon_count.open(results_dir + "_array_IN_spikemon_count_4137010141", ios::binary | ios::out);
    if(outfile__array_IN_spikemon_count.is_open())
    {
        outfile__array_IN_spikemon_count.write(reinterpret_cast<char*>(_array_IN_spikemon_count), 20*sizeof(_array_IN_spikemon_count[0]));
        outfile__array_IN_spikemon_count.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_spikemon_count." << endl;
    }
    ofstream outfile__array_IN_spikemon_N;
    outfile__array_IN_spikemon_N.open(results_dir + "_array_IN_spikemon_N_3886185682", ios::binary | ios::out);
    if(outfile__array_IN_spikemon_N.is_open())
    {
        outfile__array_IN_spikemon_N.write(reinterpret_cast<char*>(_array_IN_spikemon_N), 1*sizeof(_array_IN_spikemon_N[0]));
        outfile__array_IN_spikemon_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_IN_spikemon_N." << endl;
    }
    ofstream outfile__array_INMDA_PYso_IN_N;
    outfile__array_INMDA_PYso_IN_N.open(results_dir + "_array_INMDA_PYso_IN_N_2668924601", ios::binary | ios::out);
    if(outfile__array_INMDA_PYso_IN_N.is_open())
    {
        outfile__array_INMDA_PYso_IN_N.write(reinterpret_cast<char*>(_array_INMDA_PYso_IN_N), 1*sizeof(_array_INMDA_PYso_IN_N[0]));
        outfile__array_INMDA_PYso_IN_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_INMDA_PYso_IN_N." << endl;
    }
    ofstream outfile__array_INMDA_PYso_IN_sources;
    outfile__array_INMDA_PYso_IN_sources.open(results_dir + "_array_INMDA_PYso_IN_sources_4041238080", ios::binary | ios::out);
    if(outfile__array_INMDA_PYso_IN_sources.is_open())
    {
        outfile__array_INMDA_PYso_IN_sources.write(reinterpret_cast<char*>(_array_INMDA_PYso_IN_sources), 2000*sizeof(_array_INMDA_PYso_IN_sources[0]));
        outfile__array_INMDA_PYso_IN_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_INMDA_PYso_IN_sources." << endl;
    }
    ofstream outfile__array_INMDA_PYso_IN_targets;
    outfile__array_INMDA_PYso_IN_targets.open(results_dir + "_array_INMDA_PYso_IN_targets_2382241185", ios::binary | ios::out);
    if(outfile__array_INMDA_PYso_IN_targets.is_open())
    {
        outfile__array_INMDA_PYso_IN_targets.write(reinterpret_cast<char*>(_array_INMDA_PYso_IN_targets), 2000*sizeof(_array_INMDA_PYso_IN_targets[0]));
        outfile__array_INMDA_PYso_IN_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_INMDA_PYso_IN_targets." << endl;
    }
    ofstream outfile__array_INMDA_PYso_PYdr_N;
    outfile__array_INMDA_PYso_PYdr_N.open(results_dir + "_array_INMDA_PYso_PYdr_N_786129060", ios::binary | ios::out);
    if(outfile__array_INMDA_PYso_PYdr_N.is_open())
    {
        outfile__array_INMDA_PYso_PYdr_N.write(reinterpret_cast<char*>(_array_INMDA_PYso_PYdr_N), 1*sizeof(_array_INMDA_PYso_PYdr_N[0]));
        outfile__array_INMDA_PYso_PYdr_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_INMDA_PYso_PYdr_N." << endl;
    }
    ofstream outfile__array_INMDA_PYso_PYdr_sources;
    outfile__array_INMDA_PYso_PYdr_sources.open(results_dir + "_array_INMDA_PYso_PYdr_sources_855532748", ios::binary | ios::out);
    if(outfile__array_INMDA_PYso_PYdr_sources.is_open())
    {
        outfile__array_INMDA_PYso_PYdr_sources.write(reinterpret_cast<char*>(_array_INMDA_PYso_PYdr_sources), 2000*sizeof(_array_INMDA_PYso_PYdr_sources[0]));
        outfile__array_INMDA_PYso_PYdr_sources.close();
    } else
    {
        std::cout << "Error writing output file for _array_INMDA_PYso_PYdr_sources." << endl;
    }
    ofstream outfile__array_INMDA_PYso_PYdr_targets;
    outfile__array_INMDA_PYso_PYdr_targets.open(results_dir + "_array_INMDA_PYso_PYdr_targets_1340088109", ios::binary | ios::out);
    if(outfile__array_INMDA_PYso_PYdr_targets.is_open())
    {
        outfile__array_INMDA_PYso_PYdr_targets.write(reinterpret_cast<char*>(_array_INMDA_PYso_PYdr_targets), 2000*sizeof(_array_INMDA_PYso_PYdr_targets[0]));
        outfile__array_INMDA_PYso_PYdr_targets.close();
    } else
    {
        std::cout << "Error writing output file for _array_INMDA_PYso_PYdr_targets." << endl;
    }
    ofstream outfile__array_poissongroup_1__spikespace;
    outfile__array_poissongroup_1__spikespace.open(results_dir + "_array_poissongroup_1__spikespace_2558380132", ios::binary | ios::out);
    if(outfile__array_poissongroup_1__spikespace.is_open())
    {
        outfile__array_poissongroup_1__spikespace.write(reinterpret_cast<char*>(_array_poissongroup_1__spikespace), 21*sizeof(_array_poissongroup_1__spikespace[0]));
        outfile__array_poissongroup_1__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_1__spikespace." << endl;
    }
    ofstream outfile__array_poissongroup_1_i;
    outfile__array_poissongroup_1_i.open(results_dir + "_array_poissongroup_1_i_2566510749", ios::binary | ios::out);
    if(outfile__array_poissongroup_1_i.is_open())
    {
        outfile__array_poissongroup_1_i.write(reinterpret_cast<char*>(_array_poissongroup_1_i), 20*sizeof(_array_poissongroup_1_i[0]));
        outfile__array_poissongroup_1_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_1_i." << endl;
    }
    ofstream outfile__array_poissongroup_1_rates;
    outfile__array_poissongroup_1_rates.open(results_dir + "_array_poissongroup_1_rates_3230304882", ios::binary | ios::out);
    if(outfile__array_poissongroup_1_rates.is_open())
    {
        outfile__array_poissongroup_1_rates.write(reinterpret_cast<char*>(_array_poissongroup_1_rates), 20*sizeof(_array_poissongroup_1_rates[0]));
        outfile__array_poissongroup_1_rates.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_1_rates." << endl;
    }
    ofstream outfile__array_poissongroup_2__spikespace;
    outfile__array_poissongroup_2__spikespace.open(results_dir + "_array_poissongroup_2__spikespace_632792234", ios::binary | ios::out);
    if(outfile__array_poissongroup_2__spikespace.is_open())
    {
        outfile__array_poissongroup_2__spikespace.write(reinterpret_cast<char*>(_array_poissongroup_2__spikespace), 21*sizeof(_array_poissongroup_2__spikespace[0]));
        outfile__array_poissongroup_2__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_2__spikespace." << endl;
    }
    ofstream outfile__array_poissongroup_2_i;
    outfile__array_poissongroup_2_i.open(results_dir + "_array_poissongroup_2_i_2596234948", ios::binary | ios::out);
    if(outfile__array_poissongroup_2_i.is_open())
    {
        outfile__array_poissongroup_2_i.write(reinterpret_cast<char*>(_array_poissongroup_2_i), 20*sizeof(_array_poissongroup_2_i[0]));
        outfile__array_poissongroup_2_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_2_i." << endl;
    }
    ofstream outfile__array_poissongroup_2_rates;
    outfile__array_poissongroup_2_rates.open(results_dir + "_array_poissongroup_2_rates_4049768687", ios::binary | ios::out);
    if(outfile__array_poissongroup_2_rates.is_open())
    {
        outfile__array_poissongroup_2_rates.write(reinterpret_cast<char*>(_array_poissongroup_2_rates), 20*sizeof(_array_poissongroup_2_rates[0]));
        outfile__array_poissongroup_2_rates.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_2_rates." << endl;
    }
    ofstream outfile__array_poissongroup__spikespace;
    outfile__array_poissongroup__spikespace.open(results_dir + "_array_poissongroup__spikespace_1019000416", ios::binary | ios::out);
    if(outfile__array_poissongroup__spikespace.is_open())
    {
        outfile__array_poissongroup__spikespace.write(reinterpret_cast<char*>(_array_poissongroup__spikespace), 101*sizeof(_array_poissongroup__spikespace[0]));
        outfile__array_poissongroup__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup__spikespace." << endl;
    }
    ofstream outfile__array_poissongroup_i;
    outfile__array_poissongroup_i.open(results_dir + "_array_poissongroup_i_1277690444", ios::binary | ios::out);
    if(outfile__array_poissongroup_i.is_open())
    {
        outfile__array_poissongroup_i.write(reinterpret_cast<char*>(_array_poissongroup_i), 100*sizeof(_array_poissongroup_i[0]));
        outfile__array_poissongroup_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_i." << endl;
    }
    ofstream outfile__array_poissongroup_rates;
    outfile__array_poissongroup_rates.open(results_dir + "_array_poissongroup_rates_3353413371", ios::binary | ios::out);
    if(outfile__array_poissongroup_rates.is_open())
    {
        outfile__array_poissongroup_rates.write(reinterpret_cast<char*>(_array_poissongroup_rates), 100*sizeof(_array_poissongroup_rates[0]));
        outfile__array_poissongroup_rates.close();
    } else
    {
        std::cout << "Error writing output file for _array_poissongroup_rates." << endl;
    }
    ofstream outfile__array_PYdr_group__spikespace;
    outfile__array_PYdr_group__spikespace.open(results_dir + "_array_PYdr_group__spikespace_2781556690", ios::binary | ios::out);
    if(outfile__array_PYdr_group__spikespace.is_open())
    {
        outfile__array_PYdr_group__spikespace.write(reinterpret_cast<char*>(_array_PYdr_group__spikespace), 101*sizeof(_array_PYdr_group__spikespace[0]));
        outfile__array_PYdr_group__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group__spikespace." << endl;
    }
    ofstream outfile__array_PYdr_group_CaBuffer;
    outfile__array_PYdr_group_CaBuffer.open(results_dir + "_array_PYdr_group_CaBuffer_423175833", ios::binary | ios::out);
    if(outfile__array_PYdr_group_CaBuffer.is_open())
    {
        outfile__array_PYdr_group_CaBuffer.write(reinterpret_cast<char*>(_array_PYdr_group_CaBuffer), 100*sizeof(_array_PYdr_group_CaBuffer[0]));
        outfile__array_PYdr_group_CaBuffer.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_CaBuffer." << endl;
    }
    ofstream outfile__array_PYdr_group_gPoisson_PYdr;
    outfile__array_PYdr_group_gPoisson_PYdr.open(results_dir + "_array_PYdr_group_gPoisson_PYdr_215454924", ios::binary | ios::out);
    if(outfile__array_PYdr_group_gPoisson_PYdr.is_open())
    {
        outfile__array_PYdr_group_gPoisson_PYdr.write(reinterpret_cast<char*>(_array_PYdr_group_gPoisson_PYdr), 100*sizeof(_array_PYdr_group_gPoisson_PYdr[0]));
        outfile__array_PYdr_group_gPoisson_PYdr.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_gPoisson_PYdr." << endl;
    }
    ofstream outfile__array_PYdr_group_i;
    outfile__array_PYdr_group_i.open(results_dir + "_array_PYdr_group_i_287278988", ios::binary | ios::out);
    if(outfile__array_PYdr_group_i.is_open())
    {
        outfile__array_PYdr_group_i.write(reinterpret_cast<char*>(_array_PYdr_group_i), 100*sizeof(_array_PYdr_group_i[0]));
        outfile__array_PYdr_group_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_i." << endl;
    }
    ofstream outfile__array_PYdr_group_iAMPA_PYso_PYdr;
    outfile__array_PYdr_group_iAMPA_PYso_PYdr.open(results_dir + "_array_PYdr_group_iAMPA_PYso_PYdr_1086704099", ios::binary | ios::out);
    if(outfile__array_PYdr_group_iAMPA_PYso_PYdr.is_open())
    {
        outfile__array_PYdr_group_iAMPA_PYso_PYdr.write(reinterpret_cast<char*>(_array_PYdr_group_iAMPA_PYso_PYdr), 100*sizeof(_array_PYdr_group_iAMPA_PYso_PYdr[0]));
        outfile__array_PYdr_group_iAMPA_PYso_PYdr.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_iAMPA_PYso_PYdr." << endl;
    }
    ofstream outfile__array_PYdr_group_iAMPA_TC_PYdr;
    outfile__array_PYdr_group_iAMPA_TC_PYdr.open(results_dir + "_array_PYdr_group_iAMPA_TC_PYdr_2496847366", ios::binary | ios::out);
    if(outfile__array_PYdr_group_iAMPA_TC_PYdr.is_open())
    {
        outfile__array_PYdr_group_iAMPA_TC_PYdr.write(reinterpret_cast<char*>(_array_PYdr_group_iAMPA_TC_PYdr), 100*sizeof(_array_PYdr_group_iAMPA_TC_PYdr[0]));
        outfile__array_PYdr_group_iAMPA_TC_PYdr.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_iAMPA_TC_PYdr." << endl;
    }
    ofstream outfile__array_PYdr_group_iCOM;
    outfile__array_PYdr_group_iCOM.open(results_dir + "_array_PYdr_group_iCOM_3824444896", ios::binary | ios::out);
    if(outfile__array_PYdr_group_iCOM.is_open())
    {
        outfile__array_PYdr_group_iCOM.write(reinterpret_cast<char*>(_array_PYdr_group_iCOM), 100*sizeof(_array_PYdr_group_iCOM[0]));
        outfile__array_PYdr_group_iCOM.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_iCOM." << endl;
    }
    ofstream outfile__array_PYdr_group_iNMDA_PYso_PYdr;
    outfile__array_PYdr_group_iNMDA_PYso_PYdr.open(results_dir + "_array_PYdr_group_iNMDA_PYso_PYdr_424651670", ios::binary | ios::out);
    if(outfile__array_PYdr_group_iNMDA_PYso_PYdr.is_open())
    {
        outfile__array_PYdr_group_iNMDA_PYso_PYdr.write(reinterpret_cast<char*>(_array_PYdr_group_iNMDA_PYso_PYdr), 100*sizeof(_array_PYdr_group_iNMDA_PYso_PYdr[0]));
        outfile__array_PYdr_group_iNMDA_PYso_PYdr.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_iNMDA_PYso_PYdr." << endl;
    }
    ofstream outfile__array_PYdr_group_lastspike;
    outfile__array_PYdr_group_lastspike.open(results_dir + "_array_PYdr_group_lastspike_2977437792", ios::binary | ios::out);
    if(outfile__array_PYdr_group_lastspike.is_open())
    {
        outfile__array_PYdr_group_lastspike.write(reinterpret_cast<char*>(_array_PYdr_group_lastspike), 100*sizeof(_array_PYdr_group_lastspike[0]));
        outfile__array_PYdr_group_lastspike.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_lastspike." << endl;
    }
    ofstream outfile__array_PYdr_group_not_refractory;
    outfile__array_PYdr_group_not_refractory.open(results_dir + "_array_PYdr_group_not_refractory_1406857520", ios::binary | ios::out);
    if(outfile__array_PYdr_group_not_refractory.is_open())
    {
        outfile__array_PYdr_group_not_refractory.write(reinterpret_cast<char*>(_array_PYdr_group_not_refractory), 100*sizeof(_array_PYdr_group_not_refractory[0]));
        outfile__array_PYdr_group_not_refractory.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_not_refractory." << endl;
    }
    ofstream outfile__array_PYdr_group_v;
    outfile__array_PYdr_group_v.open(results_dir + "_array_PYdr_group_v_2618788473", ios::binary | ios::out);
    if(outfile__array_PYdr_group_v.is_open())
    {
        outfile__array_PYdr_group_v.write(reinterpret_cast<char*>(_array_PYdr_group_v), 100*sizeof(_array_PYdr_group_v[0]));
        outfile__array_PYdr_group_v.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_group_v." << endl;
    }
    ofstream outfile__array_PYdr_spikemon__source_idx;
    outfile__array_PYdr_spikemon__source_idx.open(results_dir + "_array_PYdr_spikemon__source_idx_1381534737", ios::binary | ios::out);
    if(outfile__array_PYdr_spikemon__source_idx.is_open())
    {
        outfile__array_PYdr_spikemon__source_idx.write(reinterpret_cast<char*>(_array_PYdr_spikemon__source_idx), 100*sizeof(_array_PYdr_spikemon__source_idx[0]));
        outfile__array_PYdr_spikemon__source_idx.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_spikemon__source_idx." << endl;
    }
    ofstream outfile__array_PYdr_spikemon_count;
    outfile__array_PYdr_spikemon_count.open(results_dir + "_array_PYdr_spikemon_count_3036713658", ios::binary | ios::out);
    if(outfile__array_PYdr_spikemon_count.is_open())
    {
        outfile__array_PYdr_spikemon_count.write(reinterpret_cast<char*>(_array_PYdr_spikemon_count), 100*sizeof(_array_PYdr_spikemon_count[0]));
        outfile__array_PYdr_spikemon_count.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_spikemon_count." << endl;
    }
    ofstream outfile__array_PYdr_spikemon_N;
    outfile__array_PYdr_spikemon_N.open(results_dir + "_array_PYdr_spikemon_N_1543910314", ios::binary | ios::out);
    if(outfile__array_PYdr_spikemon_N.is_open())
    {
        outfile__array_PYdr_spikemon_N.write(reinterpret_cast<char*>(_array_PYdr_spikemon_N), 1*sizeof(_array_PYdr_spikemon_N[0]));
        outfile__array_PYdr_spikemon_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYdr_spikemon_N." << endl;
    }
    ofstream outfile__array_PYso_group__spikespace;
    outfile__array_PYso_group__spikespace.open(results_dir + "_array_PYso_group__spikespace_2149615953", ios::binary | ios::out);
    if(outfile__array_PYso_group__spikespace.is_open())
    {
        outfile__array_PYso_group__spikespace.write(reinterpret_cast<char*>(_array_PYso_group__spikespace), 101*sizeof(_array_PYso_group__spikespace[0]));
        outfile__array_PYso_group__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group__spikespace." << endl;
    }
    ofstream outfile__array_PYso_group_hA;
    outfile__array_PYso_group_hA.open(results_dir + "_array_PYso_group_hA_2067833942", ios::binary | ios::out);
    if(outfile__array_PYso_group_hA.is_open())
    {
        outfile__array_PYso_group_hA.write(reinterpret_cast<char*>(_array_PYso_group_hA), 100*sizeof(_array_PYso_group_hA[0]));
        outfile__array_PYso_group_hA.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_hA." << endl;
    }
    ofstream outfile__array_PYso_group_hNa;
    outfile__array_PYso_group_hNa.open(results_dir + "_array_PYso_group_hNa_3982235369", ios::binary | ios::out);
    if(outfile__array_PYso_group_hNa.is_open())
    {
        outfile__array_PYso_group_hNa.write(reinterpret_cast<char*>(_array_PYso_group_hNa), 100*sizeof(_array_PYso_group_hNa[0]));
        outfile__array_PYso_group_hNa.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_hNa." << endl;
    }
    ofstream outfile__array_PYso_group_i;
    outfile__array_PYso_group_i.open(results_dir + "_array_PYso_group_i_2313336891", ios::binary | ios::out);
    if(outfile__array_PYso_group_i.is_open())
    {
        outfile__array_PYso_group_i.write(reinterpret_cast<char*>(_array_PYso_group_i), 100*sizeof(_array_PYso_group_i[0]));
        outfile__array_PYso_group_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_i." << endl;
    }
    ofstream outfile__array_PYso_group_iCOM;
    outfile__array_PYso_group_iCOM.open(results_dir + "_array_PYso_group_iCOM_3090622832", ios::binary | ios::out);
    if(outfile__array_PYso_group_iCOM.is_open())
    {
        outfile__array_PYso_group_iCOM.write(reinterpret_cast<char*>(_array_PYso_group_iCOM), 100*sizeof(_array_PYso_group_iCOM[0]));
        outfile__array_PYso_group_iCOM.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_iCOM." << endl;
    }
    ofstream outfile__array_PYso_group_iGABAA_IN_PYso;
    outfile__array_PYso_group_iGABAA_IN_PYso.open(results_dir + "_array_PYso_group_iGABAA_IN_PYso_102308486", ios::binary | ios::out);
    if(outfile__array_PYso_group_iGABAA_IN_PYso.is_open())
    {
        outfile__array_PYso_group_iGABAA_IN_PYso.write(reinterpret_cast<char*>(_array_PYso_group_iGABAA_IN_PYso), 100*sizeof(_array_PYso_group_iGABAA_IN_PYso[0]));
        outfile__array_PYso_group_iGABAA_IN_PYso.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_iGABAA_IN_PYso." << endl;
    }
    ofstream outfile__array_PYso_group_iKNa;
    outfile__array_PYso_group_iKNa.open(results_dir + "_array_PYso_group_iKNa_2649171306", ios::binary | ios::out);
    if(outfile__array_PYso_group_iKNa.is_open())
    {
        outfile__array_PYso_group_iKNa.write(reinterpret_cast<char*>(_array_PYso_group_iKNa), 100*sizeof(_array_PYso_group_iKNa[0]));
        outfile__array_PYso_group_iKNa.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_iKNa." << endl;
    }
    ofstream outfile__array_PYso_group_lastspike;
    outfile__array_PYso_group_lastspike.open(results_dir + "_array_PYso_group_lastspike_2065331315", ios::binary | ios::out);
    if(outfile__array_PYso_group_lastspike.is_open())
    {
        outfile__array_PYso_group_lastspike.write(reinterpret_cast<char*>(_array_PYso_group_lastspike), 100*sizeof(_array_PYso_group_lastspike[0]));
        outfile__array_PYso_group_lastspike.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_lastspike." << endl;
    }
    ofstream outfile__array_PYso_group_mKS;
    outfile__array_PYso_group_mKS.open(results_dir + "_array_PYso_group_mKS_1580691911", ios::binary | ios::out);
    if(outfile__array_PYso_group_mKS.is_open())
    {
        outfile__array_PYso_group_mKS.write(reinterpret_cast<char*>(_array_PYso_group_mKS), 100*sizeof(_array_PYso_group_mKS[0]));
        outfile__array_PYso_group_mKS.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_mKS." << endl;
    }
    ofstream outfile__array_PYso_group_nK;
    outfile__array_PYso_group_nK.open(results_dir + "_array_PYso_group_nK_3452955342", ios::binary | ios::out);
    if(outfile__array_PYso_group_nK.is_open())
    {
        outfile__array_PYso_group_nK.write(reinterpret_cast<char*>(_array_PYso_group_nK), 100*sizeof(_array_PYso_group_nK[0]));
        outfile__array_PYso_group_nK.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_nK." << endl;
    }
    ofstream outfile__array_PYso_group_not_refractory;
    outfile__array_PYso_group_not_refractory.open(results_dir + "_array_PYso_group_not_refractory_2133619795", ios::binary | ios::out);
    if(outfile__array_PYso_group_not_refractory.is_open())
    {
        outfile__array_PYso_group_not_refractory.write(reinterpret_cast<char*>(_array_PYso_group_not_refractory), 100*sizeof(_array_PYso_group_not_refractory[0]));
        outfile__array_PYso_group_not_refractory.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_not_refractory." << endl;
    }
    ofstream outfile__array_PYso_group_v;
    outfile__array_PYso_group_v.open(results_dir + "_array_PYso_group_v_82490830", ios::binary | ios::out);
    if(outfile__array_PYso_group_v.is_open())
    {
        outfile__array_PYso_group_v.write(reinterpret_cast<char*>(_array_PYso_group_v), 100*sizeof(_array_PYso_group_v[0]));
        outfile__array_PYso_group_v.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_group_v." << endl;
    }
    ofstream outfile__array_PYso_spikemon__source_idx;
    outfile__array_PYso_spikemon__source_idx.open(results_dir + "_array_PYso_spikemon__source_idx_2125334898", ios::binary | ios::out);
    if(outfile__array_PYso_spikemon__source_idx.is_open())
    {
        outfile__array_PYso_spikemon__source_idx.write(reinterpret_cast<char*>(_array_PYso_spikemon__source_idx), 100*sizeof(_array_PYso_spikemon__source_idx[0]));
        outfile__array_PYso_spikemon__source_idx.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_spikemon__source_idx." << endl;
    }
    ofstream outfile__array_PYso_spikemon_count;
    outfile__array_PYso_spikemon_count.open(results_dir + "_array_PYso_spikemon_count_1840122699", ios::binary | ios::out);
    if(outfile__array_PYso_spikemon_count.is_open())
    {
        outfile__array_PYso_spikemon_count.write(reinterpret_cast<char*>(_array_PYso_spikemon_count), 100*sizeof(_array_PYso_spikemon_count[0]));
        outfile__array_PYso_spikemon_count.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_spikemon_count." << endl;
    }
    ofstream outfile__array_PYso_spikemon_N;
    outfile__array_PYso_spikemon_N.open(results_dir + "_array_PYso_spikemon_N_130383674", ios::binary | ios::out);
    if(outfile__array_PYso_spikemon_N.is_open())
    {
        outfile__array_PYso_spikemon_N.write(reinterpret_cast<char*>(_array_PYso_spikemon_N), 1*sizeof(_array_PYso_spikemon_N[0]));
        outfile__array_PYso_spikemon_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_PYso_spikemon_N." << endl;
    }
    ofstream outfile__array_RE_group__spikespace;
    outfile__array_RE_group__spikespace.open(results_dir + "_array_RE_group__spikespace_3017101979", ios::binary | ios::out);
    if(outfile__array_RE_group__spikespace.is_open())
    {
        outfile__array_RE_group__spikespace.write(reinterpret_cast<char*>(_array_RE_group__spikespace), 21*sizeof(_array_RE_group__spikespace[0]));
        outfile__array_RE_group__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group__spikespace." << endl;
    }
    ofstream outfile__array_RE_group_gPoisson_RE;
    outfile__array_RE_group_gPoisson_RE.open(results_dir + "_array_RE_group_gPoisson_RE_2982234961", ios::binary | ios::out);
    if(outfile__array_RE_group_gPoisson_RE.is_open())
    {
        outfile__array_RE_group_gPoisson_RE.write(reinterpret_cast<char*>(_array_RE_group_gPoisson_RE), 20*sizeof(_array_RE_group_gPoisson_RE[0]));
        outfile__array_RE_group_gPoisson_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_gPoisson_RE." << endl;
    }
    ofstream outfile__array_RE_group_hNa_RE;
    outfile__array_RE_group_hNa_RE.open(results_dir + "_array_RE_group_hNa_RE_351056889", ios::binary | ios::out);
    if(outfile__array_RE_group_hNa_RE.is_open())
    {
        outfile__array_RE_group_hNa_RE.write(reinterpret_cast<char*>(_array_RE_group_hNa_RE), 20*sizeof(_array_RE_group_hNa_RE[0]));
        outfile__array_RE_group_hNa_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_hNa_RE." << endl;
    }
    ofstream outfile__array_RE_group_hT_RE;
    outfile__array_RE_group_hT_RE.open(results_dir + "_array_RE_group_hT_RE_3316204180", ios::binary | ios::out);
    if(outfile__array_RE_group_hT_RE.is_open())
    {
        outfile__array_RE_group_hT_RE.write(reinterpret_cast<char*>(_array_RE_group_hT_RE), 20*sizeof(_array_RE_group_hT_RE[0]));
        outfile__array_RE_group_hT_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_hT_RE." << endl;
    }
    ofstream outfile__array_RE_group_i;
    outfile__array_RE_group_i.open(results_dir + "_array_RE_group_i_1097777293", ios::binary | ios::out);
    if(outfile__array_RE_group_i.is_open())
    {
        outfile__array_RE_group_i.write(reinterpret_cast<char*>(_array_RE_group_i), 20*sizeof(_array_RE_group_i[0]));
        outfile__array_RE_group_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_i." << endl;
    }
    ofstream outfile__array_RE_group_iAMPA_PYso_RE;
    outfile__array_RE_group_iAMPA_PYso_RE.open(results_dir + "_array_RE_group_iAMPA_PYso_RE_1418352129", ios::binary | ios::out);
    if(outfile__array_RE_group_iAMPA_PYso_RE.is_open())
    {
        outfile__array_RE_group_iAMPA_PYso_RE.write(reinterpret_cast<char*>(_array_RE_group_iAMPA_PYso_RE), 20*sizeof(_array_RE_group_iAMPA_PYso_RE[0]));
        outfile__array_RE_group_iAMPA_PYso_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_iAMPA_PYso_RE." << endl;
    }
    ofstream outfile__array_RE_group_iAMPA_TC_RE;
    outfile__array_RE_group_iAMPA_TC_RE.open(results_dir + "_array_RE_group_iAMPA_TC_RE_3597002672", ios::binary | ios::out);
    if(outfile__array_RE_group_iAMPA_TC_RE.is_open())
    {
        outfile__array_RE_group_iAMPA_TC_RE.write(reinterpret_cast<char*>(_array_RE_group_iAMPA_TC_RE), 20*sizeof(_array_RE_group_iAMPA_TC_RE[0]));
        outfile__array_RE_group_iAMPA_TC_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_iAMPA_TC_RE." << endl;
    }
    ofstream outfile__array_RE_group_iGABAA_RE_RE;
    outfile__array_RE_group_iGABAA_RE_RE.open(results_dir + "_array_RE_group_iGABAA_RE_RE_1575339423", ios::binary | ios::out);
    if(outfile__array_RE_group_iGABAA_RE_RE.is_open())
    {
        outfile__array_RE_group_iGABAA_RE_RE.write(reinterpret_cast<char*>(_array_RE_group_iGABAA_RE_RE), 20*sizeof(_array_RE_group_iGABAA_RE_RE[0]));
        outfile__array_RE_group_iGABAA_RE_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_iGABAA_RE_RE." << endl;
    }
    ofstream outfile__array_RE_group_lastspike;
    outfile__array_RE_group_lastspike.open(results_dir + "_array_RE_group_lastspike_68371110", ios::binary | ios::out);
    if(outfile__array_RE_group_lastspike.is_open())
    {
        outfile__array_RE_group_lastspike.write(reinterpret_cast<char*>(_array_RE_group_lastspike), 20*sizeof(_array_RE_group_lastspike[0]));
        outfile__array_RE_group_lastspike.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_lastspike." << endl;
    }
    ofstream outfile__array_RE_group_mNa_RE;
    outfile__array_RE_group_mNa_RE.open(results_dir + "_array_RE_group_mNa_RE_1143022154", ios::binary | ios::out);
    if(outfile__array_RE_group_mNa_RE.is_open())
    {
        outfile__array_RE_group_mNa_RE.write(reinterpret_cast<char*>(_array_RE_group_mNa_RE), 20*sizeof(_array_RE_group_mNa_RE[0]));
        outfile__array_RE_group_mNa_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_mNa_RE." << endl;
    }
    ofstream outfile__array_RE_group_mT_RE;
    outfile__array_RE_group_mT_RE.open(results_dir + "_array_RE_group_mT_RE_222940644", ios::binary | ios::out);
    if(outfile__array_RE_group_mT_RE.is_open())
    {
        outfile__array_RE_group_mT_RE.write(reinterpret_cast<char*>(_array_RE_group_mT_RE), 20*sizeof(_array_RE_group_mT_RE[0]));
        outfile__array_RE_group_mT_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_mT_RE." << endl;
    }
    ofstream outfile__array_RE_group_nK_RE;
    outfile__array_RE_group_nK_RE.open(results_dir + "_array_RE_group_nK_RE_1116991741", ios::binary | ios::out);
    if(outfile__array_RE_group_nK_RE.is_open())
    {
        outfile__array_RE_group_nK_RE.write(reinterpret_cast<char*>(_array_RE_group_nK_RE), 20*sizeof(_array_RE_group_nK_RE[0]));
        outfile__array_RE_group_nK_RE.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_nK_RE." << endl;
    }
    ofstream outfile__array_RE_group_not_refractory;
    outfile__array_RE_group_not_refractory.open(results_dir + "_array_RE_group_not_refractory_1276083633", ios::binary | ios::out);
    if(outfile__array_RE_group_not_refractory.is_open())
    {
        outfile__array_RE_group_not_refractory.write(reinterpret_cast<char*>(_array_RE_group_not_refractory), 20*sizeof(_array_RE_group_not_refractory[0]));
        outfile__array_RE_group_not_refractory.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_not_refractory." << endl;
    }
    ofstream outfile__array_RE_group_v;
    outfile__array_RE_group_v.open(results_dir + "_array_RE_group_v_3429289336", ios::binary | ios::out);
    if(outfile__array_RE_group_v.is_open())
    {
        outfile__array_RE_group_v.write(reinterpret_cast<char*>(_array_RE_group_v), 20*sizeof(_array_RE_group_v[0]));
        outfile__array_RE_group_v.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_group_v." << endl;
    }
    ofstream outfile__array_RE_spikemon__source_idx;
    outfile__array_RE_spikemon__source_idx.open(results_dir + "_array_RE_spikemon__source_idx_1301086352", ios::binary | ios::out);
    if(outfile__array_RE_spikemon__source_idx.is_open())
    {
        outfile__array_RE_spikemon__source_idx.write(reinterpret_cast<char*>(_array_RE_spikemon__source_idx), 20*sizeof(_array_RE_spikemon__source_idx[0]));
        outfile__array_RE_spikemon__source_idx.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_spikemon__source_idx." << endl;
    }
    ofstream outfile__array_RE_spikemon_count;
    outfile__array_RE_spikemon_count.open(results_dir + "_array_RE_spikemon_count_244671751", ios::binary | ios::out);
    if(outfile__array_RE_spikemon_count.is_open())
    {
        outfile__array_RE_spikemon_count.write(reinterpret_cast<char*>(_array_RE_spikemon_count), 20*sizeof(_array_RE_spikemon_count[0]));
        outfile__array_RE_spikemon_count.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_spikemon_count." << endl;
    }
    ofstream outfile__array_RE_spikemon_N;
    outfile__array_RE_spikemon_N.open(results_dir + "_array_RE_spikemon_N_3321807269", ios::binary | ios::out);
    if(outfile__array_RE_spikemon_N.is_open())
    {
        outfile__array_RE_spikemon_N.write(reinterpret_cast<char*>(_array_RE_spikemon_N), 1*sizeof(_array_RE_spikemon_N[0]));
        outfile__array_RE_spikemon_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_RE_spikemon_N." << endl;
    }
    ofstream outfile__array_synapses_1_N;
    outfile__array_synapses_1_N.open(results_dir + "_array_synapses_1_N_1771729519", ios::binary | ios::out);
    if(outfile__array_synapses_1_N.is_open())
    {
        outfile__array_synapses_1_N.write(reinterpret_cast<char*>(_array_synapses_1_N), 1*sizeof(_array_synapses_1_N[0]));
        outfile__array_synapses_1_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_synapses_1_N." << endl;
    }
    ofstream outfile__array_synapses_2_N;
    outfile__array_synapses_2_N.open(results_dir + "_array_synapses_2_N_1809632310", ios::binary | ios::out);
    if(outfile__array_synapses_2_N.is_open())
    {
        outfile__array_synapses_2_N.write(reinterpret_cast<char*>(_array_synapses_2_N), 1*sizeof(_array_synapses_2_N[0]));
        outfile__array_synapses_2_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_synapses_2_N." << endl;
    }
    ofstream outfile__array_synapses_N;
    outfile__array_synapses_N.open(results_dir + "_array_synapses_N_483293785", ios::binary | ios::out);
    if(outfile__array_synapses_N.is_open())
    {
        outfile__array_synapses_N.write(reinterpret_cast<char*>(_array_synapses_N), 1*sizeof(_array_synapses_N[0]));
        outfile__array_synapses_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_synapses_N." << endl;
    }
    ofstream outfile__array_TC_group__spikespace;
    outfile__array_TC_group__spikespace.open(results_dir + "_array_TC_group__spikespace_2854929501", ios::binary | ios::out);
    if(outfile__array_TC_group__spikespace.is_open())
    {
        outfile__array_TC_group__spikespace.write(reinterpret_cast<char*>(_array_TC_group__spikespace), 21*sizeof(_array_TC_group__spikespace[0]));
        outfile__array_TC_group__spikespace.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group__spikespace." << endl;
    }
    ofstream outfile__array_TC_group_Ca_TC;
    outfile__array_TC_group_Ca_TC.open(results_dir + "_array_TC_group_Ca_TC_1397998282", ios::binary | ios::out);
    if(outfile__array_TC_group_Ca_TC.is_open())
    {
        outfile__array_TC_group_Ca_TC.write(reinterpret_cast<char*>(_array_TC_group_Ca_TC), 20*sizeof(_array_TC_group_Ca_TC[0]));
        outfile__array_TC_group_Ca_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_Ca_TC." << endl;
    }
    ofstream outfile__array_TC_group_dCa_extra;
    outfile__array_TC_group_dCa_extra.open(results_dir + "_array_TC_group_dCa_extra_2457723365", ios::binary | ios::out);
    if(outfile__array_TC_group_dCa_extra.is_open())
    {
        outfile__array_TC_group_dCa_extra.write(reinterpret_cast<char*>(_array_TC_group_dCa_extra), 20*sizeof(_array_TC_group_dCa_extra[0]));
        outfile__array_TC_group_dCa_extra.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_dCa_extra." << endl;
    }
    ofstream outfile__array_TC_group_gPoisson_TC;
    outfile__array_TC_group_gPoisson_TC.open(results_dir + "_array_TC_group_gPoisson_TC_386391844", ios::binary | ios::out);
    if(outfile__array_TC_group_gPoisson_TC.is_open())
    {
        outfile__array_TC_group_gPoisson_TC.write(reinterpret_cast<char*>(_array_TC_group_gPoisson_TC), 20*sizeof(_array_TC_group_gPoisson_TC[0]));
        outfile__array_TC_group_gPoisson_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_gPoisson_TC." << endl;
    }
    ofstream outfile__array_TC_group_hNa_TC;
    outfile__array_TC_group_hNa_TC.open(results_dir + "_array_TC_group_hNa_TC_1911369230", ios::binary | ios::out);
    if(outfile__array_TC_group_hNa_TC.is_open())
    {
        outfile__array_TC_group_hNa_TC.write(reinterpret_cast<char*>(_array_TC_group_hNa_TC), 20*sizeof(_array_TC_group_hNa_TC[0]));
        outfile__array_TC_group_hNa_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_hNa_TC." << endl;
    }
    ofstream outfile__array_TC_group_hT_TC;
    outfile__array_TC_group_hT_TC.open(results_dir + "_array_TC_group_hT_TC_582472780", ios::binary | ios::out);
    if(outfile__array_TC_group_hT_TC.is_open())
    {
        outfile__array_TC_group_hT_TC.write(reinterpret_cast<char*>(_array_TC_group_hT_TC), 20*sizeof(_array_TC_group_hT_TC[0]));
        outfile__array_TC_group_hT_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_hT_TC." << endl;
    }
    ofstream outfile__array_TC_group_i;
    outfile__array_TC_group_i.open(results_dir + "_array_TC_group_i_1063962944", ios::binary | ios::out);
    if(outfile__array_TC_group_i.is_open())
    {
        outfile__array_TC_group_i.write(reinterpret_cast<char*>(_array_TC_group_i), 20*sizeof(_array_TC_group_i[0]));
        outfile__array_TC_group_i.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_i." << endl;
    }
    ofstream outfile__array_TC_group_iAMPA_PYso_TC;
    outfile__array_TC_group_iAMPA_PYso_TC.open(results_dir + "_array_TC_group_iAMPA_PYso_TC_1870597108", ios::binary | ios::out);
    if(outfile__array_TC_group_iAMPA_PYso_TC.is_open())
    {
        outfile__array_TC_group_iAMPA_PYso_TC.write(reinterpret_cast<char*>(_array_TC_group_iAMPA_PYso_TC), 20*sizeof(_array_TC_group_iAMPA_PYso_TC[0]));
        outfile__array_TC_group_iAMPA_PYso_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_iAMPA_PYso_TC." << endl;
    }
    ofstream outfile__array_TC_group_iGABAA_RE_TC;
    outfile__array_TC_group_iGABAA_RE_TC.open(results_dir + "_array_TC_group_iGABAA_RE_TC_2428657471", ios::binary | ios::out);
    if(outfile__array_TC_group_iGABAA_RE_TC.is_open())
    {
        outfile__array_TC_group_iGABAA_RE_TC.write(reinterpret_cast<char*>(_array_TC_group_iGABAA_RE_TC), 20*sizeof(_array_TC_group_iGABAA_RE_TC[0]));
        outfile__array_TC_group_iGABAA_RE_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_iGABAA_RE_TC." << endl;
    }
    ofstream outfile__array_TC_group_iGABAB_RE_TC;
    outfile__array_TC_group_iGABAB_RE_TC.open(results_dir + "_array_TC_group_iGABAB_RE_TC_2703902114", ios::binary | ios::out);
    if(outfile__array_TC_group_iGABAB_RE_TC.is_open())
    {
        outfile__array_TC_group_iGABAB_RE_TC.write(reinterpret_cast<char*>(_array_TC_group_iGABAB_RE_TC), 20*sizeof(_array_TC_group_iGABAB_RE_TC[0]));
        outfile__array_TC_group_iGABAB_RE_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_iGABAB_RE_TC." << endl;
    }
    ofstream outfile__array_TC_group_lastspike;
    outfile__array_TC_group_lastspike.open(results_dir + "_array_TC_group_lastspike_2908440875", ios::binary | ios::out);
    if(outfile__array_TC_group_lastspike.is_open())
    {
        outfile__array_TC_group_lastspike.write(reinterpret_cast<char*>(_array_TC_group_lastspike), 20*sizeof(_array_TC_group_lastspike[0]));
        outfile__array_TC_group_lastspike.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_lastspike." << endl;
    }
    ofstream outfile__array_TC_group_mNa_TC;
    outfile__array_TC_group_mNa_TC.open(results_dir + "_array_TC_group_mNa_TC_555794365", ios::binary | ios::out);
    if(outfile__array_TC_group_mNa_TC.is_open())
    {
        outfile__array_TC_group_mNa_TC.write(reinterpret_cast<char*>(_array_TC_group_mNa_TC), 20*sizeof(_array_TC_group_mNa_TC[0]));
        outfile__array_TC_group_mNa_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_mNa_TC." << endl;
    }
    ofstream outfile__array_TC_group_nK_TC;
    outfile__array_TC_group_nK_TC.open(results_dir + "_array_TC_group_nK_TC_2777507365", ios::binary | ios::out);
    if(outfile__array_TC_group_nK_TC.is_open())
    {
        outfile__array_TC_group_nK_TC.write(reinterpret_cast<char*>(_array_TC_group_nK_TC), 20*sizeof(_array_TC_group_nK_TC[0]));
        outfile__array_TC_group_nK_TC.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_nK_TC." << endl;
    }
    ofstream outfile__array_TC_group_not_refractory;
    outfile__array_TC_group_not_refractory.open(results_dir + "_array_TC_group_not_refractory_3543422283", ios::binary | ios::out);
    if(outfile__array_TC_group_not_refractory.is_open())
    {
        outfile__array_TC_group_not_refractory.write(reinterpret_cast<char*>(_array_TC_group_not_refractory), 20*sizeof(_array_TC_group_not_refractory[0]));
        outfile__array_TC_group_not_refractory.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_not_refractory." << endl;
    }
    ofstream outfile__array_TC_group_Open;
    outfile__array_TC_group_Open.open(results_dir + "_array_TC_group_Open_1538063267", ios::binary | ios::out);
    if(outfile__array_TC_group_Open.is_open())
    {
        outfile__array_TC_group_Open.write(reinterpret_cast<char*>(_array_TC_group_Open), 20*sizeof(_array_TC_group_Open[0]));
        outfile__array_TC_group_Open.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_Open." << endl;
    }
    ofstream outfile__array_TC_group_OpenLocked;
    outfile__array_TC_group_OpenLocked.open(results_dir + "_array_TC_group_OpenLocked_709273671", ios::binary | ios::out);
    if(outfile__array_TC_group_OpenLocked.is_open())
    {
        outfile__array_TC_group_OpenLocked.write(reinterpret_cast<char*>(_array_TC_group_OpenLocked), 20*sizeof(_array_TC_group_OpenLocked[0]));
        outfile__array_TC_group_OpenLocked.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_OpenLocked." << endl;
    }
    ofstream outfile__array_TC_group_Pone;
    outfile__array_TC_group_Pone.open(results_dir + "_array_TC_group_Pone_814372964", ios::binary | ios::out);
    if(outfile__array_TC_group_Pone.is_open())
    {
        outfile__array_TC_group_Pone.write(reinterpret_cast<char*>(_array_TC_group_Pone), 20*sizeof(_array_TC_group_Pone[0]));
        outfile__array_TC_group_Pone.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_Pone." << endl;
    }
    ofstream outfile__array_TC_group_v;
    outfile__array_TC_group_v.open(results_dir + "_array_TC_group_v_2992817333", ios::binary | ios::out);
    if(outfile__array_TC_group_v.is_open())
    {
        outfile__array_TC_group_v.write(reinterpret_cast<char*>(_array_TC_group_v), 20*sizeof(_array_TC_group_v[0]));
        outfile__array_TC_group_v.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_group_v." << endl;
    }
    ofstream outfile__array_TC_spikemon__source_idx;
    outfile__array_TC_spikemon__source_idx.open(results_dir + "_array_TC_spikemon__source_idx_3535153258", ios::binary | ios::out);
    if(outfile__array_TC_spikemon__source_idx.is_open())
    {
        outfile__array_TC_spikemon__source_idx.write(reinterpret_cast<char*>(_array_TC_spikemon__source_idx), 20*sizeof(_array_TC_spikemon__source_idx[0]));
        outfile__array_TC_spikemon__source_idx.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_spikemon__source_idx." << endl;
    }
    ofstream outfile__array_TC_spikemon_count;
    outfile__array_TC_spikemon_count.open(results_dir + "_array_TC_spikemon_count_4197497835", ios::binary | ios::out);
    if(outfile__array_TC_spikemon_count.is_open())
    {
        outfile__array_TC_spikemon_count.write(reinterpret_cast<char*>(_array_TC_spikemon_count), 20*sizeof(_array_TC_spikemon_count[0]));
        outfile__array_TC_spikemon_count.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_spikemon_count." << endl;
    }
    ofstream outfile__array_TC_spikemon_N;
    outfile__array_TC_spikemon_N.open(results_dir + "_array_TC_spikemon_N_2317794200", ios::binary | ios::out);
    if(outfile__array_TC_spikemon_N.is_open())
    {
        outfile__array_TC_spikemon_N.write(reinterpret_cast<char*>(_array_TC_spikemon_N), 1*sizeof(_array_TC_spikemon_N[0]));
        outfile__array_TC_spikemon_N.close();
    } else
    {
        std::cout << "Error writing output file for _array_TC_spikemon_N." << endl;
    }

    ofstream outfile__dynamic_array_IAMPA_PYso_IN__synaptic_post;
    outfile__dynamic_array_IAMPA_PYso_IN__synaptic_post.open(results_dir + "_dynamic_array_IAMPA_PYso_IN__synaptic_post_1740671229", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN__synaptic_post.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN__synaptic_post.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN__synaptic_post[0]), _dynamic_array_IAMPA_PYso_IN__synaptic_post.size()*sizeof(_dynamic_array_IAMPA_PYso_IN__synaptic_post[0]));
            outfile__dynamic_array_IAMPA_PYso_IN__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN__synaptic_pre;
    outfile__dynamic_array_IAMPA_PYso_IN__synaptic_pre.open(results_dir + "_dynamic_array_IAMPA_PYso_IN__synaptic_pre_420518020", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN__synaptic_pre[0]), _dynamic_array_IAMPA_PYso_IN__synaptic_pre.size()*sizeof(_dynamic_array_IAMPA_PYso_IN__synaptic_pre[0]));
            outfile__dynamic_array_IAMPA_PYso_IN__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_alpha_AMPA;
    outfile__dynamic_array_IAMPA_PYso_IN_alpha_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_alpha_AMPA_2372597723", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_alpha_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_alpha_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_alpha_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_alpha_AMPA[0]), _dynamic_array_IAMPA_PYso_IN_alpha_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_alpha_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_alpha_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_alpha_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_delay;
    outfile__dynamic_array_IAMPA_PYso_IN_delay.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_delay_2420578933", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_delay.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_delay.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_delay.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_delay[0]), _dynamic_array_IAMPA_PYso_IN_delay.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_delay[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_delay." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_delay_1;
    outfile__dynamic_array_IAMPA_PYso_IN_delay_1.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_delay_1_2182002640", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_delay_1.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_delay_1.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_delay_1[0]), _dynamic_array_IAMPA_PYso_IN_delay_1.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_delay_1[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_delay_1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_delay_1." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_deprFactor;
    outfile__dynamic_array_IAMPA_PYso_IN_deprFactor.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_deprFactor_3329853613", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_deprFactor.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_deprFactor.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_deprFactor.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_deprFactor[0]), _dynamic_array_IAMPA_PYso_IN_deprFactor.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_deprFactor[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_deprFactor.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_deprFactor." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_EAMPA;
    outfile__dynamic_array_IAMPA_PYso_IN_EAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_EAMPA_2489759571", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_EAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_EAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_EAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_EAMPA[0]), _dynamic_array_IAMPA_PYso_IN_EAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_EAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_EAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_EAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_gAMPA;
    outfile__dynamic_array_IAMPA_PYso_IN_gAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_gAMPA_795331383", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_gAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_gAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_gAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_gAMPA[0]), _dynamic_array_IAMPA_PYso_IN_gAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_gAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_gAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_gAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_N_incoming;
    outfile__dynamic_array_IAMPA_PYso_IN_N_incoming.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_N_incoming_2384803100", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_N_incoming.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_N_incoming.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_N_incoming[0]), _dynamic_array_IAMPA_PYso_IN_N_incoming.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_N_incoming[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_N_outgoing;
    outfile__dynamic_array_IAMPA_PYso_IN_N_outgoing.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_N_outgoing_2839102918", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_N_outgoing.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_N_outgoing.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_N_outgoing[0]), _dynamic_array_IAMPA_PYso_IN_N_outgoing.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_N_outgoing[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_Npost;
    outfile__dynamic_array_IAMPA_PYso_IN_Npost.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_Npost_2029140042", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_Npost.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_Npost.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_Npost[0]), _dynamic_array_IAMPA_PYso_IN_Npost.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_Npost[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_Npre;
    outfile__dynamic_array_IAMPA_PYso_IN_Npre.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_Npre_4045619543", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_Npre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_Npre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_Npre[0]), _dynamic_array_IAMPA_PYso_IN_Npre.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_Npre[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_radius;
    outfile__dynamic_array_IAMPA_PYso_IN_radius.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_radius_2095873215", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_radius.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_radius.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_radius.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_radius[0]), _dynamic_array_IAMPA_PYso_IN_radius.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_radius[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_radius." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool;
    outfile__dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool_2404868049", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool[0]), _dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_res_AMPA;
    outfile__dynamic_array_IAMPA_PYso_IN_res_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_res_AMPA_1427517787", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_res_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_res_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_res_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_res_AMPA[0]), _dynamic_array_IAMPA_PYso_IN_res_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_res_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_res_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_res_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_sAMPA;
    outfile__dynamic_array_IAMPA_PYso_IN_sAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_sAMPA_3121085045", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_sAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_sAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_sAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_sAMPA[0]), _dynamic_array_IAMPA_PYso_IN_sAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_sAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_sAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_sAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_spike_activity;
    outfile__dynamic_array_IAMPA_PYso_IN_spike_activity.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_spike_activity_2767746912", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_spike_activity.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_spike_activity.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_spike_activity.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_spike_activity[0]), _dynamic_array_IAMPA_PYso_IN_spike_activity.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_spike_activity[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_spike_activity.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_spike_activity." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_tauAMPA;
    outfile__dynamic_array_IAMPA_PYso_IN_tauAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_tauAMPA_1199313786", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_tauAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_tauAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_tauAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_tauAMPA[0]), _dynamic_array_IAMPA_PYso_IN_tauAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_tauAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_tauAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_tauAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_IN_tauRes_AMPA;
    outfile__dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_IN_tauRes_AMPA_3767562965", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_IN_tauRes_AMPA[0]), _dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_IN_tauRes_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_IN_tauRes_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_IN_tauRes_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_post;
    outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_post.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr__synaptic_post_1744184240", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_post.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr__synaptic_post.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr__synaptic_post[0]), _dynamic_array_IAMPA_PYso_PYdr__synaptic_post.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr__synaptic_post[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_pre;
    outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_pre.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr__synaptic_pre_800839044", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr__synaptic_pre[0]), _dynamic_array_IAMPA_PYso_PYdr__synaptic_pre.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr__synaptic_pre[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA;
    outfile__dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA_337350468", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA[0]), _dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_alpha_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_delay;
    outfile__dynamic_array_IAMPA_PYso_PYdr_delay.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_delay_1952554892", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_delay.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_delay.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_delay.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_delay[0]), _dynamic_array_IAMPA_PYso_PYdr_delay.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_delay[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_delay." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_delay_1;
    outfile__dynamic_array_IAMPA_PYso_PYdr_delay_1.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_delay_1_3668082349", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_delay_1.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_delay_1.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_delay_1[0]), _dynamic_array_IAMPA_PYso_PYdr_delay_1.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_delay_1[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_delay_1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_delay_1." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_deprFactor;
    outfile__dynamic_array_IAMPA_PYso_PYdr_deprFactor.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_deprFactor_1594421298", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_deprFactor.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_deprFactor.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_deprFactor.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_deprFactor[0]), _dynamic_array_IAMPA_PYso_PYdr_deprFactor.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_deprFactor[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_deprFactor.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_deprFactor." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_EAMPA;
    outfile__dynamic_array_IAMPA_PYso_PYdr_EAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_EAMPA_1883255466", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_EAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_EAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_EAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_EAMPA[0]), _dynamic_array_IAMPA_PYso_PYdr_EAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_EAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_EAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_EAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_gAMPA;
    outfile__dynamic_array_IAMPA_PYso_PYdr_gAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_gAMPA_3410054862", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_gAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_gAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_gAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_gAMPA[0]), _dynamic_array_IAMPA_PYso_PYdr_gAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_gAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_gAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_gAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_N_incoming;
    outfile__dynamic_array_IAMPA_PYso_PYdr_N_incoming.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_N_incoming_391400835", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_N_incoming.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_N_incoming.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_N_incoming[0]), _dynamic_array_IAMPA_PYso_PYdr_N_incoming.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_N_incoming[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_N_outgoing;
    outfile__dynamic_array_IAMPA_PYso_PYdr_N_outgoing.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_N_outgoing_810040665", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_N_outgoing.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_N_outgoing.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_N_outgoing[0]), _dynamic_array_IAMPA_PYso_PYdr_N_outgoing.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_N_outgoing[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_Npost;
    outfile__dynamic_array_IAMPA_PYso_PYdr_Npost.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_Npost_2631187891", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_Npost.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_Npost.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_Npost[0]), _dynamic_array_IAMPA_PYso_PYdr_Npost.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_Npost[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_Npre;
    outfile__dynamic_array_IAMPA_PYso_PYdr_Npre.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_Npre_3647583195", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_Npre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_Npre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_Npre[0]), _dynamic_array_IAMPA_PYso_PYdr_Npre.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_Npre[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_radius;
    outfile__dynamic_array_IAMPA_PYso_PYdr_radius.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_radius_3093893250", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_radius.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_radius.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_radius.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_radius[0]), _dynamic_array_IAMPA_PYso_PYdr_radius.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_radius[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_radius." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool;
    outfile__dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool_4147242851", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool[0]), _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
    outfile__dynamic_array_IAMPA_PYso_PYdr_res_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_res_AMPA_2080018451", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_res_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_res_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_res_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_res_AMPA[0]), _dynamic_array_IAMPA_PYso_PYdr_res_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_res_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_res_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_res_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_sAMPA;
    outfile__dynamic_array_IAMPA_PYso_PYdr_sAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_sAMPA_1579253644", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_sAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_sAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_sAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_sAMPA[0]), _dynamic_array_IAMPA_PYso_PYdr_sAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_sAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_sAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_sAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_spike_activity;
    outfile__dynamic_array_IAMPA_PYso_PYdr_spike_activity.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_spike_activity_2765020205", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_spike_activity.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_spike_activity.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_spike_activity.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_spike_activity[0]), _dynamic_array_IAMPA_PYso_PYdr_spike_activity.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_spike_activity[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_spike_activity.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_spike_activity." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_tauAMPA;
    outfile__dynamic_array_IAMPA_PYso_PYdr_tauAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_tauAMPA_533780999", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_tauAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_tauAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_tauAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_tauAMPA[0]), _dynamic_array_IAMPA_PYso_PYdr_tauAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_tauAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_tauAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_tauAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA;
    outfile__dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA_2159644024", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA[0]), _dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_PYdr_tauRes_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE__synaptic_post;
    outfile__dynamic_array_IAMPA_PYso_RE__synaptic_post.open(results_dir + "_dynamic_array_IAMPA_PYso_RE__synaptic_post_2680172071", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE__synaptic_post.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE__synaptic_post.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE__synaptic_post[0]), _dynamic_array_IAMPA_PYso_RE__synaptic_post.size()*sizeof(_dynamic_array_IAMPA_PYso_RE__synaptic_post[0]));
            outfile__dynamic_array_IAMPA_PYso_RE__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE__synaptic_pre;
    outfile__dynamic_array_IAMPA_PYso_RE__synaptic_pre.open(results_dir + "_dynamic_array_IAMPA_PYso_RE__synaptic_pre_2065181529", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE__synaptic_pre[0]), _dynamic_array_IAMPA_PYso_RE__synaptic_pre.size()*sizeof(_dynamic_array_IAMPA_PYso_RE__synaptic_pre[0]));
            outfile__dynamic_array_IAMPA_PYso_RE__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE;
    outfile__dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE_1942244292", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE[0]), _dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_EAMPA_PYso_RE." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE;
    outfile__dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE_2699474298", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE[0]), _dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_gAMPA_PYso_RE." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_N_incoming;
    outfile__dynamic_array_IAMPA_PYso_RE_N_incoming.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_N_incoming_2893645931", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_N_incoming.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_N_incoming.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_N_incoming[0]), _dynamic_array_IAMPA_PYso_RE_N_incoming.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_N_incoming[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_N_outgoing;
    outfile__dynamic_array_IAMPA_PYso_RE_N_outgoing.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_N_outgoing_2338715825", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_N_outgoing.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_N_outgoing.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_N_outgoing[0]), _dynamic_array_IAMPA_PYso_RE_N_outgoing.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_N_outgoing[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_Npost;
    outfile__dynamic_array_IAMPA_PYso_RE_Npost.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_Npost_1529102797", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_Npost.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_Npost.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_Npost[0]), _dynamic_array_IAMPA_PYso_RE_Npost.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_Npost[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_Npre;
    outfile__dynamic_array_IAMPA_PYso_RE_Npre.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_Npre_4265441696", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_Npre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_Npre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_Npre[0]), _dynamic_array_IAMPA_PYso_RE_Npre.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_Npre[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_P_AMPA;
    outfile__dynamic_array_IAMPA_PYso_RE_P_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_P_AMPA_3555339976", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_P_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_P_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_P_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_P_AMPA[0]), _dynamic_array_IAMPA_PYso_RE_P_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_P_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_P_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_P_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_radius;
    outfile__dynamic_array_IAMPA_PYso_RE_radius.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_radius_252951621", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_radius.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_radius.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_radius.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_radius[0]), _dynamic_array_IAMPA_PYso_RE_radius.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_radius[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_radius." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool;
    outfile__dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool_208916321", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool[0]), _dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
    outfile__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE_3562567990", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE[0]), _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_RE_tauAMPA;
    outfile__dynamic_array_IAMPA_PYso_RE_tauAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_RE_tauAMPA_443021240", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_RE_tauAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_RE_tauAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_RE_tauAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_RE_tauAMPA[0]), _dynamic_array_IAMPA_PYso_RE_tauAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_RE_tauAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_RE_tauAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_RE_tauAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC__synaptic_post;
    outfile__dynamic_array_IAMPA_PYso_TC__synaptic_post.open(results_dir + "_dynamic_array_IAMPA_PYso_TC__synaptic_post_1801842891", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC__synaptic_post.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC__synaptic_post.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC__synaptic_post[0]), _dynamic_array_IAMPA_PYso_TC__synaptic_post.size()*sizeof(_dynamic_array_IAMPA_PYso_TC__synaptic_post[0]));
            outfile__dynamic_array_IAMPA_PYso_TC__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC__synaptic_pre;
    outfile__dynamic_array_IAMPA_PYso_TC__synaptic_pre.open(results_dir + "_dynamic_array_IAMPA_PYso_TC__synaptic_pre_175347279", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC__synaptic_pre[0]), _dynamic_array_IAMPA_PYso_TC__synaptic_pre.size()*sizeof(_dynamic_array_IAMPA_PYso_TC__synaptic_pre[0]));
            outfile__dynamic_array_IAMPA_PYso_TC__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_EAMPA;
    outfile__dynamic_array_IAMPA_PYso_TC_EAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_EAMPA_303034985", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_EAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_EAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_EAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_EAMPA[0]), _dynamic_array_IAMPA_PYso_TC_EAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_EAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_EAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_EAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_gAMPA;
    outfile__dynamic_array_IAMPA_PYso_TC_gAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_gAMPA_2836303373", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_gAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_gAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_gAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_gAMPA[0]), _dynamic_array_IAMPA_PYso_TC_gAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_gAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_gAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_gAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_N_incoming;
    outfile__dynamic_array_IAMPA_PYso_TC_N_incoming.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_N_incoming_3819041366", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_N_incoming.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_N_incoming.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_N_incoming[0]), _dynamic_array_IAMPA_PYso_TC_N_incoming.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_N_incoming[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_N_outgoing;
    outfile__dynamic_array_IAMPA_PYso_TC_N_outgoing.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_N_outgoing_3300754060", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_N_outgoing.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_N_outgoing.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_N_outgoing[0]), _dynamic_array_IAMPA_PYso_TC_N_outgoing.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_N_outgoing[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_Npost;
    outfile__dynamic_array_IAMPA_PYso_TC_Npost.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_Npost_4271603056", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_Npost.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_Npost.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_Npost[0]), _dynamic_array_IAMPA_PYso_TC_Npost.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_Npost[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_Npre;
    outfile__dynamic_array_IAMPA_PYso_TC_Npre.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_Npre_1270132615", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_Npre.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_Npre.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_Npre[0]), _dynamic_array_IAMPA_PYso_TC_Npre.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_Npre[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_P_AMPA;
    outfile__dynamic_array_IAMPA_PYso_TC_P_AMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_P_AMPA_1721720517", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_P_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_P_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_P_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_P_AMPA[0]), _dynamic_array_IAMPA_PYso_TC_P_AMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_P_AMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_P_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_P_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_radius;
    outfile__dynamic_array_IAMPA_PYso_TC_radius.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_radius_3127299144", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_radius.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_radius.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_radius.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_radius[0]), _dynamic_array_IAMPA_PYso_TC_radius.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_radius[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_radius." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool;
    outfile__dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool_1367666507", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool[0]), _dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_sAMPA;
    outfile__dynamic_array_IAMPA_PYso_TC_sAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_sAMPA_1013890895", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_sAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_sAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_sAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_sAMPA[0]), _dynamic_array_IAMPA_PYso_TC_sAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_sAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_sAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_sAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_PYso_TC_tauAMPA;
    outfile__dynamic_array_IAMPA_PYso_TC_tauAMPA.open(results_dir + "_dynamic_array_IAMPA_PYso_TC_tauAMPA_1684274805", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_PYso_TC_tauAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_PYso_TC_tauAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_PYso_TC_tauAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_PYso_TC_tauAMPA[0]), _dynamic_array_IAMPA_PYso_TC_tauAMPA.size()*sizeof(_dynamic_array_IAMPA_PYso_TC_tauAMPA[0]));
            outfile__dynamic_array_IAMPA_PYso_TC_tauAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_PYso_TC_tauAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN__synaptic_post;
    outfile__dynamic_array_IAMPA_TC_IN__synaptic_post.open(results_dir + "_dynamic_array_IAMPA_TC_IN__synaptic_post_4157283575", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN__synaptic_post.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN__synaptic_post.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN__synaptic_post[0]), _dynamic_array_IAMPA_TC_IN__synaptic_post.size()*sizeof(_dynamic_array_IAMPA_TC_IN__synaptic_post[0]));
            outfile__dynamic_array_IAMPA_TC_IN__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN__synaptic_pre;
    outfile__dynamic_array_IAMPA_TC_IN__synaptic_pre.open(results_dir + "_dynamic_array_IAMPA_TC_IN__synaptic_pre_2912486283", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN__synaptic_pre[0]), _dynamic_array_IAMPA_TC_IN__synaptic_pre.size()*sizeof(_dynamic_array_IAMPA_TC_IN__synaptic_pre[0]));
            outfile__dynamic_array_IAMPA_TC_IN__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_EAMPA;
    outfile__dynamic_array_IAMPA_TC_IN_EAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_IN_EAMPA_3630334459", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_EAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_EAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_EAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_EAMPA[0]), _dynamic_array_IAMPA_TC_IN_EAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_IN_EAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_IN_EAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_EAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_gAMPA;
    outfile__dynamic_array_IAMPA_TC_IN_gAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_IN_gAMPA_1667497375", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_gAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_gAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_gAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_gAMPA[0]), _dynamic_array_IAMPA_TC_IN_gAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_IN_gAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_IN_gAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_gAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_N_incoming;
    outfile__dynamic_array_IAMPA_TC_IN_N_incoming.open(results_dir + "_dynamic_array_IAMPA_TC_IN_N_incoming_1089990371", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_N_incoming.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_N_incoming.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_N_incoming[0]), _dynamic_array_IAMPA_TC_IN_N_incoming.size()*sizeof(_dynamic_array_IAMPA_TC_IN_N_incoming[0]));
            outfile__dynamic_array_IAMPA_TC_IN_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_N_outgoing;
    outfile__dynamic_array_IAMPA_TC_IN_N_outgoing.open(results_dir + "_dynamic_array_IAMPA_TC_IN_N_outgoing_1743494713", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_N_outgoing.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_N_outgoing.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_N_outgoing[0]), _dynamic_array_IAMPA_TC_IN_N_outgoing.size()*sizeof(_dynamic_array_IAMPA_TC_IN_N_outgoing[0]));
            outfile__dynamic_array_IAMPA_TC_IN_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_Npost;
    outfile__dynamic_array_IAMPA_TC_IN_Npost.open(results_dir + "_dynamic_array_IAMPA_TC_IN_Npost_888540898", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_Npost.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_Npost.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_Npost[0]), _dynamic_array_IAMPA_TC_IN_Npost.size()*sizeof(_dynamic_array_IAMPA_TC_IN_Npost[0]));
            outfile__dynamic_array_IAMPA_TC_IN_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_Npre;
    outfile__dynamic_array_IAMPA_TC_IN_Npre.open(results_dir + "_dynamic_array_IAMPA_TC_IN_Npre_2617620342", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_Npre.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_Npre.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_Npre[0]), _dynamic_array_IAMPA_TC_IN_Npre.size()*sizeof(_dynamic_array_IAMPA_TC_IN_Npre[0]));
            outfile__dynamic_array_IAMPA_TC_IN_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_P_AMPA;
    outfile__dynamic_array_IAMPA_TC_IN_P_AMPA.open(results_dir + "_dynamic_array_IAMPA_TC_IN_P_AMPA_2018820574", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_P_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_P_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_P_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_P_AMPA[0]), _dynamic_array_IAMPA_TC_IN_P_AMPA.size()*sizeof(_dynamic_array_IAMPA_TC_IN_P_AMPA[0]));
            outfile__dynamic_array_IAMPA_TC_IN_P_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_P_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_radius;
    outfile__dynamic_array_IAMPA_TC_IN_radius.open(results_dir + "_dynamic_array_IAMPA_TC_IN_radius_2762822483", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_radius.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_radius.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_radius.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_radius[0]), _dynamic_array_IAMPA_TC_IN_radius.size()*sizeof(_dynamic_array_IAMPA_TC_IN_radius[0]));
            outfile__dynamic_array_IAMPA_TC_IN_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_radius." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_remove_recurrent_bool;
    outfile__dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.open(results_dir + "_dynamic_array_IAMPA_TC_IN_remove_recurrent_bool_1732281471", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_remove_recurrent_bool[0]), _dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.size()*sizeof(_dynamic_array_IAMPA_TC_IN_remove_recurrent_bool[0]));
            outfile__dynamic_array_IAMPA_TC_IN_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_sAMPA;
    outfile__dynamic_array_IAMPA_TC_IN_sAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_IN_sAMPA_4127444189", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_sAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_sAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_sAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_sAMPA[0]), _dynamic_array_IAMPA_TC_IN_sAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_IN_sAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_IN_sAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_sAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_IN_tauAMPA;
    outfile__dynamic_array_IAMPA_TC_IN_tauAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_IN_tauAMPA_3994614790", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_IN_tauAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_IN_tauAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_IN_tauAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_IN_tauAMPA[0]), _dynamic_array_IAMPA_TC_IN_tauAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_IN_tauAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_IN_tauAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_IN_tauAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_post;
    outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_post.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr__synaptic_post_1786449004", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_post.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr__synaptic_post.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr__synaptic_post[0]), _dynamic_array_IAMPA_TC_PYdr__synaptic_post.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr__synaptic_post[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_pre;
    outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_pre.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr__synaptic_pre_3489543950", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr__synaptic_pre[0]), _dynamic_array_IAMPA_TC_PYdr__synaptic_pre.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr__synaptic_pre[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_EAMPA;
    outfile__dynamic_array_IAMPA_TC_PYdr_EAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_EAMPA_3643071958", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_EAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_EAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_EAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_EAMPA[0]), _dynamic_array_IAMPA_TC_PYdr_EAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_EAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_EAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_EAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_gAMPA;
    outfile__dynamic_array_IAMPA_TC_PYdr_gAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_gAMPA_1646633394", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_gAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_gAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_gAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_gAMPA[0]), _dynamic_array_IAMPA_TC_PYdr_gAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_gAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_gAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_gAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_N_incoming;
    outfile__dynamic_array_IAMPA_TC_PYdr_N_incoming.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_N_incoming_2507240831", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_N_incoming.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_N_incoming.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_N_incoming[0]), _dynamic_array_IAMPA_TC_PYdr_N_incoming.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_N_incoming[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_N_outgoing;
    outfile__dynamic_array_IAMPA_TC_PYdr_N_outgoing.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_N_outgoing_2993489317", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_N_outgoing.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_N_outgoing.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_N_outgoing[0]), _dynamic_array_IAMPA_TC_PYdr_N_outgoing.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_N_outgoing[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_Npost;
    outfile__dynamic_array_IAMPA_TC_PYdr_Npost.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_Npost_900747983", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_Npost.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_Npost.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_Npost[0]), _dynamic_array_IAMPA_TC_PYdr_Npost.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_Npost[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_Npre;
    outfile__dynamic_array_IAMPA_TC_PYdr_Npre.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_Npre_19438647", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_Npre.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_Npre.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_Npre[0]), _dynamic_array_IAMPA_TC_PYdr_Npre.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_Npre[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_P_AMPA;
    outfile__dynamic_array_IAMPA_TC_PYdr_P_AMPA.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_P_AMPA_1032510455", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_P_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_P_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_P_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_P_AMPA[0]), _dynamic_array_IAMPA_TC_PYdr_P_AMPA.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_P_AMPA[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_P_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_P_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_radius;
    outfile__dynamic_array_IAMPA_TC_PYdr_radius.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_radius_3782429050", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_radius.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_radius.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_radius.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_radius[0]), _dynamic_array_IAMPA_TC_PYdr_radius.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_radius[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_radius." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool;
    outfile__dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool_3602158927", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool[0]), _dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_sAMPA;
    outfile__dynamic_array_IAMPA_TC_PYdr_sAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_sAMPA_4148531440", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_sAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_sAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_sAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_sAMPA[0]), _dynamic_array_IAMPA_TC_PYdr_sAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_sAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_sAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_sAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_PYdr_tauAMPA;
    outfile__dynamic_array_IAMPA_TC_PYdr_tauAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_PYdr_tauAMPA_2901391984", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_PYdr_tauAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_PYdr_tauAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_PYdr_tauAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_PYdr_tauAMPA[0]), _dynamic_array_IAMPA_TC_PYdr_tauAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_PYdr_tauAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_PYdr_tauAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_PYdr_tauAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE__synaptic_post;
    outfile__dynamic_array_IAMPA_TC_RE__synaptic_post.open(results_dir + "_dynamic_array_IAMPA_TC_RE__synaptic_post_265009709", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE__synaptic_post.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE__synaptic_post.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE__synaptic_post[0]), _dynamic_array_IAMPA_TC_RE__synaptic_post.size()*sizeof(_dynamic_array_IAMPA_TC_RE__synaptic_post[0]));
            outfile__dynamic_array_IAMPA_TC_RE__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE__synaptic_pre;
    outfile__dynamic_array_IAMPA_TC_RE__synaptic_pre.open(results_dir + "_dynamic_array_IAMPA_TC_RE__synaptic_pre_3482431574", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE__synaptic_pre[0]), _dynamic_array_IAMPA_TC_RE__synaptic_pre.size()*sizeof(_dynamic_array_IAMPA_TC_RE__synaptic_pre[0]));
            outfile__dynamic_array_IAMPA_TC_RE__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_EAMPA;
    outfile__dynamic_array_IAMPA_TC_RE_EAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_RE_EAMPA_4222941308", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_EAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_EAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_EAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_EAMPA[0]), _dynamic_array_IAMPA_TC_RE_EAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_RE_EAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_RE_EAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_EAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_gAMPA;
    outfile__dynamic_array_IAMPA_TC_RE_gAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_RE_gAMPA_1085637656", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_gAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_gAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_gAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_gAMPA[0]), _dynamic_array_IAMPA_TC_RE_gAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_RE_gAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_RE_gAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_gAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_N_incoming;
    outfile__dynamic_array_IAMPA_TC_RE_N_incoming.open(results_dir + "_dynamic_array_IAMPA_TC_RE_N_incoming_1655395220", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_N_incoming.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_N_incoming.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_N_incoming[0]), _dynamic_array_IAMPA_TC_RE_N_incoming.size()*sizeof(_dynamic_array_IAMPA_TC_RE_N_incoming[0]));
            outfile__dynamic_array_IAMPA_TC_RE_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_N_outgoing;
    outfile__dynamic_array_IAMPA_TC_RE_N_outgoing.open(results_dir + "_dynamic_array_IAMPA_TC_RE_N_outgoing_1169630030", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_N_outgoing.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_N_outgoing.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_N_outgoing[0]), _dynamic_array_IAMPA_TC_RE_N_outgoing.size()*sizeof(_dynamic_array_IAMPA_TC_RE_N_outgoing[0]));
            outfile__dynamic_array_IAMPA_TC_RE_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_Npost;
    outfile__dynamic_array_IAMPA_TC_RE_Npost.open(results_dir + "_dynamic_array_IAMPA_TC_RE_Npost_388003685", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_Npost.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_Npost.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_Npost[0]), _dynamic_array_IAMPA_TC_RE_Npost.size()*sizeof(_dynamic_array_IAMPA_TC_RE_Npost[0]));
            outfile__dynamic_array_IAMPA_TC_RE_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_Npre;
    outfile__dynamic_array_IAMPA_TC_RE_Npre.open(results_dir + "_dynamic_array_IAMPA_TC_RE_Npre_2468056961", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_Npre.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_Npre.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_Npre[0]), _dynamic_array_IAMPA_TC_RE_Npre.size()*sizeof(_dynamic_array_IAMPA_TC_RE_Npre[0]));
            outfile__dynamic_array_IAMPA_TC_RE_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_P_AMPA;
    outfile__dynamic_array_IAMPA_TC_RE_P_AMPA.open(results_dir + "_dynamic_array_IAMPA_TC_RE_P_AMPA_195756324", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_P_AMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_P_AMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_P_AMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_P_AMPA[0]), _dynamic_array_IAMPA_TC_RE_P_AMPA.size()*sizeof(_dynamic_array_IAMPA_TC_RE_P_AMPA[0]));
            outfile__dynamic_array_IAMPA_TC_RE_P_AMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_P_AMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_sAMPA;
    outfile__dynamic_array_IAMPA_TC_RE_sAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_RE_sAMPA_3587552602", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_sAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_sAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_sAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_sAMPA[0]), _dynamic_array_IAMPA_TC_RE_sAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_RE_sAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_RE_sAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_sAMPA." << endl;
    }
    ofstream outfile__dynamic_array_IAMPA_TC_RE_tauAMPA;
    outfile__dynamic_array_IAMPA_TC_RE_tauAMPA.open(results_dir + "_dynamic_array_IAMPA_TC_RE_tauAMPA_3003324612", ios::binary | ios::out);
    if(outfile__dynamic_array_IAMPA_TC_RE_tauAMPA.is_open())
    {
        if (! _dynamic_array_IAMPA_TC_RE_tauAMPA.empty() )
        {
            outfile__dynamic_array_IAMPA_TC_RE_tauAMPA.write(reinterpret_cast<char*>(&_dynamic_array_IAMPA_TC_RE_tauAMPA[0]), _dynamic_array_IAMPA_TC_RE_tauAMPA.size()*sizeof(_dynamic_array_IAMPA_TC_RE_tauAMPA[0]));
            outfile__dynamic_array_IAMPA_TC_RE_tauAMPA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IAMPA_TC_RE_tauAMPA." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_post;
    outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_post.open(results_dir + "_dynamic_array_ICOM_PYdr_PYso__synaptic_post_1315199411", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_post.is_open())
    {
        if (! _dynamic_array_ICOM_PYdr_PYso__synaptic_post.empty() )
        {
            outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYdr_PYso__synaptic_post[0]), _dynamic_array_ICOM_PYdr_PYso__synaptic_post.size()*sizeof(_dynamic_array_ICOM_PYdr_PYso__synaptic_post[0]));
            outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYdr_PYso__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_pre;
    outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_pre.open(results_dir + "_dynamic_array_ICOM_PYdr_PYso__synaptic_pre_1678266109", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_pre.is_open())
    {
        if (! _dynamic_array_ICOM_PYdr_PYso__synaptic_pre.empty() )
        {
            outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYdr_PYso__synaptic_pre[0]), _dynamic_array_ICOM_PYdr_PYso__synaptic_pre.size()*sizeof(_dynamic_array_ICOM_PYdr_PYso__synaptic_pre[0]));
            outfile__dynamic_array_ICOM_PYdr_PYso__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYdr_PYso__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYdr_PYso_gCOM;
    outfile__dynamic_array_ICOM_PYdr_PYso_gCOM.open(results_dir + "_dynamic_array_ICOM_PYdr_PYso_gCOM_4273950301", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYdr_PYso_gCOM.is_open())
    {
        if (! _dynamic_array_ICOM_PYdr_PYso_gCOM.empty() )
        {
            outfile__dynamic_array_ICOM_PYdr_PYso_gCOM.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYdr_PYso_gCOM[0]), _dynamic_array_ICOM_PYdr_PYso_gCOM.size()*sizeof(_dynamic_array_ICOM_PYdr_PYso_gCOM[0]));
            outfile__dynamic_array_ICOM_PYdr_PYso_gCOM.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYdr_PYso_gCOM." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYdr_PYso_N_incoming;
    outfile__dynamic_array_ICOM_PYdr_PYso_N_incoming.open(results_dir + "_dynamic_array_ICOM_PYdr_PYso_N_incoming_390757840", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYdr_PYso_N_incoming.is_open())
    {
        if (! _dynamic_array_ICOM_PYdr_PYso_N_incoming.empty() )
        {
            outfile__dynamic_array_ICOM_PYdr_PYso_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYdr_PYso_N_incoming[0]), _dynamic_array_ICOM_PYdr_PYso_N_incoming.size()*sizeof(_dynamic_array_ICOM_PYdr_PYso_N_incoming[0]));
            outfile__dynamic_array_ICOM_PYdr_PYso_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYdr_PYso_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYdr_PYso_N_outgoing;
    outfile__dynamic_array_ICOM_PYdr_PYso_N_outgoing.open(results_dir + "_dynamic_array_ICOM_PYdr_PYso_N_outgoing_810945802", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYdr_PYso_N_outgoing.is_open())
    {
        if (! _dynamic_array_ICOM_PYdr_PYso_N_outgoing.empty() )
        {
            outfile__dynamic_array_ICOM_PYdr_PYso_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYdr_PYso_N_outgoing[0]), _dynamic_array_ICOM_PYdr_PYso_N_outgoing.size()*sizeof(_dynamic_array_ICOM_PYdr_PYso_N_outgoing[0]));
            outfile__dynamic_array_ICOM_PYdr_PYso_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYdr_PYso_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_post;
    outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_post.open(results_dir + "_dynamic_array_ICOM_PYso_PYdr__synaptic_post_3117492476", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_post.is_open())
    {
        if (! _dynamic_array_ICOM_PYso_PYdr__synaptic_post.empty() )
        {
            outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYso_PYdr__synaptic_post[0]), _dynamic_array_ICOM_PYso_PYdr__synaptic_post.size()*sizeof(_dynamic_array_ICOM_PYso_PYdr__synaptic_post[0]));
            outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYso_PYdr__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_pre;
    outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_pre.open(results_dir + "_dynamic_array_ICOM_PYso_PYdr__synaptic_pre_3006159977", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_pre.is_open())
    {
        if (! _dynamic_array_ICOM_PYso_PYdr__synaptic_pre.empty() )
        {
            outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYso_PYdr__synaptic_pre[0]), _dynamic_array_ICOM_PYso_PYdr__synaptic_pre.size()*sizeof(_dynamic_array_ICOM_PYso_PYdr__synaptic_pre[0]));
            outfile__dynamic_array_ICOM_PYso_PYdr__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYso_PYdr__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYso_PYdr_gCOM;
    outfile__dynamic_array_ICOM_PYso_PYdr_gCOM.open(results_dir + "_dynamic_array_ICOM_PYso_PYdr_gCOM_4268519655", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYso_PYdr_gCOM.is_open())
    {
        if (! _dynamic_array_ICOM_PYso_PYdr_gCOM.empty() )
        {
            outfile__dynamic_array_ICOM_PYso_PYdr_gCOM.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYso_PYdr_gCOM[0]), _dynamic_array_ICOM_PYso_PYdr_gCOM.size()*sizeof(_dynamic_array_ICOM_PYso_PYdr_gCOM[0]));
            outfile__dynamic_array_ICOM_PYso_PYdr_gCOM.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYso_PYdr_gCOM." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYso_PYdr_N_incoming;
    outfile__dynamic_array_ICOM_PYso_PYdr_N_incoming.open(results_dir + "_dynamic_array_ICOM_PYso_PYdr_N_incoming_2263571795", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYso_PYdr_N_incoming.is_open())
    {
        if (! _dynamic_array_ICOM_PYso_PYdr_N_incoming.empty() )
        {
            outfile__dynamic_array_ICOM_PYso_PYdr_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYso_PYdr_N_incoming[0]), _dynamic_array_ICOM_PYso_PYdr_N_incoming.size()*sizeof(_dynamic_array_ICOM_PYso_PYdr_N_incoming[0]));
            outfile__dynamic_array_ICOM_PYso_PYdr_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYso_PYdr_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_ICOM_PYso_PYdr_N_outgoing;
    outfile__dynamic_array_ICOM_PYso_PYdr_N_outgoing.open(results_dir + "_dynamic_array_ICOM_PYso_PYdr_N_outgoing_2717330825", ios::binary | ios::out);
    if(outfile__dynamic_array_ICOM_PYso_PYdr_N_outgoing.is_open())
    {
        if (! _dynamic_array_ICOM_PYso_PYdr_N_outgoing.empty() )
        {
            outfile__dynamic_array_ICOM_PYso_PYdr_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_ICOM_PYso_PYdr_N_outgoing[0]), _dynamic_array_ICOM_PYso_PYdr_N_outgoing.size()*sizeof(_dynamic_array_ICOM_PYso_PYdr_N_outgoing[0]));
            outfile__dynamic_array_ICOM_PYso_PYdr_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_ICOM_PYso_PYdr_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN__synaptic_post;
    outfile__dynamic_array_IGABAA_IN_IN__synaptic_post.open(results_dir + "_dynamic_array_IGABAA_IN_IN__synaptic_post_2377424692", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN__synaptic_post.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN__synaptic_post.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN__synaptic_post[0]), _dynamic_array_IGABAA_IN_IN__synaptic_post.size()*sizeof(_dynamic_array_IGABAA_IN_IN__synaptic_post[0]));
            outfile__dynamic_array_IGABAA_IN_IN__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN__synaptic_pre;
    outfile__dynamic_array_IGABAA_IN_IN__synaptic_pre.open(results_dir + "_dynamic_array_IGABAA_IN_IN__synaptic_pre_3091228672", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN__synaptic_pre[0]), _dynamic_array_IGABAA_IN_IN__synaptic_pre.size()*sizeof(_dynamic_array_IGABAA_IN_IN__synaptic_pre[0]));
            outfile__dynamic_array_IGABAA_IN_IN__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_alpha_GABAA;
    outfile__dynamic_array_IGABAA_IN_IN_alpha_GABAA.open(results_dir + "_dynamic_array_IGABAA_IN_IN_alpha_GABAA_3314396914", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_alpha_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_alpha_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_alpha_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_alpha_GABAA[0]), _dynamic_array_IGABAA_IN_IN_alpha_GABAA.size()*sizeof(_dynamic_array_IGABAA_IN_IN_alpha_GABAA[0]));
            outfile__dynamic_array_IGABAA_IN_IN_alpha_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_alpha_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_delay;
    outfile__dynamic_array_IGABAA_IN_IN_delay.open(results_dir + "_dynamic_array_IGABAA_IN_IN_delay_1250362046", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_delay.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_delay.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_delay.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_delay[0]), _dynamic_array_IGABAA_IN_IN_delay.size()*sizeof(_dynamic_array_IGABAA_IN_IN_delay[0]));
            outfile__dynamic_array_IGABAA_IN_IN_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_delay." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_delay_1;
    outfile__dynamic_array_IGABAA_IN_IN_delay_1.open(results_dir + "_dynamic_array_IGABAA_IN_IN_delay_1_3110927299", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_delay_1.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_delay_1.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_delay_1[0]), _dynamic_array_IGABAA_IN_IN_delay_1.size()*sizeof(_dynamic_array_IGABAA_IN_IN_delay_1[0]));
            outfile__dynamic_array_IGABAA_IN_IN_delay_1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_delay_1." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_deprFactor;
    outfile__dynamic_array_IGABAA_IN_IN_deprFactor.open(results_dir + "_dynamic_array_IGABAA_IN_IN_deprFactor_661231042", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_deprFactor.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_deprFactor.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_deprFactor.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_deprFactor[0]), _dynamic_array_IGABAA_IN_IN_deprFactor.size()*sizeof(_dynamic_array_IGABAA_IN_IN_deprFactor[0]));
            outfile__dynamic_array_IGABAA_IN_IN_deprFactor.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_deprFactor." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_EGABAA;
    outfile__dynamic_array_IGABAA_IN_IN_EGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_IN_EGABAA_4151778900", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_EGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_EGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_EGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_EGABAA[0]), _dynamic_array_IGABAA_IN_IN_EGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_IN_EGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_IN_EGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_EGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_gGABAA;
    outfile__dynamic_array_IGABAA_IN_IN_gGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_IN_gGABAA_3172186729", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_gGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_gGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_gGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_gGABAA[0]), _dynamic_array_IGABAA_IN_IN_gGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_IN_gGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_IN_gGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_gGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_N_incoming;
    outfile__dynamic_array_IGABAA_IN_IN_N_incoming.open(results_dir + "_dynamic_array_IGABAA_IN_IN_N_incoming_1865754739", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_N_incoming.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_N_incoming.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_N_incoming[0]), _dynamic_array_IGABAA_IN_IN_N_incoming.size()*sizeof(_dynamic_array_IGABAA_IN_IN_N_incoming[0]));
            outfile__dynamic_array_IGABAA_IN_IN_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_N_outgoing;
    outfile__dynamic_array_IGABAA_IN_IN_N_outgoing.open(results_dir + "_dynamic_array_IGABAA_IN_IN_N_outgoing_1210669225", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_N_outgoing.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_N_outgoing.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_N_outgoing[0]), _dynamic_array_IGABAA_IN_IN_N_outgoing.size()*sizeof(_dynamic_array_IGABAA_IN_IN_N_outgoing[0]));
            outfile__dynamic_array_IGABAA_IN_IN_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_Npost;
    outfile__dynamic_array_IGABAA_IN_IN_Npost.open(results_dir + "_dynamic_array_IGABAA_IN_IN_Npost_2721194113", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_Npost.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_Npost.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_Npost[0]), _dynamic_array_IGABAA_IN_IN_Npost.size()*sizeof(_dynamic_array_IGABAA_IN_IN_Npost[0]));
            outfile__dynamic_array_IGABAA_IN_IN_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_Npre;
    outfile__dynamic_array_IGABAA_IN_IN_Npre.open(results_dir + "_dynamic_array_IGABAA_IN_IN_Npre_1367562812", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_Npre.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_Npre.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_Npre[0]), _dynamic_array_IGABAA_IN_IN_Npre.size()*sizeof(_dynamic_array_IGABAA_IN_IN_Npre[0]));
            outfile__dynamic_array_IGABAA_IN_IN_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_propoCondMult;
    outfile__dynamic_array_IGABAA_IN_IN_propoCondMult.open(results_dir + "_dynamic_array_IGABAA_IN_IN_propoCondMult_2207902782", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_propoCondMult.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_propoCondMult.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_propoCondMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_propoCondMult[0]), _dynamic_array_IGABAA_IN_IN_propoCondMult.size()*sizeof(_dynamic_array_IGABAA_IN_IN_propoCondMult[0]));
            outfile__dynamic_array_IGABAA_IN_IN_propoCondMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_propoCondMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_propoTauMult;
    outfile__dynamic_array_IGABAA_IN_IN_propoTauMult.open(results_dir + "_dynamic_array_IGABAA_IN_IN_propoTauMult_3043920180", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_propoTauMult.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_propoTauMult.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_propoTauMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_propoTauMult[0]), _dynamic_array_IGABAA_IN_IN_propoTauMult.size()*sizeof(_dynamic_array_IGABAA_IN_IN_propoTauMult[0]));
            outfile__dynamic_array_IGABAA_IN_IN_propoTauMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_propoTauMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_radius;
    outfile__dynamic_array_IGABAA_IN_IN_radius.open(results_dir + "_dynamic_array_IGABAA_IN_IN_radius_1887479711", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_radius.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_radius.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_radius.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_radius[0]), _dynamic_array_IGABAA_IN_IN_radius.size()*sizeof(_dynamic_array_IGABAA_IN_IN_radius[0]));
            outfile__dynamic_array_IGABAA_IN_IN_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_radius." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_remove_recurrent_bool;
    outfile__dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.open(results_dir + "_dynamic_array_IGABAA_IN_IN_remove_recurrent_bool_2106675449", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_remove_recurrent_bool[0]), _dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.size()*sizeof(_dynamic_array_IGABAA_IN_IN_remove_recurrent_bool[0]));
            outfile__dynamic_array_IGABAA_IN_IN_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_res_GABAA;
    outfile__dynamic_array_IGABAA_IN_IN_res_GABAA.open(results_dir + "_dynamic_array_IGABAA_IN_IN_res_GABAA_3951315595", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_res_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_res_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_res_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_res_GABAA[0]), _dynamic_array_IGABAA_IN_IN_res_GABAA.size()*sizeof(_dynamic_array_IGABAA_IN_IN_res_GABAA[0]));
            outfile__dynamic_array_IGABAA_IN_IN_res_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_res_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_sGABAA;
    outfile__dynamic_array_IGABAA_IN_IN_sGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_IN_sGABAA_626326244", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_sGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_sGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_sGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_sGABAA[0]), _dynamic_array_IGABAA_IN_IN_sGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_IN_sGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_IN_sGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_sGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_spike_activity;
    outfile__dynamic_array_IGABAA_IN_IN_spike_activity.open(results_dir + "_dynamic_array_IGABAA_IN_IN_spike_activity_1317822633", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_spike_activity.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_spike_activity.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_spike_activity.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_spike_activity[0]), _dynamic_array_IGABAA_IN_IN_spike_activity.size()*sizeof(_dynamic_array_IGABAA_IN_IN_spike_activity[0]));
            outfile__dynamic_array_IGABAA_IN_IN_spike_activity.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_spike_activity." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_tauGABAA;
    outfile__dynamic_array_IGABAA_IN_IN_tauGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_IN_tauGABAA_1040147502", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_tauGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_tauGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_tauGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_tauGABAA[0]), _dynamic_array_IGABAA_IN_IN_tauGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_IN_tauGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_IN_tauGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_tauGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_IN_tauRes_GABAA;
    outfile__dynamic_array_IGABAA_IN_IN_tauRes_GABAA.open(results_dir + "_dynamic_array_IGABAA_IN_IN_tauRes_GABAA_2128080645", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_IN_tauRes_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_IN_tauRes_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_IN_tauRes_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_IN_tauRes_GABAA[0]), _dynamic_array_IGABAA_IN_IN_tauRes_GABAA.size()*sizeof(_dynamic_array_IGABAA_IN_IN_tauRes_GABAA[0]));
            outfile__dynamic_array_IGABAA_IN_IN_tauRes_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_IN_tauRes_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso__synaptic_post;
    outfile__dynamic_array_IGABAA_IN_PYso__synaptic_post.open(results_dir + "_dynamic_array_IGABAA_IN_PYso__synaptic_post_649312757", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso__synaptic_post.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso__synaptic_post.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso__synaptic_post[0]), _dynamic_array_IGABAA_IN_PYso__synaptic_post.size()*sizeof(_dynamic_array_IGABAA_IN_PYso__synaptic_post[0]));
            outfile__dynamic_array_IGABAA_IN_PYso__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso__synaptic_pre;
    outfile__dynamic_array_IGABAA_IN_PYso__synaptic_pre.open(results_dir + "_dynamic_array_IGABAA_IN_PYso__synaptic_pre_1851519023", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso__synaptic_pre[0]), _dynamic_array_IGABAA_IN_PYso__synaptic_pre.size()*sizeof(_dynamic_array_IGABAA_IN_PYso__synaptic_pre[0]));
            outfile__dynamic_array_IGABAA_IN_PYso__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_alpha_GABAA;
    outfile__dynamic_array_IGABAA_IN_PYso_alpha_GABAA.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_alpha_GABAA_692886791", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_alpha_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_alpha_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_alpha_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_alpha_GABAA[0]), _dynamic_array_IGABAA_IN_PYso_alpha_GABAA.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_alpha_GABAA[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_alpha_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_alpha_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_delay;
    outfile__dynamic_array_IGABAA_IN_PYso_delay.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_delay_1064572579", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_delay.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_delay.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_delay.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_delay[0]), _dynamic_array_IGABAA_IN_PYso_delay.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_delay[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_delay." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_delay_1;
    outfile__dynamic_array_IGABAA_IN_PYso_delay_1.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_delay_1_1557804625", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_delay_1.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_delay_1.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_delay_1[0]), _dynamic_array_IGABAA_IN_PYso_delay_1.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_delay_1[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_delay_1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_delay_1." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_deprFactor;
    outfile__dynamic_array_IGABAA_IN_PYso_deprFactor.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_deprFactor_2242397699", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_deprFactor.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_deprFactor.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_deprFactor.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_deprFactor[0]), _dynamic_array_IGABAA_IN_PYso_deprFactor.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_deprFactor[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_deprFactor.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_deprFactor." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_EGABAA;
    outfile__dynamic_array_IGABAA_IN_PYso_EGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_EGABAA_2483323289", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_EGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_EGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_EGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_EGABAA[0]), _dynamic_array_IGABAA_IN_PYso_EGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_EGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_EGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_EGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_gGABAA;
    outfile__dynamic_array_IGABAA_IN_PYso_gGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_gGABAA_3730843044", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_gGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_gGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_gGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_gGABAA[0]), _dynamic_array_IGABAA_IN_PYso_gGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_gGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_gGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_gGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_N_incoming;
    outfile__dynamic_array_IGABAA_IN_PYso_N_incoming.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_N_incoming_3455384498", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_N_incoming.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_N_incoming.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_N_incoming[0]), _dynamic_array_IGABAA_IN_PYso_N_incoming.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_N_incoming[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_N_outgoing;
    outfile__dynamic_array_IGABAA_IN_PYso_N_outgoing.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_N_outgoing_3941108584", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_N_outgoing.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_N_outgoing.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_N_outgoing[0]), _dynamic_array_IGABAA_IN_PYso_N_outgoing.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_N_outgoing[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_Npost;
    outfile__dynamic_array_IGABAA_IN_PYso_Npost.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_Npost_3619763356", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_Npost.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_Npost.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_Npost[0]), _dynamic_array_IGABAA_IN_PYso_Npost.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_Npost[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_Npre;
    outfile__dynamic_array_IGABAA_IN_PYso_Npre.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_Npre_3358907390", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_Npre.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_Npre.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_Npre[0]), _dynamic_array_IGABAA_IN_PYso_Npre.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_Npre[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_propoCondMult;
    outfile__dynamic_array_IGABAA_IN_PYso_propoCondMult.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_propoCondMult_1434612753", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_propoCondMult.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_propoCondMult.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_propoCondMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_propoCondMult[0]), _dynamic_array_IGABAA_IN_PYso_propoCondMult.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_propoCondMult[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_propoCondMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_propoCondMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_propoTauMult;
    outfile__dynamic_array_IGABAA_IN_PYso_propoTauMult.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_propoTauMult_2018856596", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_propoTauMult.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_propoTauMult.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_propoTauMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_propoTauMult[0]), _dynamic_array_IGABAA_IN_PYso_propoTauMult.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_propoTauMult[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_propoTauMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_propoTauMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_radius;
    outfile__dynamic_array_IGABAA_IN_PYso_radius.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_radius_334706770", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_radius.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_radius.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_radius.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_radius[0]), _dynamic_array_IGABAA_IN_PYso_radius.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_radius[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_radius." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool;
    outfile__dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool_3997458681", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool[0]), _dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_res_GABAA;
    outfile__dynamic_array_IGABAA_IN_PYso_res_GABAA.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_res_GABAA_1294399912", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_res_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_res_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_res_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_res_GABAA[0]), _dynamic_array_IGABAA_IN_PYso_res_GABAA.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_res_GABAA[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_res_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_res_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_sGABAA;
    outfile__dynamic_array_IGABAA_IN_PYso_sGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_sGABAA_1176986921", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_sGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_sGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_sGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_sGABAA[0]), _dynamic_array_IGABAA_IN_PYso_sGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_sGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_sGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_sGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_spike_activity;
    outfile__dynamic_array_IGABAA_IN_PYso_spike_activity.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_spike_activity_3851109992", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_spike_activity.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_spike_activity.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_spike_activity.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_spike_activity[0]), _dynamic_array_IGABAA_IN_PYso_spike_activity.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_spike_activity[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_spike_activity.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_spike_activity." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_tauGABAA;
    outfile__dynamic_array_IGABAA_IN_PYso_tauGABAA.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_tauGABAA_588980363", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_tauGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_tauGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_tauGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_tauGABAA[0]), _dynamic_array_IGABAA_IN_PYso_tauGABAA.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_tauGABAA[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_tauGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_tauGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_IN_PYso_tauRes_GABAA;
    outfile__dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.open(results_dir + "_dynamic_array_IGABAA_IN_PYso_tauRes_GABAA_3018598565", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_IN_PYso_tauRes_GABAA[0]), _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.size()*sizeof(_dynamic_array_IGABAA_IN_PYso_tauRes_GABAA[0]));
            outfile__dynamic_array_IGABAA_IN_PYso_tauRes_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_IN_PYso_tauRes_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE__synaptic_post;
    outfile__dynamic_array_IGABAA_RE_RE__synaptic_post.open(results_dir + "_dynamic_array_IGABAA_RE_RE__synaptic_post_234328735", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE__synaptic_post.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE__synaptic_post.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE__synaptic_post[0]), _dynamic_array_IGABAA_RE_RE__synaptic_post.size()*sizeof(_dynamic_array_IGABAA_RE_RE__synaptic_post[0]));
            outfile__dynamic_array_IGABAA_RE_RE__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE__synaptic_pre;
    outfile__dynamic_array_IGABAA_RE_RE__synaptic_pre.open(results_dir + "_dynamic_array_IGABAA_RE_RE__synaptic_pre_2653554837", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE__synaptic_pre[0]), _dynamic_array_IGABAA_RE_RE__synaptic_pre.size()*sizeof(_dynamic_array_IGABAA_RE_RE__synaptic_pre[0]));
            outfile__dynamic_array_IGABAA_RE_RE__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_EGABAA;
    outfile__dynamic_array_IGABAA_RE_RE_EGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_RE_EGABAA_3659730968", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_EGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_EGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_EGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_EGABAA[0]), _dynamic_array_IGABAA_RE_RE_EGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_RE_EGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_RE_EGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_EGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_gGABAA;
    outfile__dynamic_array_IGABAA_RE_RE_gGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_RE_gGABAA_2420618277", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_gGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_gGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_gGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_gGABAA[0]), _dynamic_array_IGABAA_RE_RE_gGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_RE_gGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_RE_gGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_gGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_N_incoming;
    outfile__dynamic_array_IGABAA_RE_RE_N_incoming.open(results_dir + "_dynamic_array_IGABAA_RE_RE_N_incoming_794896089", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_N_incoming.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_N_incoming.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_N_incoming[0]), _dynamic_array_IGABAA_RE_RE_N_incoming.size()*sizeof(_dynamic_array_IGABAA_RE_RE_N_incoming[0]));
            outfile__dynamic_array_IGABAA_RE_RE_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_N_outgoing;
    outfile__dynamic_array_IGABAA_RE_RE_N_outgoing.open(results_dir + "_dynamic_array_IGABAA_RE_RE_N_outgoing_142431747", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_N_outgoing.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_N_outgoing.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_N_outgoing[0]), _dynamic_array_IGABAA_RE_RE_N_outgoing.size()*sizeof(_dynamic_array_IGABAA_RE_RE_N_outgoing[0]));
            outfile__dynamic_array_IGABAA_RE_RE_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_Npost;
    outfile__dynamic_array_IGABAA_RE_RE_Npost.open(results_dir + "_dynamic_array_IGABAA_RE_RE_Npost_4107526526", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_Npost.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_Npost.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_Npost[0]), _dynamic_array_IGABAA_RE_RE_Npost.size()*sizeof(_dynamic_array_IGABAA_RE_RE_Npost[0]));
            outfile__dynamic_array_IGABAA_RE_RE_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_Npre;
    outfile__dynamic_array_IGABAA_RE_RE_Npre.open(results_dir + "_dynamic_array_IGABAA_RE_RE_Npre_59147785", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_Npre.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_Npre.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_Npre[0]), _dynamic_array_IGABAA_RE_RE_Npre.size()*sizeof(_dynamic_array_IGABAA_RE_RE_Npre[0]));
            outfile__dynamic_array_IGABAA_RE_RE_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_P_GABAA;
    outfile__dynamic_array_IGABAA_RE_RE_P_GABAA.open(results_dir + "_dynamic_array_IGABAA_RE_RE_P_GABAA_2243949125", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_P_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_P_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_P_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_P_GABAA[0]), _dynamic_array_IGABAA_RE_RE_P_GABAA.size()*sizeof(_dynamic_array_IGABAA_RE_RE_P_GABAA[0]));
            outfile__dynamic_array_IGABAA_RE_RE_P_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_P_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_propoCondMult;
    outfile__dynamic_array_IGABAA_RE_RE_propoCondMult.open(results_dir + "_dynamic_array_IGABAA_RE_RE_propoCondMult_2784201899", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_propoCondMult.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_propoCondMult.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_propoCondMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_propoCondMult[0]), _dynamic_array_IGABAA_RE_RE_propoCondMult.size()*sizeof(_dynamic_array_IGABAA_RE_RE_propoCondMult[0]));
            outfile__dynamic_array_IGABAA_RE_RE_propoCondMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_propoCondMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_propoTauMult;
    outfile__dynamic_array_IGABAA_RE_RE_propoTauMult.open(results_dir + "_dynamic_array_IGABAA_RE_RE_propoTauMult_104771588", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_propoTauMult.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_propoTauMult.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_propoTauMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_propoTauMult[0]), _dynamic_array_IGABAA_RE_RE_propoTauMult.size()*sizeof(_dynamic_array_IGABAA_RE_RE_propoTauMult[0]));
            outfile__dynamic_array_IGABAA_RE_RE_propoTauMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_propoTauMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_sGABAA;
    outfile__dynamic_array_IGABAA_RE_RE_sGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_RE_sGABAA_134280360", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_sGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_sGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_sGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_sGABAA[0]), _dynamic_array_IGABAA_RE_RE_sGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_RE_sGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_RE_sGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_sGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_RE_tauGABAA;
    outfile__dynamic_array_IGABAA_RE_RE_tauGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_RE_tauGABAA_2257282164", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_RE_tauGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_RE_tauGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_RE_tauGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_RE_tauGABAA[0]), _dynamic_array_IGABAA_RE_RE_tauGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_RE_tauGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_RE_tauGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_RE_tauGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC__synaptic_post;
    outfile__dynamic_array_IGABAA_RE_TC__synaptic_post.open(results_dir + "_dynamic_array_IGABAA_RE_TC__synaptic_post_4182920307", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC__synaptic_post.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC__synaptic_post.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC__synaptic_post[0]), _dynamic_array_IGABAA_RE_TC__synaptic_post.size()*sizeof(_dynamic_array_IGABAA_RE_TC__synaptic_post[0]));
            outfile__dynamic_array_IGABAA_RE_TC__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC__synaptic_pre;
    outfile__dynamic_array_IGABAA_RE_TC__synaptic_pre.open(results_dir + "_dynamic_array_IGABAA_RE_TC__synaptic_pre_4014060931", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC__synaptic_pre[0]), _dynamic_array_IGABAA_RE_TC__synaptic_pre.size()*sizeof(_dynamic_array_IGABAA_RE_TC__synaptic_pre[0]));
            outfile__dynamic_array_IGABAA_RE_TC__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_EGABAA;
    outfile__dynamic_array_IGABAA_RE_TC_EGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_TC_EGABAA_1867931669", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_EGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_EGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_EGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_EGABAA[0]), _dynamic_array_IGABAA_RE_TC_EGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_TC_EGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_TC_EGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_EGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_gGABAA;
    outfile__dynamic_array_IGABAA_RE_TC_gGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_TC_gGABAA_624084008", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_gGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_gGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_gGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_gGABAA[0]), _dynamic_array_IGABAA_RE_TC_gGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_TC_gGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_TC_gGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_gGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_N_incoming;
    outfile__dynamic_array_IGABAA_RE_TC_N_incoming.open(results_dir + "_dynamic_array_IGABAA_RE_TC_N_incoming_1622760676", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_N_incoming.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_N_incoming.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_N_incoming[0]), _dynamic_array_IGABAA_RE_TC_N_incoming.size()*sizeof(_dynamic_array_IGABAA_RE_TC_N_incoming[0]));
            outfile__dynamic_array_IGABAA_RE_TC_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_N_outgoing;
    outfile__dynamic_array_IGABAA_RE_TC_N_outgoing.open(results_dir + "_dynamic_array_IGABAA_RE_TC_N_outgoing_1202007102", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_N_outgoing.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_N_outgoing.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_N_outgoing[0]), _dynamic_array_IGABAA_RE_TC_N_outgoing.size()*sizeof(_dynamic_array_IGABAA_RE_TC_N_outgoing[0]));
            outfile__dynamic_array_IGABAA_RE_TC_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_Npost;
    outfile__dynamic_array_IGABAA_RE_TC_Npost.open(results_dir + "_dynamic_array_IGABAA_RE_TC_Npost_1366089155", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_Npost.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_Npost.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_Npost[0]), _dynamic_array_IGABAA_RE_TC_Npost.size()*sizeof(_dynamic_array_IGABAA_RE_TC_Npost[0]));
            outfile__dynamic_array_IGABAA_RE_TC_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_Npre;
    outfile__dynamic_array_IGABAA_RE_TC_Npre.open(results_dir + "_dynamic_array_IGABAA_RE_TC_Npre_3054456878", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_Npre.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_Npre.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_Npre[0]), _dynamic_array_IGABAA_RE_TC_Npre.size()*sizeof(_dynamic_array_IGABAA_RE_TC_Npre[0]));
            outfile__dynamic_array_IGABAA_RE_TC_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_P_GABAA;
    outfile__dynamic_array_IGABAA_RE_TC_P_GABAA.open(results_dir + "_dynamic_array_IGABAA_RE_TC_P_GABAA_4223400840", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_P_GABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_P_GABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_P_GABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_P_GABAA[0]), _dynamic_array_IGABAA_RE_TC_P_GABAA.size()*sizeof(_dynamic_array_IGABAA_RE_TC_P_GABAA[0]));
            outfile__dynamic_array_IGABAA_RE_TC_P_GABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_P_GABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_propoCondMult;
    outfile__dynamic_array_IGABAA_RE_TC_propoCondMult.open(results_dir + "_dynamic_array_IGABAA_RE_TC_propoCondMult_3566743997", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_propoCondMult.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_propoCondMult.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_propoCondMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_propoCondMult[0]), _dynamic_array_IGABAA_RE_TC_propoCondMult.size()*sizeof(_dynamic_array_IGABAA_RE_TC_propoCondMult[0]));
            outfile__dynamic_array_IGABAA_RE_TC_propoCondMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_propoCondMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_propoTauMult;
    outfile__dynamic_array_IGABAA_RE_TC_propoTauMult.open(results_dir + "_dynamic_array_IGABAA_RE_TC_propoTauMult_3691392832", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_propoTauMult.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_propoTauMult.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_propoTauMult.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_propoTauMult[0]), _dynamic_array_IGABAA_RE_TC_propoTauMult.size()*sizeof(_dynamic_array_IGABAA_RE_TC_propoTauMult[0]));
            outfile__dynamic_array_IGABAA_RE_TC_propoTauMult.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_propoTauMult." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_sGABAA;
    outfile__dynamic_array_IGABAA_RE_TC_sGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_TC_sGABAA_3178595493", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_sGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_sGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_sGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_sGABAA[0]), _dynamic_array_IGABAA_RE_TC_sGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_TC_sGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_TC_sGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_sGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAA_RE_TC_tauGABAA;
    outfile__dynamic_array_IGABAA_RE_TC_tauGABAA.open(results_dir + "_dynamic_array_IGABAA_RE_TC_tauGABAA_1663098480", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAA_RE_TC_tauGABAA.is_open())
    {
        if (! _dynamic_array_IGABAA_RE_TC_tauGABAA.empty() )
        {
            outfile__dynamic_array_IGABAA_RE_TC_tauGABAA.write(reinterpret_cast<char*>(&_dynamic_array_IGABAA_RE_TC_tauGABAA[0]), _dynamic_array_IGABAA_RE_TC_tauGABAA.size()*sizeof(_dynamic_array_IGABAA_RE_TC_tauGABAA[0]));
            outfile__dynamic_array_IGABAA_RE_TC_tauGABAA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAA_RE_TC_tauGABAA." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC__synaptic_post;
    outfile__dynamic_array_IGABAB_RE_TC__synaptic_post.open(results_dir + "_dynamic_array_IGABAB_RE_TC__synaptic_post_3224035635", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC__synaptic_post.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC__synaptic_post.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC__synaptic_post[0]), _dynamic_array_IGABAB_RE_TC__synaptic_post.size()*sizeof(_dynamic_array_IGABAB_RE_TC__synaptic_post[0]));
            outfile__dynamic_array_IGABAB_RE_TC__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC__synaptic_pre;
    outfile__dynamic_array_IGABAB_RE_TC__synaptic_pre.open(results_dir + "_dynamic_array_IGABAB_RE_TC__synaptic_pre_2489268064", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC__synaptic_pre[0]), _dynamic_array_IGABAB_RE_TC__synaptic_pre.size()*sizeof(_dynamic_array_IGABAB_RE_TC__synaptic_pre[0]));
            outfile__dynamic_array_IGABAB_RE_TC__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_EGABAB;
    outfile__dynamic_array_IGABAB_RE_TC_EGABAB.open(results_dir + "_dynamic_array_IGABAB_RE_TC_EGABAB_2319324276", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_EGABAB.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_EGABAB.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_EGABAB.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_EGABAB[0]), _dynamic_array_IGABAB_RE_TC_EGABAB.size()*sizeof(_dynamic_array_IGABAB_RE_TC_EGABAB[0]));
            outfile__dynamic_array_IGABAB_RE_TC_EGABAB.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_EGABAB." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_gGABAB;
    outfile__dynamic_array_IGABAB_RE_TC_gGABAB.open(results_dir + "_dynamic_array_IGABAB_RE_TC_gGABAB_3227168841", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_gGABAB.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_gGABAB.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_gGABAB.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_gGABAB[0]), _dynamic_array_IGABAB_RE_TC_gGABAB.size()*sizeof(_dynamic_array_IGABAB_RE_TC_gGABAB[0]));
            outfile__dynamic_array_IGABAB_RE_TC_gGABAB.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_gGABAB." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_K1;
    outfile__dynamic_array_IGABAB_RE_TC_K1.open(results_dir + "_dynamic_array_IGABAB_RE_TC_K1_2309617870", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_K1.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_K1.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_K1.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_K1[0]), _dynamic_array_IGABAB_RE_TC_K1.size()*sizeof(_dynamic_array_IGABAB_RE_TC_K1[0]));
            outfile__dynamic_array_IGABAB_RE_TC_K1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_K1." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_K2;
    outfile__dynamic_array_IGABAB_RE_TC_K2.open(results_dir + "_dynamic_array_IGABAB_RE_TC_K2_278964596", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_K2.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_K2.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_K2.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_K2[0]), _dynamic_array_IGABAB_RE_TC_K2.size()*sizeof(_dynamic_array_IGABAB_RE_TC_K2[0]));
            outfile__dynamic_array_IGABAB_RE_TC_K2.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_K2." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_K3;
    outfile__dynamic_array_IGABAB_RE_TC_K3.open(results_dir + "_dynamic_array_IGABAB_RE_TC_K3_1739037154", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_K3.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_K3.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_K3.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_K3[0]), _dynamic_array_IGABAB_RE_TC_K3.size()*sizeof(_dynamic_array_IGABAB_RE_TC_K3[0]));
            outfile__dynamic_array_IGABAB_RE_TC_K3.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_K3." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_K4;
    outfile__dynamic_array_IGABAB_RE_TC_K4.open(results_dir + "_dynamic_array_IGABAB_RE_TC_K4_4190309441", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_K4.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_K4.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_K4.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_K4[0]), _dynamic_array_IGABAB_RE_TC_K4.size()*sizeof(_dynamic_array_IGABAB_RE_TC_K4[0]));
            outfile__dynamic_array_IGABAB_RE_TC_K4.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_K4." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_N_incoming;
    outfile__dynamic_array_IGABAB_RE_TC_N_incoming.open(results_dir + "_dynamic_array_IGABAB_RE_TC_N_incoming_2319417734", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_N_incoming.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_N_incoming.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_N_incoming[0]), _dynamic_array_IGABAB_RE_TC_N_incoming.size()*sizeof(_dynamic_array_IGABAB_RE_TC_N_incoming[0]));
            outfile__dynamic_array_IGABAB_RE_TC_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_N_outgoing;
    outfile__dynamic_array_IGABAB_RE_TC_N_outgoing.open(results_dir + "_dynamic_array_IGABAB_RE_TC_N_outgoing_2904813916", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_N_outgoing.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_N_outgoing.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_N_outgoing[0]), _dynamic_array_IGABAB_RE_TC_N_outgoing.size()*sizeof(_dynamic_array_IGABAB_RE_TC_N_outgoing[0]));
            outfile__dynamic_array_IGABAB_RE_TC_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_Npost;
    outfile__dynamic_array_IGABAB_RE_TC_Npost.open(results_dir + "_dynamic_array_IGABAB_RE_TC_Npost_3970347277", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_Npost.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_Npost.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_Npost.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_Npost[0]), _dynamic_array_IGABAB_RE_TC_Npost.size()*sizeof(_dynamic_array_IGABAB_RE_TC_Npost[0]));
            outfile__dynamic_array_IGABAB_RE_TC_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_Npost." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_Npre;
    outfile__dynamic_array_IGABAB_RE_TC_Npre.open(results_dir + "_dynamic_array_IGABAB_RE_TC_Npre_3247538910", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_Npre.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_Npre.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_Npre.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_Npre[0]), _dynamic_array_IGABAB_RE_TC_Npre.size()*sizeof(_dynamic_array_IGABAB_RE_TC_Npre[0]));
            outfile__dynamic_array_IGABAB_RE_TC_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_Npre." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_rGABAB;
    outfile__dynamic_array_IGABAB_RE_TC_rGABAB.open(results_dir + "_dynamic_array_IGABAB_RE_TC_rGABAB_2470519649", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_rGABAB.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_rGABAB.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_rGABAB.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_rGABAB[0]), _dynamic_array_IGABAB_RE_TC_rGABAB.size()*sizeof(_dynamic_array_IGABAB_RE_TC_rGABAB[0]));
            outfile__dynamic_array_IGABAB_RE_TC_rGABAB.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_rGABAB." << endl;
    }
    ofstream outfile__dynamic_array_IGABAB_RE_TC_sGABAB;
    outfile__dynamic_array_IGABAB_RE_TC_sGABAB.open(results_dir + "_dynamic_array_IGABAB_RE_TC_sGABAB_1478357188", ios::binary | ios::out);
    if(outfile__dynamic_array_IGABAB_RE_TC_sGABAB.is_open())
    {
        if (! _dynamic_array_IGABAB_RE_TC_sGABAB.empty() )
        {
            outfile__dynamic_array_IGABAB_RE_TC_sGABAB.write(reinterpret_cast<char*>(&_dynamic_array_IGABAB_RE_TC_sGABAB[0]), _dynamic_array_IGABAB_RE_TC_sGABAB.size()*sizeof(_dynamic_array_IGABAB_RE_TC_sGABAB[0]));
            outfile__dynamic_array_IGABAB_RE_TC_sGABAB.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IGABAB_RE_TC_sGABAB." << endl;
    }
    ofstream outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_post;
    outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_post.open(results_dir + "_dynamic_array_IKNa_PYso_PYdr__synaptic_post_714788998", ios::binary | ios::out);
    if(outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_post.is_open())
    {
        if (! _dynamic_array_IKNa_PYso_PYdr__synaptic_post.empty() )
        {
            outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_IKNa_PYso_PYdr__synaptic_post[0]), _dynamic_array_IKNa_PYso_PYdr__synaptic_post.size()*sizeof(_dynamic_array_IKNa_PYso_PYdr__synaptic_post[0]));
            outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IKNa_PYso_PYdr__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_pre;
    outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_pre.open(results_dir + "_dynamic_array_IKNa_PYso_PYdr__synaptic_pre_4046044132", ios::binary | ios::out);
    if(outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_pre.is_open())
    {
        if (! _dynamic_array_IKNa_PYso_PYdr__synaptic_pre.empty() )
        {
            outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_IKNa_PYso_PYdr__synaptic_pre[0]), _dynamic_array_IKNa_PYso_PYdr__synaptic_pre.size()*sizeof(_dynamic_array_IKNa_PYso_PYdr__synaptic_pre[0]));
            outfile__dynamic_array_IKNa_PYso_PYdr__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IKNa_PYso_PYdr__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_IKNa_PYso_PYdr_concNa;
    outfile__dynamic_array_IKNa_PYso_PYdr_concNa.open(results_dir + "_dynamic_array_IKNa_PYso_PYdr_concNa_3939991475", ios::binary | ios::out);
    if(outfile__dynamic_array_IKNa_PYso_PYdr_concNa.is_open())
    {
        if (! _dynamic_array_IKNa_PYso_PYdr_concNa.empty() )
        {
            outfile__dynamic_array_IKNa_PYso_PYdr_concNa.write(reinterpret_cast<char*>(&_dynamic_array_IKNa_PYso_PYdr_concNa[0]), _dynamic_array_IKNa_PYso_PYdr_concNa.size()*sizeof(_dynamic_array_IKNa_PYso_PYdr_concNa[0]));
            outfile__dynamic_array_IKNa_PYso_PYdr_concNa.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IKNa_PYso_PYdr_concNa." << endl;
    }
    ofstream outfile__dynamic_array_IKNa_PYso_PYdr_hNa_local;
    outfile__dynamic_array_IKNa_PYso_PYdr_hNa_local.open(results_dir + "_dynamic_array_IKNa_PYso_PYdr_hNa_local_2476603758", ios::binary | ios::out);
    if(outfile__dynamic_array_IKNa_PYso_PYdr_hNa_local.is_open())
    {
        if (! _dynamic_array_IKNa_PYso_PYdr_hNa_local.empty() )
        {
            outfile__dynamic_array_IKNa_PYso_PYdr_hNa_local.write(reinterpret_cast<char*>(&_dynamic_array_IKNa_PYso_PYdr_hNa_local[0]), _dynamic_array_IKNa_PYso_PYdr_hNa_local.size()*sizeof(_dynamic_array_IKNa_PYso_PYdr_hNa_local[0]));
            outfile__dynamic_array_IKNa_PYso_PYdr_hNa_local.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IKNa_PYso_PYdr_hNa_local." << endl;
    }
    ofstream outfile__dynamic_array_IKNa_PYso_PYdr_N_incoming;
    outfile__dynamic_array_IKNa_PYso_PYdr_N_incoming.open(results_dir + "_dynamic_array_IKNa_PYso_PYdr_N_incoming_2733157626", ios::binary | ios::out);
    if(outfile__dynamic_array_IKNa_PYso_PYdr_N_incoming.is_open())
    {
        if (! _dynamic_array_IKNa_PYso_PYdr_N_incoming.empty() )
        {
            outfile__dynamic_array_IKNa_PYso_PYdr_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_IKNa_PYso_PYdr_N_incoming[0]), _dynamic_array_IKNa_PYso_PYdr_N_incoming.size()*sizeof(_dynamic_array_IKNa_PYso_PYdr_N_incoming[0]));
            outfile__dynamic_array_IKNa_PYso_PYdr_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IKNa_PYso_PYdr_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_IKNa_PYso_PYdr_N_outgoing;
    outfile__dynamic_array_IKNa_PYso_PYdr_N_outgoing.open(results_dir + "_dynamic_array_IKNa_PYso_PYdr_N_outgoing_2247416864", ios::binary | ios::out);
    if(outfile__dynamic_array_IKNa_PYso_PYdr_N_outgoing.is_open())
    {
        if (! _dynamic_array_IKNa_PYso_PYdr_N_outgoing.empty() )
        {
            outfile__dynamic_array_IKNa_PYso_PYdr_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_IKNa_PYso_PYdr_N_outgoing[0]), _dynamic_array_IKNa_PYso_PYdr_N_outgoing.size()*sizeof(_dynamic_array_IKNa_PYso_PYdr_N_outgoing[0]));
            outfile__dynamic_array_IKNa_PYso_PYdr_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IKNa_PYso_PYdr_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_IN_spikemon_i;
    outfile__dynamic_array_IN_spikemon_i.open(results_dir + "_dynamic_array_IN_spikemon_i_3810621142", ios::binary | ios::out);
    if(outfile__dynamic_array_IN_spikemon_i.is_open())
    {
        if (! _dynamic_array_IN_spikemon_i.empty() )
        {
            outfile__dynamic_array_IN_spikemon_i.write(reinterpret_cast<char*>(&_dynamic_array_IN_spikemon_i[0]), _dynamic_array_IN_spikemon_i.size()*sizeof(_dynamic_array_IN_spikemon_i[0]));
            outfile__dynamic_array_IN_spikemon_i.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IN_spikemon_i." << endl;
    }
    ofstream outfile__dynamic_array_IN_spikemon_t;
    outfile__dynamic_array_IN_spikemon_t.open(results_dir + "_dynamic_array_IN_spikemon_t_2150046223", ios::binary | ios::out);
    if(outfile__dynamic_array_IN_spikemon_t.is_open())
    {
        if (! _dynamic_array_IN_spikemon_t.empty() )
        {
            outfile__dynamic_array_IN_spikemon_t.write(reinterpret_cast<char*>(&_dynamic_array_IN_spikemon_t[0]), _dynamic_array_IN_spikemon_t.size()*sizeof(_dynamic_array_IN_spikemon_t[0]));
            outfile__dynamic_array_IN_spikemon_t.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_IN_spikemon_t." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN__synaptic_post;
    outfile__dynamic_array_INMDA_PYso_IN__synaptic_post.open(results_dir + "_dynamic_array_INMDA_PYso_IN__synaptic_post_3560797540", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN__synaptic_post.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN__synaptic_post.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN__synaptic_post[0]), _dynamic_array_INMDA_PYso_IN__synaptic_post.size()*sizeof(_dynamic_array_INMDA_PYso_IN__synaptic_post[0]));
            outfile__dynamic_array_INMDA_PYso_IN__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN__synaptic_pre;
    outfile__dynamic_array_INMDA_PYso_IN__synaptic_pre.open(results_dir + "_dynamic_array_INMDA_PYso_IN__synaptic_pre_2190946172", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN__synaptic_pre.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN__synaptic_pre.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN__synaptic_pre[0]), _dynamic_array_INMDA_PYso_IN__synaptic_pre.size()*sizeof(_dynamic_array_INMDA_PYso_IN__synaptic_pre[0]));
            outfile__dynamic_array_INMDA_PYso_IN__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_alphaS_NMDA;
    outfile__dynamic_array_INMDA_PYso_IN_alphaS_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_alphaS_NMDA_121155512", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_alphaS_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_alphaS_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_alphaS_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_alphaS_NMDA[0]), _dynamic_array_INMDA_PYso_IN_alphaS_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_alphaS_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_alphaS_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_alphaS_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_alphaX_NMDA;
    outfile__dynamic_array_INMDA_PYso_IN_alphaX_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_alphaX_NMDA_1845452667", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_alphaX_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_alphaX_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_alphaX_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_alphaX_NMDA[0]), _dynamic_array_INMDA_PYso_IN_alphaX_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_alphaX_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_alphaX_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_alphaX_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_delay;
    outfile__dynamic_array_INMDA_PYso_IN_delay.open(results_dir + "_dynamic_array_INMDA_PYso_IN_delay_1412353028", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_delay.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_delay.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_delay.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_delay[0]), _dynamic_array_INMDA_PYso_IN_delay.size()*sizeof(_dynamic_array_INMDA_PYso_IN_delay[0]));
            outfile__dynamic_array_INMDA_PYso_IN_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_delay." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_delay_1;
    outfile__dynamic_array_INMDA_PYso_IN_delay_1.open(results_dir + "_dynamic_array_INMDA_PYso_IN_delay_1_3684322543", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_delay_1.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_delay_1.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_delay_1[0]), _dynamic_array_INMDA_PYso_IN_delay_1.size()*sizeof(_dynamic_array_INMDA_PYso_IN_delay_1[0]));
            outfile__dynamic_array_INMDA_PYso_IN_delay_1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_delay_1." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_deprFactor;
    outfile__dynamic_array_INMDA_PYso_IN_deprFactor.open(results_dir + "_dynamic_array_INMDA_PYso_IN_deprFactor_3915662706", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_deprFactor.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_deprFactor.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_deprFactor.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_deprFactor[0]), _dynamic_array_INMDA_PYso_IN_deprFactor.size()*sizeof(_dynamic_array_INMDA_PYso_IN_deprFactor[0]));
            outfile__dynamic_array_INMDA_PYso_IN_deprFactor.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_deprFactor." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_ENMDA;
    outfile__dynamic_array_INMDA_PYso_IN_ENMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_ENMDA_650290721", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_ENMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_ENMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_ENMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_ENMDA[0]), _dynamic_array_INMDA_PYso_IN_ENMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_ENMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_ENMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_ENMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_gNMDA;
    outfile__dynamic_array_INMDA_PYso_IN_gNMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_gNMDA_2646858309", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_gNMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_gNMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_gNMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_gNMDA[0]), _dynamic_array_INMDA_PYso_IN_gNMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_gNMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_gNMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_gNMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_N_incoming;
    outfile__dynamic_array_INMDA_PYso_IN_N_incoming.open(results_dir + "_dynamic_array_INMDA_PYso_IN_N_incoming_2704865475", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_N_incoming.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_N_incoming.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_N_incoming[0]), _dynamic_array_INMDA_PYso_IN_N_incoming.size()*sizeof(_dynamic_array_INMDA_PYso_IN_N_incoming[0]));
            outfile__dynamic_array_INMDA_PYso_IN_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_N_outgoing;
    outfile__dynamic_array_INMDA_PYso_IN_N_outgoing.open(results_dir + "_dynamic_array_INMDA_PYso_IN_N_outgoing_2250541081", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_N_outgoing.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_N_outgoing.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_N_outgoing[0]), _dynamic_array_INMDA_PYso_IN_N_outgoing.size()*sizeof(_dynamic_array_INMDA_PYso_IN_N_outgoing[0]));
            outfile__dynamic_array_INMDA_PYso_IN_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_Npost;
    outfile__dynamic_array_INMDA_PYso_IN_Npost.open(results_dir + "_dynamic_array_INMDA_PYso_IN_Npost_3164335675", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_Npost.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_Npost.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_Npost.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_Npost[0]), _dynamic_array_INMDA_PYso_IN_Npost.size()*sizeof(_dynamic_array_INMDA_PYso_IN_Npost[0]));
            outfile__dynamic_array_INMDA_PYso_IN_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_Npost." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_Npre;
    outfile__dynamic_array_INMDA_PYso_IN_Npre.open(results_dir + "_dynamic_array_INMDA_PYso_IN_Npre_4190079150", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_Npre.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_Npre.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_Npre.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_Npre[0]), _dynamic_array_INMDA_PYso_IN_Npre.size()*sizeof(_dynamic_array_INMDA_PYso_IN_Npre[0]));
            outfile__dynamic_array_INMDA_PYso_IN_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_Npre." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_radius;
    outfile__dynamic_array_INMDA_PYso_IN_radius.open(results_dir + "_dynamic_array_INMDA_PYso_IN_radius_1529500867", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_radius.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_radius.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_radius.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_radius[0]), _dynamic_array_INMDA_PYso_IN_radius.size()*sizeof(_dynamic_array_INMDA_PYso_IN_radius[0]));
            outfile__dynamic_array_INMDA_PYso_IN_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_radius." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_remove_recurrent_bool;
    outfile__dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.open(results_dir + "_dynamic_array_INMDA_PYso_IN_remove_recurrent_bool_3815248737", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_remove_recurrent_bool[0]), _dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.size()*sizeof(_dynamic_array_INMDA_PYso_IN_remove_recurrent_bool[0]));
            outfile__dynamic_array_INMDA_PYso_IN_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_res_NMDA;
    outfile__dynamic_array_INMDA_PYso_IN_res_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_res_NMDA_2514767774", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_res_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_res_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_res_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_res_NMDA[0]), _dynamic_array_INMDA_PYso_IN_res_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_res_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_res_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_res_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_sNMDA;
    outfile__dynamic_array_INMDA_PYso_IN_sNMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_sNMDA_144960263", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_sNMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_sNMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_sNMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_sNMDA[0]), _dynamic_array_INMDA_PYso_IN_sNMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_sNMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_sNMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_sNMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_spike_activity;
    outfile__dynamic_array_INMDA_PYso_IN_spike_activity.open(results_dir + "_dynamic_array_INMDA_PYso_IN_spike_activity_386239225", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_spike_activity.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_spike_activity.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_spike_activity.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_spike_activity[0]), _dynamic_array_INMDA_PYso_IN_spike_activity.size()*sizeof(_dynamic_array_INMDA_PYso_IN_spike_activity[0]));
            outfile__dynamic_array_INMDA_PYso_IN_spike_activity.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_spike_activity." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_tauRes_NMDA;
    outfile__dynamic_array_INMDA_PYso_IN_tauRes_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_tauRes_NMDA_2149476190", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_tauRes_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_tauRes_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_tauRes_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_tauRes_NMDA[0]), _dynamic_array_INMDA_PYso_IN_tauRes_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_tauRes_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_tauRes_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_tauRes_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_tauS_NMDA;
    outfile__dynamic_array_INMDA_PYso_IN_tauS_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_tauS_NMDA_438999200", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_tauS_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_tauS_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_tauS_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_tauS_NMDA[0]), _dynamic_array_INMDA_PYso_IN_tauS_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_tauS_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_tauS_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_tauS_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_tauX_NMDA;
    outfile__dynamic_array_INMDA_PYso_IN_tauX_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_tauX_NMDA_1894606947", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_tauX_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_tauX_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_tauX_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_tauX_NMDA[0]), _dynamic_array_INMDA_PYso_IN_tauX_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_tauX_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_tauX_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_tauX_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_IN_xNMDA;
    outfile__dynamic_array_INMDA_PYso_IN_xNMDA.open(results_dir + "_dynamic_array_INMDA_PYso_IN_xNMDA_2138298902", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_IN_xNMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_IN_xNMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_IN_xNMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_IN_xNMDA[0]), _dynamic_array_INMDA_PYso_IN_xNMDA.size()*sizeof(_dynamic_array_INMDA_PYso_IN_xNMDA[0]));
            outfile__dynamic_array_INMDA_PYso_IN_xNMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_IN_xNMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_post;
    outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_post.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr__synaptic_post_75050943", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_post.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr__synaptic_post.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr__synaptic_post[0]), _dynamic_array_INMDA_PYso_PYdr__synaptic_post.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr__synaptic_post[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_pre;
    outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_pre.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr__synaptic_pre_2799371161", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_pre.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr__synaptic_pre.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr__synaptic_pre[0]), _dynamic_array_INMDA_PYso_PYdr__synaptic_pre.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr__synaptic_pre[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA_3939841381", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA[0]), _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_alphaS_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA_2148713894", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA[0]), _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_alphaX_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_delay;
    outfile__dynamic_array_INMDA_PYso_PYdr_delay.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_delay_771056819", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_delay.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_delay.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_delay.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_delay[0]), _dynamic_array_INMDA_PYso_PYdr_delay.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_delay[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_delay." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_delay_1;
    outfile__dynamic_array_INMDA_PYso_PYdr_delay_1.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_delay_1_2819827345", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_delay_1.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_delay_1.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_delay_1[0]), _dynamic_array_INMDA_PYso_PYdr_delay_1.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_delay_1[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_delay_1.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_delay_1." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_deprFactor;
    outfile__dynamic_array_INMDA_PYso_PYdr_deprFactor.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_deprFactor_3162355954", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_deprFactor.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_deprFactor.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_deprFactor.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_deprFactor[0]), _dynamic_array_INMDA_PYso_PYdr_deprFactor.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_deprFactor[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_deprFactor.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_deprFactor." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_ENMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_ENMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_ENMDA_1595477654", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_ENMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_ENMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_ENMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_ENMDA[0]), _dynamic_array_INMDA_PYso_PYdr_ENMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_ENMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_ENMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_ENMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_gNMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_gNMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_gNMDA_3826807538", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_gNMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_gNMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_gNMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_gNMDA[0]), _dynamic_array_INMDA_PYso_PYdr_gNMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_gNMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_gNMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_gNMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_N_incoming;
    outfile__dynamic_array_INMDA_PYso_PYdr_N_incoming.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_N_incoming_4095804739", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_N_incoming.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_N_incoming.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_N_incoming[0]), _dynamic_array_INMDA_PYso_PYdr_N_incoming.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_N_incoming[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_N_outgoing;
    outfile__dynamic_array_INMDA_PYso_PYdr_N_outgoing.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_N_outgoing_3544020377", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_N_outgoing.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_N_outgoing.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_N_outgoing[0]), _dynamic_array_INMDA_PYso_PYdr_N_outgoing.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_N_outgoing[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_Npost;
    outfile__dynamic_array_INMDA_PYso_PYdr_Npost.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_Npost_3309323916", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_Npost.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_Npost.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_Npost.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_Npost[0]), _dynamic_array_INMDA_PYso_PYdr_Npost.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_Npost[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_Npost.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_Npost." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_Npre;
    outfile__dynamic_array_INMDA_PYso_PYdr_Npre.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_Npre_4272921511", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_Npre.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_Npre.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_Npre.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_Npre[0]), _dynamic_array_INMDA_PYso_PYdr_Npre.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_Npre[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_Npre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_Npre." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_radius;
    outfile__dynamic_array_INMDA_PYso_PYdr_radius.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_radius_240560452", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_radius.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_radius.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_radius.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_radius[0]), _dynamic_array_INMDA_PYso_PYdr_radius.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_radius[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_radius.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_radius." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool;
    outfile__dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool_1822304351", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool[0]), _dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_remove_recurrent_bool." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_res_NMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_res_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_res_NMDA_573217487", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_res_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_res_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_res_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_res_NMDA[0]), _dynamic_array_INMDA_PYso_PYdr_res_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_res_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_res_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_res_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_sNMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_sNMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_sNMDA_1903707056", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_sNMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_sNMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_sNMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_sNMDA[0]), _dynamic_array_INMDA_PYso_PYdr_sNMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_sNMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_sNMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_sNMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_spike_activity;
    outfile__dynamic_array_INMDA_PYso_PYdr_spike_activity.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_spike_activity_3342979106", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_spike_activity.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_spike_activity.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_spike_activity.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_spike_activity[0]), _dynamic_array_INMDA_PYso_PYdr_spike_activity.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_spike_activity[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_spike_activity.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_spike_activity." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA_1844706691", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA[0]), _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_tauRes_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_tauS_NMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_tauS_NMDA_116471575", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_tauS_NMDA[0]), _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_tauS_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_tauS_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_tauS_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_tauX_NMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_tauX_NMDA_1815529428", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_tauX_NMDA[0]), _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_tauX_NMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_tauX_NMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_tauX_NMDA." << endl;
    }
    ofstream outfile__dynamic_array_INMDA_PYso_PYdr_xNMDA;
    outfile__dynamic_array_INMDA_PYso_PYdr_xNMDA.open(results_dir + "_dynamic_array_INMDA_PYso_PYdr_xNMDA_111702689", ios::binary | ios::out);
    if(outfile__dynamic_array_INMDA_PYso_PYdr_xNMDA.is_open())
    {
        if (! _dynamic_array_INMDA_PYso_PYdr_xNMDA.empty() )
        {
            outfile__dynamic_array_INMDA_PYso_PYdr_xNMDA.write(reinterpret_cast<char*>(&_dynamic_array_INMDA_PYso_PYdr_xNMDA[0]), _dynamic_array_INMDA_PYso_PYdr_xNMDA.size()*sizeof(_dynamic_array_INMDA_PYso_PYdr_xNMDA[0]));
            outfile__dynamic_array_INMDA_PYso_PYdr_xNMDA.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_INMDA_PYso_PYdr_xNMDA." << endl;
    }
    ofstream outfile__dynamic_array_PYdr_spikemon_i;
    outfile__dynamic_array_PYdr_spikemon_i.open(results_dir + "_dynamic_array_PYdr_spikemon_i_1085799997", ios::binary | ios::out);
    if(outfile__dynamic_array_PYdr_spikemon_i.is_open())
    {
        if (! _dynamic_array_PYdr_spikemon_i.empty() )
        {
            outfile__dynamic_array_PYdr_spikemon_i.write(reinterpret_cast<char*>(&_dynamic_array_PYdr_spikemon_i[0]), _dynamic_array_PYdr_spikemon_i.size()*sizeof(_dynamic_array_PYdr_spikemon_i[0]));
            outfile__dynamic_array_PYdr_spikemon_i.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_PYdr_spikemon_i." << endl;
    }
    ofstream outfile__dynamic_array_PYdr_spikemon_t;
    outfile__dynamic_array_PYdr_spikemon_t.open(results_dir + "_dynamic_array_PYdr_spikemon_t_598840036", ios::binary | ios::out);
    if(outfile__dynamic_array_PYdr_spikemon_t.is_open())
    {
        if (! _dynamic_array_PYdr_spikemon_t.empty() )
        {
            outfile__dynamic_array_PYdr_spikemon_t.write(reinterpret_cast<char*>(&_dynamic_array_PYdr_spikemon_t[0]), _dynamic_array_PYdr_spikemon_t.size()*sizeof(_dynamic_array_PYdr_spikemon_t[0]));
            outfile__dynamic_array_PYdr_spikemon_t.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_PYdr_spikemon_t." << endl;
    }
    ofstream outfile__dynamic_array_PYso_spikemon_i;
    outfile__dynamic_array_PYso_spikemon_i.open(results_dir + "_dynamic_array_PYso_spikemon_i_460632749", ios::binary | ios::out);
    if(outfile__dynamic_array_PYso_spikemon_i.is_open())
    {
        if (! _dynamic_array_PYso_spikemon_i.empty() )
        {
            outfile__dynamic_array_PYso_spikemon_i.write(reinterpret_cast<char*>(&_dynamic_array_PYso_spikemon_i[0]), _dynamic_array_PYso_spikemon_i.size()*sizeof(_dynamic_array_PYso_spikemon_i[0]));
            outfile__dynamic_array_PYso_spikemon_i.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_PYso_spikemon_i." << endl;
    }
    ofstream outfile__dynamic_array_PYso_spikemon_t;
    outfile__dynamic_array_PYso_spikemon_t.open(results_dir + "_dynamic_array_PYso_spikemon_t_2020793972", ios::binary | ios::out);
    if(outfile__dynamic_array_PYso_spikemon_t.is_open())
    {
        if (! _dynamic_array_PYso_spikemon_t.empty() )
        {
            outfile__dynamic_array_PYso_spikemon_t.write(reinterpret_cast<char*>(&_dynamic_array_PYso_spikemon_t[0]), _dynamic_array_PYso_spikemon_t.size()*sizeof(_dynamic_array_PYso_spikemon_t[0]));
            outfile__dynamic_array_PYso_spikemon_t.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_PYso_spikemon_t." << endl;
    }
    ofstream outfile__dynamic_array_RE_spikemon_i;
    outfile__dynamic_array_RE_spikemon_i.open(results_dir + "_dynamic_array_RE_spikemon_i_3246246817", ios::binary | ios::out);
    if(outfile__dynamic_array_RE_spikemon_i.is_open())
    {
        if (! _dynamic_array_RE_spikemon_i.empty() )
        {
            outfile__dynamic_array_RE_spikemon_i.write(reinterpret_cast<char*>(&_dynamic_array_RE_spikemon_i[0]), _dynamic_array_RE_spikemon_i.size()*sizeof(_dynamic_array_RE_spikemon_i[0]));
            outfile__dynamic_array_RE_spikemon_i.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_RE_spikemon_i." << endl;
    }
    ofstream outfile__dynamic_array_RE_spikemon_t;
    outfile__dynamic_array_RE_spikemon_t.open(results_dir + "_dynamic_array_RE_spikemon_t_2726012792", ios::binary | ios::out);
    if(outfile__dynamic_array_RE_spikemon_t.is_open())
    {
        if (! _dynamic_array_RE_spikemon_t.empty() )
        {
            outfile__dynamic_array_RE_spikemon_t.write(reinterpret_cast<char*>(&_dynamic_array_RE_spikemon_t[0]), _dynamic_array_RE_spikemon_t.size()*sizeof(_dynamic_array_RE_spikemon_t[0]));
            outfile__dynamic_array_RE_spikemon_t.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_RE_spikemon_t." << endl;
    }
    ofstream outfile__dynamic_array_synapses_1__synaptic_post;
    outfile__dynamic_array_synapses_1__synaptic_post.open(results_dir + "_dynamic_array_synapses_1__synaptic_post_1999337987", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_1__synaptic_post.is_open())
    {
        if (! _dynamic_array_synapses_1__synaptic_post.empty() )
        {
            outfile__dynamic_array_synapses_1__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_synapses_1__synaptic_post[0]), _dynamic_array_synapses_1__synaptic_post.size()*sizeof(_dynamic_array_synapses_1__synaptic_post[0]));
            outfile__dynamic_array_synapses_1__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_1__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_synapses_1__synaptic_pre;
    outfile__dynamic_array_synapses_1__synaptic_pre.open(results_dir + "_dynamic_array_synapses_1__synaptic_pre_681065502", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_1__synaptic_pre.is_open())
    {
        if (! _dynamic_array_synapses_1__synaptic_pre.empty() )
        {
            outfile__dynamic_array_synapses_1__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_synapses_1__synaptic_pre[0]), _dynamic_array_synapses_1__synaptic_pre.size()*sizeof(_dynamic_array_synapses_1__synaptic_pre[0]));
            outfile__dynamic_array_synapses_1__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_1__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_synapses_1_delay;
    outfile__dynamic_array_synapses_1_delay.open(results_dir + "_dynamic_array_synapses_1_delay_2373823482", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_1_delay.is_open())
    {
        if (! _dynamic_array_synapses_1_delay.empty() )
        {
            outfile__dynamic_array_synapses_1_delay.write(reinterpret_cast<char*>(&_dynamic_array_synapses_1_delay[0]), _dynamic_array_synapses_1_delay.size()*sizeof(_dynamic_array_synapses_1_delay[0]));
            outfile__dynamic_array_synapses_1_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_1_delay." << endl;
    }
    ofstream outfile__dynamic_array_synapses_1_N_incoming;
    outfile__dynamic_array_synapses_1_N_incoming.open(results_dir + "_dynamic_array_synapses_1_N_incoming_3469555706", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_1_N_incoming.is_open())
    {
        if (! _dynamic_array_synapses_1_N_incoming.empty() )
        {
            outfile__dynamic_array_synapses_1_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_synapses_1_N_incoming[0]), _dynamic_array_synapses_1_N_incoming.size()*sizeof(_dynamic_array_synapses_1_N_incoming[0]));
            outfile__dynamic_array_synapses_1_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_1_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_synapses_1_N_outgoing;
    outfile__dynamic_array_synapses_1_N_outgoing.open(results_dir + "_dynamic_array_synapses_1_N_outgoing_3922806560", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_1_N_outgoing.is_open())
    {
        if (! _dynamic_array_synapses_1_N_outgoing.empty() )
        {
            outfile__dynamic_array_synapses_1_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_synapses_1_N_outgoing[0]), _dynamic_array_synapses_1_N_outgoing.size()*sizeof(_dynamic_array_synapses_1_N_outgoing[0]));
            outfile__dynamic_array_synapses_1_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_1_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_synapses_2__synaptic_post;
    outfile__dynamic_array_synapses_2__synaptic_post.open(results_dir + "_dynamic_array_synapses_2__synaptic_post_1591987953", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_2__synaptic_post.is_open())
    {
        if (! _dynamic_array_synapses_2__synaptic_post.empty() )
        {
            outfile__dynamic_array_synapses_2__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_synapses_2__synaptic_post[0]), _dynamic_array_synapses_2__synaptic_post.size()*sizeof(_dynamic_array_synapses_2__synaptic_post[0]));
            outfile__dynamic_array_synapses_2__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_2__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_synapses_2__synaptic_pre;
    outfile__dynamic_array_synapses_2__synaptic_pre.open(results_dir + "_dynamic_array_synapses_2__synaptic_pre_971331175", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_2__synaptic_pre.is_open())
    {
        if (! _dynamic_array_synapses_2__synaptic_pre.empty() )
        {
            outfile__dynamic_array_synapses_2__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_synapses_2__synaptic_pre[0]), _dynamic_array_synapses_2__synaptic_pre.size()*sizeof(_dynamic_array_synapses_2__synaptic_pre[0]));
            outfile__dynamic_array_synapses_2__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_2__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_synapses_2_delay;
    outfile__dynamic_array_synapses_2_delay.open(results_dir + "_dynamic_array_synapses_2_delay_3163926887", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_2_delay.is_open())
    {
        if (! _dynamic_array_synapses_2_delay.empty() )
        {
            outfile__dynamic_array_synapses_2_delay.write(reinterpret_cast<char*>(&_dynamic_array_synapses_2_delay[0]), _dynamic_array_synapses_2_delay.size()*sizeof(_dynamic_array_synapses_2_delay[0]));
            outfile__dynamic_array_synapses_2_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_2_delay." << endl;
    }
    ofstream outfile__dynamic_array_synapses_2_N_incoming;
    outfile__dynamic_array_synapses_2_N_incoming.open(results_dir + "_dynamic_array_synapses_2_N_incoming_3109283082", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_2_N_incoming.is_open())
    {
        if (! _dynamic_array_synapses_2_N_incoming.empty() )
        {
            outfile__dynamic_array_synapses_2_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_synapses_2_N_incoming[0]), _dynamic_array_synapses_2_N_incoming.size()*sizeof(_dynamic_array_synapses_2_N_incoming[0]));
            outfile__dynamic_array_synapses_2_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_2_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_synapses_2_N_outgoing;
    outfile__dynamic_array_synapses_2_N_outgoing.open(results_dir + "_dynamic_array_synapses_2_N_outgoing_2656015824", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_2_N_outgoing.is_open())
    {
        if (! _dynamic_array_synapses_2_N_outgoing.empty() )
        {
            outfile__dynamic_array_synapses_2_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_synapses_2_N_outgoing[0]), _dynamic_array_synapses_2_N_outgoing.size()*sizeof(_dynamic_array_synapses_2_N_outgoing[0]));
            outfile__dynamic_array_synapses_2_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_2_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_synapses__synaptic_post;
    outfile__dynamic_array_synapses__synaptic_post.open(results_dir + "_dynamic_array_synapses__synaptic_post_1801389495", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses__synaptic_post.is_open())
    {
        if (! _dynamic_array_synapses__synaptic_post.empty() )
        {
            outfile__dynamic_array_synapses__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_synapses__synaptic_post[0]), _dynamic_array_synapses__synaptic_post.size()*sizeof(_dynamic_array_synapses__synaptic_post[0]));
            outfile__dynamic_array_synapses__synaptic_post.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses__synaptic_post." << endl;
    }
    ofstream outfile__dynamic_array_synapses__synaptic_pre;
    outfile__dynamic_array_synapses__synaptic_pre.open(results_dir + "_dynamic_array_synapses__synaptic_pre_814148175", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses__synaptic_pre.is_open())
    {
        if (! _dynamic_array_synapses__synaptic_pre.empty() )
        {
            outfile__dynamic_array_synapses__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_synapses__synaptic_pre[0]), _dynamic_array_synapses__synaptic_pre.size()*sizeof(_dynamic_array_synapses__synaptic_pre[0]));
            outfile__dynamic_array_synapses__synaptic_pre.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses__synaptic_pre." << endl;
    }
    ofstream outfile__dynamic_array_synapses_delay;
    outfile__dynamic_array_synapses_delay.open(results_dir + "_dynamic_array_synapses_delay_3246960869", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_delay.is_open())
    {
        if (! _dynamic_array_synapses_delay.empty() )
        {
            outfile__dynamic_array_synapses_delay.write(reinterpret_cast<char*>(&_dynamic_array_synapses_delay[0]), _dynamic_array_synapses_delay.size()*sizeof(_dynamic_array_synapses_delay[0]));
            outfile__dynamic_array_synapses_delay.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_delay." << endl;
    }
    ofstream outfile__dynamic_array_synapses_N_incoming;
    outfile__dynamic_array_synapses_N_incoming.open(results_dir + "_dynamic_array_synapses_N_incoming_1151751685", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_N_incoming.is_open())
    {
        if (! _dynamic_array_synapses_N_incoming.empty() )
        {
            outfile__dynamic_array_synapses_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_synapses_N_incoming[0]), _dynamic_array_synapses_N_incoming.size()*sizeof(_dynamic_array_synapses_N_incoming[0]));
            outfile__dynamic_array_synapses_N_incoming.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_N_incoming." << endl;
    }
    ofstream outfile__dynamic_array_synapses_N_outgoing;
    outfile__dynamic_array_synapses_N_outgoing.open(results_dir + "_dynamic_array_synapses_N_outgoing_1673144031", ios::binary | ios::out);
    if(outfile__dynamic_array_synapses_N_outgoing.is_open())
    {
        if (! _dynamic_array_synapses_N_outgoing.empty() )
        {
            outfile__dynamic_array_synapses_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_synapses_N_outgoing[0]), _dynamic_array_synapses_N_outgoing.size()*sizeof(_dynamic_array_synapses_N_outgoing[0]));
            outfile__dynamic_array_synapses_N_outgoing.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_synapses_N_outgoing." << endl;
    }
    ofstream outfile__dynamic_array_TC_spikemon_i;
    outfile__dynamic_array_TC_spikemon_i.open(results_dir + "_dynamic_array_TC_spikemon_i_2393223580", ios::binary | ios::out);
    if(outfile__dynamic_array_TC_spikemon_i.is_open())
    {
        if (! _dynamic_array_TC_spikemon_i.empty() )
        {
            outfile__dynamic_array_TC_spikemon_i.write(reinterpret_cast<char*>(&_dynamic_array_TC_spikemon_i[0]), _dynamic_array_TC_spikemon_i.size()*sizeof(_dynamic_array_TC_spikemon_i[0]));
            outfile__dynamic_array_TC_spikemon_i.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_TC_spikemon_i." << endl;
    }
    ofstream outfile__dynamic_array_TC_spikemon_t;
    outfile__dynamic_array_TC_spikemon_t.open(results_dir + "_dynamic_array_TC_spikemon_t_3986939205", ios::binary | ios::out);
    if(outfile__dynamic_array_TC_spikemon_t.is_open())
    {
        if (! _dynamic_array_TC_spikemon_t.empty() )
        {
            outfile__dynamic_array_TC_spikemon_t.write(reinterpret_cast<char*>(&_dynamic_array_TC_spikemon_t[0]), _dynamic_array_TC_spikemon_t.size()*sizeof(_dynamic_array_TC_spikemon_t[0]));
            outfile__dynamic_array_TC_spikemon_t.close();
        }
    } else
    {
        std::cout << "Error writing output file for _dynamic_array_TC_spikemon_t." << endl;
    }


    // Write spike queue states to disk
    ofstream outfile_IAMPA_PYso_IN_down;
    outfile_IAMPA_PYso_IN_down.open(results_dir + "IAMPA_PYso_IN_down_queue", ios::out);
    if (outfile_IAMPA_PYso_IN_down.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IAMPA_PYso_IN_down << *IAMPA_PYso_IN_down.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IAMPA_PYso_IN_down' for file" << std::endl;
    }
    ofstream outfile_IAMPA_PYso_IN_up;
    outfile_IAMPA_PYso_IN_up.open(results_dir + "IAMPA_PYso_IN_up_queue", ios::out);
    if (outfile_IAMPA_PYso_IN_up.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IAMPA_PYso_IN_up << *IAMPA_PYso_IN_up.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IAMPA_PYso_IN_up' for file" << std::endl;
    }
    ofstream outfile_IAMPA_PYso_PYdr_down;
    outfile_IAMPA_PYso_PYdr_down.open(results_dir + "IAMPA_PYso_PYdr_down_queue", ios::out);
    if (outfile_IAMPA_PYso_PYdr_down.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IAMPA_PYso_PYdr_down << *IAMPA_PYso_PYdr_down.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IAMPA_PYso_PYdr_down' for file" << std::endl;
    }
    ofstream outfile_IAMPA_PYso_PYdr_up;
    outfile_IAMPA_PYso_PYdr_up.open(results_dir + "IAMPA_PYso_PYdr_up_queue", ios::out);
    if (outfile_IAMPA_PYso_PYdr_up.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IAMPA_PYso_PYdr_up << *IAMPA_PYso_PYdr_up.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IAMPA_PYso_PYdr_up' for file" << std::endl;
    }
    ofstream outfile_IGABAA_IN_IN_down;
    outfile_IGABAA_IN_IN_down.open(results_dir + "IGABAA_IN_IN_down_queue", ios::out);
    if (outfile_IGABAA_IN_IN_down.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IGABAA_IN_IN_down << *IGABAA_IN_IN_down.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IGABAA_IN_IN_down' for file" << std::endl;
    }
    ofstream outfile_IGABAA_IN_IN_up;
    outfile_IGABAA_IN_IN_up.open(results_dir + "IGABAA_IN_IN_up_queue", ios::out);
    if (outfile_IGABAA_IN_IN_up.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IGABAA_IN_IN_up << *IGABAA_IN_IN_up.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IGABAA_IN_IN_up' for file" << std::endl;
    }
    ofstream outfile_IGABAA_IN_PYso_down;
    outfile_IGABAA_IN_PYso_down.open(results_dir + "IGABAA_IN_PYso_down_queue", ios::out);
    if (outfile_IGABAA_IN_PYso_down.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IGABAA_IN_PYso_down << *IGABAA_IN_PYso_down.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IGABAA_IN_PYso_down' for file" << std::endl;
    }
    ofstream outfile_IGABAA_IN_PYso_up;
    outfile_IGABAA_IN_PYso_up.open(results_dir + "IGABAA_IN_PYso_up_queue", ios::out);
    if (outfile_IGABAA_IN_PYso_up.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_IGABAA_IN_PYso_up << *IGABAA_IN_PYso_up.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'IGABAA_IN_PYso_up' for file" << std::endl;
    }
    ofstream outfile_INMDA_PYso_IN_down;
    outfile_INMDA_PYso_IN_down.open(results_dir + "INMDA_PYso_IN_down_queue", ios::out);
    if (outfile_INMDA_PYso_IN_down.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_INMDA_PYso_IN_down << *INMDA_PYso_IN_down.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'INMDA_PYso_IN_down' for file" << std::endl;
    }
    ofstream outfile_INMDA_PYso_IN_up;
    outfile_INMDA_PYso_IN_up.open(results_dir + "INMDA_PYso_IN_up_queue", ios::out);
    if (outfile_INMDA_PYso_IN_up.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_INMDA_PYso_IN_up << *INMDA_PYso_IN_up.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'INMDA_PYso_IN_up' for file" << std::endl;
    }
    ofstream outfile_INMDA_PYso_PYdr_down;
    outfile_INMDA_PYso_PYdr_down.open(results_dir + "INMDA_PYso_PYdr_down_queue", ios::out);
    if (outfile_INMDA_PYso_PYdr_down.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_INMDA_PYso_PYdr_down << *INMDA_PYso_PYdr_down.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'INMDA_PYso_PYdr_down' for file" << std::endl;
    }
    ofstream outfile_INMDA_PYso_PYdr_up;
    outfile_INMDA_PYso_PYdr_up.open(results_dir + "INMDA_PYso_PYdr_up_queue", ios::out);
    if (outfile_INMDA_PYso_PYdr_up.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_INMDA_PYso_PYdr_up << *INMDA_PYso_PYdr_up.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'INMDA_PYso_PYdr_up' for file" << std::endl;
    }
    ofstream outfile_synapses_pre;
    outfile_synapses_pre.open(results_dir + "synapses_pre_queue", ios::out);
    if (outfile_synapses_pre.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_synapses_pre << *synapses_pre.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'synapses_pre' for file" << std::endl;
    }
    ofstream outfile_synapses_1_pre;
    outfile_synapses_1_pre.open(results_dir + "synapses_1_pre_queue", ios::out);
    if (outfile_synapses_1_pre.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_synapses_1_pre << *synapses_1_pre.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'synapses_1_pre' for file" << std::endl;
    }
    ofstream outfile_synapses_2_pre;
    outfile_synapses_2_pre.open(results_dir + "synapses_2_pre_queue", ios::out);
    if (outfile_synapses_2_pre.is_open()) {
        for (int i=0; i<2; i++) {
            outfile_synapses_2_pre << *synapses_2_pre.queue[i] << "\n";
        }
    } else {
        std::cout << "Error writing spike queue state for 'synapses_2_pre' for file" << std::endl;
    }

    // Write random generator state to disk
    ofstream random_generator_state;
    random_generator_state.open(results_dir + "random_generator_state", ios::out);
    if (random_generator_state.is_open()) {
        for (int i=0; i<2; i++)
            random_generator_state << _random_generators[i] << "\n";
    } else {
        std::cout << "Error writing random generator state to file." << std::endl;
    }

    // Write last run info to disk
    ofstream outfile_last_run_info;
    outfile_last_run_info.open(results_dir + "last_run_info.txt", ios::out);
    if(outfile_last_run_info.is_open())
    {
        outfile_last_run_info << (Network::_last_run_time) << " " << (Network::_last_run_completed_fraction) << std::endl;
        outfile_last_run_info.close();
    } else
    {
        std::cout << "Error writing last run info to file." << std::endl;
    }
}

void _dealloc_arrays()
{
    using namespace brian;


    // static arrays
    if(_static_array__array_IAMPA_PYso_IN_sources!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_IN_sources;
        _static_array__array_IAMPA_PYso_IN_sources = 0;
    }
    if(_static_array__array_IAMPA_PYso_IN_targets!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_IN_targets;
        _static_array__array_IAMPA_PYso_IN_targets = 0;
    }
    if(_static_array__array_IAMPA_PYso_PYdr_sources!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_PYdr_sources;
        _static_array__array_IAMPA_PYso_PYdr_sources = 0;
    }
    if(_static_array__array_IAMPA_PYso_PYdr_targets!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_PYdr_targets;
        _static_array__array_IAMPA_PYso_PYdr_targets = 0;
    }
    if(_static_array__array_IAMPA_PYso_RE_sources!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_RE_sources;
        _static_array__array_IAMPA_PYso_RE_sources = 0;
    }
    if(_static_array__array_IAMPA_PYso_RE_targets!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_RE_targets;
        _static_array__array_IAMPA_PYso_RE_targets = 0;
    }
    if(_static_array__array_IAMPA_PYso_TC_sources!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_TC_sources;
        _static_array__array_IAMPA_PYso_TC_sources = 0;
    }
    if(_static_array__array_IAMPA_PYso_TC_targets!=0)
    {
        delete [] _static_array__array_IAMPA_PYso_TC_targets;
        _static_array__array_IAMPA_PYso_TC_targets = 0;
    }
    if(_static_array__array_IAMPA_TC_IN_sources!=0)
    {
        delete [] _static_array__array_IAMPA_TC_IN_sources;
        _static_array__array_IAMPA_TC_IN_sources = 0;
    }
    if(_static_array__array_IAMPA_TC_IN_targets!=0)
    {
        delete [] _static_array__array_IAMPA_TC_IN_targets;
        _static_array__array_IAMPA_TC_IN_targets = 0;
    }
    if(_static_array__array_IAMPA_TC_PYdr_sources!=0)
    {
        delete [] _static_array__array_IAMPA_TC_PYdr_sources;
        _static_array__array_IAMPA_TC_PYdr_sources = 0;
    }
    if(_static_array__array_IAMPA_TC_PYdr_targets!=0)
    {
        delete [] _static_array__array_IAMPA_TC_PYdr_targets;
        _static_array__array_IAMPA_TC_PYdr_targets = 0;
    }
    if(_static_array__array_IGABAA_IN_IN_sources!=0)
    {
        delete [] _static_array__array_IGABAA_IN_IN_sources;
        _static_array__array_IGABAA_IN_IN_sources = 0;
    }
    if(_static_array__array_IGABAA_IN_IN_targets!=0)
    {
        delete [] _static_array__array_IGABAA_IN_IN_targets;
        _static_array__array_IGABAA_IN_IN_targets = 0;
    }
    if(_static_array__array_IGABAA_IN_PYso_sources!=0)
    {
        delete [] _static_array__array_IGABAA_IN_PYso_sources;
        _static_array__array_IGABAA_IN_PYso_sources = 0;
    }
    if(_static_array__array_IGABAA_IN_PYso_targets!=0)
    {
        delete [] _static_array__array_IGABAA_IN_PYso_targets;
        _static_array__array_IGABAA_IN_PYso_targets = 0;
    }
    if(_static_array__array_INMDA_PYso_IN_sources!=0)
    {
        delete [] _static_array__array_INMDA_PYso_IN_sources;
        _static_array__array_INMDA_PYso_IN_sources = 0;
    }
    if(_static_array__array_INMDA_PYso_IN_targets!=0)
    {
        delete [] _static_array__array_INMDA_PYso_IN_targets;
        _static_array__array_INMDA_PYso_IN_targets = 0;
    }
    if(_static_array__array_INMDA_PYso_PYdr_sources!=0)
    {
        delete [] _static_array__array_INMDA_PYso_PYdr_sources;
        _static_array__array_INMDA_PYso_PYdr_sources = 0;
    }
    if(_static_array__array_INMDA_PYso_PYdr_targets!=0)
    {
        delete [] _static_array__array_INMDA_PYso_PYdr_targets;
        _static_array__array_INMDA_PYso_PYdr_targets = 0;
    }
    if(_static_array__array_IN_group_hNa_IN!=0)
    {
        delete [] _static_array__array_IN_group_hNa_IN;
        _static_array__array_IN_group_hNa_IN = 0;
    }
    if(_static_array__array_IN_group_nK_IN!=0)
    {
        delete [] _static_array__array_IN_group_nK_IN;
        _static_array__array_IN_group_nK_IN = 0;
    }
    if(_static_array__array_IN_group_v!=0)
    {
        delete [] _static_array__array_IN_group_v;
        _static_array__array_IN_group_v = 0;
    }
    if(_static_array__array_PYdr_group_CaBuffer!=0)
    {
        delete [] _static_array__array_PYdr_group_CaBuffer;
        _static_array__array_PYdr_group_CaBuffer = 0;
    }
    if(_static_array__array_PYdr_group_v!=0)
    {
        delete [] _static_array__array_PYdr_group_v;
        _static_array__array_PYdr_group_v = 0;
    }
    if(_static_array__array_PYso_group_hA!=0)
    {
        delete [] _static_array__array_PYso_group_hA;
        _static_array__array_PYso_group_hA = 0;
    }
    if(_static_array__array_PYso_group_hNa!=0)
    {
        delete [] _static_array__array_PYso_group_hNa;
        _static_array__array_PYso_group_hNa = 0;
    }
    if(_static_array__array_PYso_group_mKS!=0)
    {
        delete [] _static_array__array_PYso_group_mKS;
        _static_array__array_PYso_group_mKS = 0;
    }
    if(_static_array__array_PYso_group_nK!=0)
    {
        delete [] _static_array__array_PYso_group_nK;
        _static_array__array_PYso_group_nK = 0;
    }
    if(_static_array__array_PYso_group_v!=0)
    {
        delete [] _static_array__array_PYso_group_v;
        _static_array__array_PYso_group_v = 0;
    }
    if(_static_array__array_RE_group_hNa_RE!=0)
    {
        delete [] _static_array__array_RE_group_hNa_RE;
        _static_array__array_RE_group_hNa_RE = 0;
    }
    if(_static_array__array_RE_group_hT_RE!=0)
    {
        delete [] _static_array__array_RE_group_hT_RE;
        _static_array__array_RE_group_hT_RE = 0;
    }
    if(_static_array__array_RE_group_mNa_RE!=0)
    {
        delete [] _static_array__array_RE_group_mNa_RE;
        _static_array__array_RE_group_mNa_RE = 0;
    }
    if(_static_array__array_RE_group_mT_RE!=0)
    {
        delete [] _static_array__array_RE_group_mT_RE;
        _static_array__array_RE_group_mT_RE = 0;
    }
    if(_static_array__array_RE_group_nK_RE!=0)
    {
        delete [] _static_array__array_RE_group_nK_RE;
        _static_array__array_RE_group_nK_RE = 0;
    }
    if(_static_array__array_RE_group_v!=0)
    {
        delete [] _static_array__array_RE_group_v;
        _static_array__array_RE_group_v = 0;
    }
    if(_static_array__array_TC_group_Ca_TC!=0)
    {
        delete [] _static_array__array_TC_group_Ca_TC;
        _static_array__array_TC_group_Ca_TC = 0;
    }
    if(_static_array__array_TC_group_Open!=0)
    {
        delete [] _static_array__array_TC_group_Open;
        _static_array__array_TC_group_Open = 0;
    }
    if(_static_array__array_TC_group_OpenLocked!=0)
    {
        delete [] _static_array__array_TC_group_OpenLocked;
        _static_array__array_TC_group_OpenLocked = 0;
    }
    if(_static_array__array_TC_group_Pone!=0)
    {
        delete [] _static_array__array_TC_group_Pone;
        _static_array__array_TC_group_Pone = 0;
    }
    if(_static_array__array_TC_group_hNa_TC!=0)
    {
        delete [] _static_array__array_TC_group_hNa_TC;
        _static_array__array_TC_group_hNa_TC = 0;
    }
    if(_static_array__array_TC_group_hT_TC!=0)
    {
        delete [] _static_array__array_TC_group_hT_TC;
        _static_array__array_TC_group_hT_TC = 0;
    }
    if(_static_array__array_TC_group_mNa_TC!=0)
    {
        delete [] _static_array__array_TC_group_mNa_TC;
        _static_array__array_TC_group_mNa_TC = 0;
    }
    if(_static_array__array_TC_group_nK_TC!=0)
    {
        delete [] _static_array__array_TC_group_nK_TC;
        _static_array__array_TC_group_nK_TC = 0;
    }
    if(_static_array__array_TC_group_v!=0)
    {
        delete [] _static_array__array_TC_group_v;
        _static_array__array_TC_group_v = 0;
    }
    if(_static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA;
        _static_array__dynamic_array_IAMPA_PYso_IN_res_AMPA = 0;
    }
    if(_static_array__dynamic_array_IAMPA_PYso_IN_sAMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_PYso_IN_sAMPA;
        _static_array__dynamic_array_IAMPA_PYso_IN_sAMPA = 0;
    }
    if(_static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA;
        _static_array__dynamic_array_IAMPA_PYso_PYdr_res_AMPA = 0;
    }
    if(_static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA;
        _static_array__dynamic_array_IAMPA_PYso_PYdr_sAMPA = 0;
    }
    if(_static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE;
        _static_array__dynamic_array_IAMPA_PYso_RE_sAMPA_PYso_RE = 0;
    }
    if(_static_array__dynamic_array_IAMPA_PYso_TC_sAMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_PYso_TC_sAMPA;
        _static_array__dynamic_array_IAMPA_PYso_TC_sAMPA = 0;
    }
    if(_static_array__dynamic_array_IAMPA_TC_IN_sAMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_TC_IN_sAMPA;
        _static_array__dynamic_array_IAMPA_TC_IN_sAMPA = 0;
    }
    if(_static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA;
        _static_array__dynamic_array_IAMPA_TC_PYdr_sAMPA = 0;
    }
    if(_static_array__dynamic_array_IAMPA_TC_RE_sAMPA!=0)
    {
        delete [] _static_array__dynamic_array_IAMPA_TC_RE_sAMPA;
        _static_array__dynamic_array_IAMPA_TC_RE_sAMPA = 0;
    }
    if(_static_array__dynamic_array_ICOM_PYdr_PYso_gCOM!=0)
    {
        delete [] _static_array__dynamic_array_ICOM_PYdr_PYso_gCOM;
        _static_array__dynamic_array_ICOM_PYdr_PYso_gCOM = 0;
    }
    if(_static_array__dynamic_array_ICOM_PYso_PYdr_gCOM!=0)
    {
        delete [] _static_array__dynamic_array_ICOM_PYso_PYdr_gCOM;
        _static_array__dynamic_array_ICOM_PYso_PYdr_gCOM = 0;
    }
    if(_static_array__dynamic_array_IGABAA_IN_IN_res_GABAA!=0)
    {
        delete [] _static_array__dynamic_array_IGABAA_IN_IN_res_GABAA;
        _static_array__dynamic_array_IGABAA_IN_IN_res_GABAA = 0;
    }
    if(_static_array__dynamic_array_IGABAA_IN_IN_sGABAA!=0)
    {
        delete [] _static_array__dynamic_array_IGABAA_IN_IN_sGABAA;
        _static_array__dynamic_array_IGABAA_IN_IN_sGABAA = 0;
    }
    if(_static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA!=0)
    {
        delete [] _static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA;
        _static_array__dynamic_array_IGABAA_IN_PYso_res_GABAA = 0;
    }
    if(_static_array__dynamic_array_IGABAA_IN_PYso_sGABAA!=0)
    {
        delete [] _static_array__dynamic_array_IGABAA_IN_PYso_sGABAA;
        _static_array__dynamic_array_IGABAA_IN_PYso_sGABAA = 0;
    }
    if(_static_array__dynamic_array_IGABAA_RE_RE_sGABAA!=0)
    {
        delete [] _static_array__dynamic_array_IGABAA_RE_RE_sGABAA;
        _static_array__dynamic_array_IGABAA_RE_RE_sGABAA = 0;
    }
    if(_static_array__dynamic_array_IGABAA_RE_TC_sGABAA!=0)
    {
        delete [] _static_array__dynamic_array_IGABAA_RE_TC_sGABAA;
        _static_array__dynamic_array_IGABAA_RE_TC_sGABAA = 0;
    }
    if(_static_array__dynamic_array_IGABAB_RE_TC_rGABAB!=0)
    {
        delete [] _static_array__dynamic_array_IGABAB_RE_TC_rGABAB;
        _static_array__dynamic_array_IGABAB_RE_TC_rGABAB = 0;
    }
    if(_static_array__dynamic_array_IGABAB_RE_TC_sGABAB!=0)
    {
        delete [] _static_array__dynamic_array_IGABAB_RE_TC_sGABAB;
        _static_array__dynamic_array_IGABAB_RE_TC_sGABAB = 0;
    }
    if(_static_array__dynamic_array_IKNa_PYso_PYdr_concNa!=0)
    {
        delete [] _static_array__dynamic_array_IKNa_PYso_PYdr_concNa;
        _static_array__dynamic_array_IKNa_PYso_PYdr_concNa = 0;
    }
    if(_static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local!=0)
    {
        delete [] _static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local;
        _static_array__dynamic_array_IKNa_PYso_PYdr_hNa_local = 0;
    }
    if(_static_array__dynamic_array_INMDA_PYso_IN_res_NMDA!=0)
    {
        delete [] _static_array__dynamic_array_INMDA_PYso_IN_res_NMDA;
        _static_array__dynamic_array_INMDA_PYso_IN_res_NMDA = 0;
    }
    if(_static_array__dynamic_array_INMDA_PYso_IN_sNMDA!=0)
    {
        delete [] _static_array__dynamic_array_INMDA_PYso_IN_sNMDA;
        _static_array__dynamic_array_INMDA_PYso_IN_sNMDA = 0;
    }
    if(_static_array__dynamic_array_INMDA_PYso_IN_xNMDA!=0)
    {
        delete [] _static_array__dynamic_array_INMDA_PYso_IN_xNMDA;
        _static_array__dynamic_array_INMDA_PYso_IN_xNMDA = 0;
    }
    if(_static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA!=0)
    {
        delete [] _static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA;
        _static_array__dynamic_array_INMDA_PYso_PYdr_res_NMDA = 0;
    }
    if(_static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA!=0)
    {
        delete [] _static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA;
        _static_array__dynamic_array_INMDA_PYso_PYdr_sNMDA = 0;
    }
    if(_static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA!=0)
    {
        delete [] _static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA;
        _static_array__dynamic_array_INMDA_PYso_PYdr_xNMDA = 0;
    }
}

