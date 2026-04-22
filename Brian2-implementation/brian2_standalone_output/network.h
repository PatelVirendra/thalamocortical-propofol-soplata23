
#ifndef _BRIAN_NETWORK_H
#define _BRIAN_NETWORK_H

#include<vector>
#include<utility>
#include<set>
#include "brianlib/clocks.h"

typedef void (*codeobj_func)();

class Network
{
    std::set<BaseClock*> clocks, curclocks;
    void compute_clocks();
    BaseClock* next_clocks();
public:
    std::vector< std::pair< BaseClock*, codeobj_func > > objects;
    double t;
    static double _last_run_time;
    static double _last_run_completed_fraction;
    static bool _globally_stopped;
    static bool _globally_running;

    Network();
    void clear();
    void add(BaseClock *clock, codeobj_func func);
    void run(const double duration, void (*report_func)(const double, const double, const double, const double), const double report_period);
};

#endif

