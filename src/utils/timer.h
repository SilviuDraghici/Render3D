#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>

class Timer {
    std::string timer_name;
    std::chrono::steady_clock::time_point start_t, end_t;

   public:
    Timer();
    Timer(const char *name);
    Timer(const std::string& name);
    
    void start();
    void end();
    
    void print_elapsed_time();
    void print_elapsed_time(std::ostream &strm);
};

#endif
