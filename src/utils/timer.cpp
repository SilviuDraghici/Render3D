#include "timer.h"

#include <iostream>

//std::string timer_name;
//std::chrono::steady_clock::time_point start_t, end_t;

Timer::Timer(){}

Timer::Timer(const char *name) {
    timer_name = name;
}

Timer::Timer(const std::string& name) {
    timer_name = name;
}

void Timer::start() {
    start_t = std::chrono::steady_clock::now();
}

void Timer::end() {
    end_t = std::chrono::steady_clock::now();
}

void Timer::print_elapsed_time(std::ostream &strm) {
    if(!timer_name.empty()){
        strm << timer_name << " ";
    }

    auto time = end_t - start_t;
    auto mins = std::chrono::duration_cast<std::chrono::minutes>(time).count();
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(time).count() - mins * 60;
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(time).count() - secs * 1000;
    strm << "Time: min:" << mins << " sec:" << secs << " milli:" << milli << std::endl;
}

void Timer::print_elapsed_time(){
    print_elapsed_time(std::cout);
}

