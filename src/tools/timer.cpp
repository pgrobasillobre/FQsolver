#include "timer.hpp"
#include "output.hpp"

#include <iomanip>
#include <stdexcept>
#include <ctime>
#include <cmath>
#include <sys/resource.h>

namespace
{
double current_cpu_seconds()
{
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0)
    {
        return 0.0;
    }

    const double user_seconds = static_cast<double>(usage.ru_utime.tv_sec) +
                                static_cast<double>(usage.ru_utime.tv_usec) * 1.0e-6;
    const double system_seconds = static_cast<double>(usage.ru_stime.tv_sec) +
                                  static_cast<double>(usage.ru_stime.tv_usec) * 1.0e-6;
    return user_seconds + system_seconds;
}
}

//----------------------------------------------------------------------
// Initialize default timers (e.g., "total")
void Timer::initialize()
{
    timers["total"] = TimeData{};
}
//----------------------------------------------------------------------
// Start the timer with the given name
void Timer::start(const std::string &name)
{
    auto it = timers.find(name);
    if (it != timers.end())
    {
        it->second.start_time = std::chrono::high_resolution_clock::now();
        it->second.start_cpu_time = current_cpu_seconds();
        it->second.started = true;
    }
    else
    {
        throw std::runtime_error("Timer \"" + name + "\" is not initialized");
    }
}
//----------------------------------------------------------------------
// Finish the timer with the given name
void Timer::finish(const std::string &name)
{
    auto it = timers.find(name);
    if (it != timers.end() && it->second.started)
    {
        it->second.end_time = std::chrono::high_resolution_clock::now();
        it->second.end_cpu_time = current_cpu_seconds();
        it->second.finished = true;
    }
    else
    {
        throw std::runtime_error("Timer \"" + name + "\" was not started or initialized");
    }
}
//----------------------------------------------------------------------
// Print timer results and termination message
void Timer::conclude(const Output &out)
{
    using namespace std::chrono;

    const auto now = system_clock::now();
    const std::time_t end_time = system_clock::to_time_t(now);
    const std::tm *tm_ptr = std::localtime(&end_time);

    // Hardcoded total time only
    const std::string timer_name = "total";

    long long elapsed_seconds = 0;
    long long cpu_seconds = 0;
    if (timers.count(timer_name) && timers[timer_name].started && timers[timer_name].finished)
    {
        elapsed_seconds = duration_cast<seconds>(
                              timers[timer_name].end_time - timers[timer_name].start_time)
                              .count();
        cpu_seconds = static_cast<long long>(std::llround(timers[timer_name].end_cpu_time -
                                                          timers[timer_name].start_cpu_time));
    }

    int cpu_hours = static_cast<int>(cpu_seconds / 3600);
    int cpu_minutes = static_cast<int>((cpu_seconds % 3600) / 60);
    int cpu_secs = static_cast<int>(cpu_seconds % 60);

    int elapsed_hours = static_cast<int>(elapsed_seconds / 3600);
    int elapsed_minutes = static_cast<int>((elapsed_seconds % 3600) / 60);
    int elapsed_secs = static_cast<int>(elapsed_seconds % 60);

    // Header
    // out.stream() << " " << out.sticks << "\n\n";
    out.stream() << std::setw(37) << " " << "We should write this code in C++.\n\n";
    out.stream() << std::setw(53) << " " << "-- P. Grobas Illobre\n\n";
    out.stream() << " " << out.sticks << "\n\n";

    // Timing Info
    out.stream() << std::setw(42) << " " << "CPU Time:    "
                 << std::setw(3) << cpu_hours << " h "
                 << std::setw(2) << cpu_minutes << " min "
                 << std::setw(2) << cpu_secs << " sec\n";

    out.stream() << std::setw(42) << " " << "Elapsed Time:"
                 << std::setw(3) << elapsed_hours << " h "
                 << std::setw(2) << elapsed_minutes << " min "
                 << std::setw(2) << elapsed_secs << " sec\n";

    out.stream() << "\n " << out.sticks << "\n\n";

    // Date/Time Stamp
    out.stream() << std::setw(4) << " "
                 << "Normal Termination of FQSolver program in date "
                 << std::put_time(tm_ptr, "%d/%m/%Y at %H:%M:%S") << "\n\n";

    out.stream() << " " << out.sticks << "\n";
}
//----------------------------------------------------------------------
