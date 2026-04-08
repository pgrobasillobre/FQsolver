#include "timer.hpp"
#include "output.hpp"

#include <iomanip>
#include <stdexcept>
#include <ctime>

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
    if (timers.count(timer_name) && timers[timer_name].started && timers[timer_name].finished)
    {
        elapsed_seconds = duration_cast<seconds>(
                              timers[timer_name].end_time - timers[timer_name].start_time)
                              .count();
    }

    int hours = static_cast<int>(elapsed_seconds / 3600);
    int minutes = static_cast<int>((elapsed_seconds % 3600) / 60);
    int seconds = static_cast<int>(elapsed_seconds % 60);

    // Header
    // out.stream() << " " << out.sticks << "\n\n";
    out.stream() << std::setw(37) << " " << "We should write this code in C++.\n\n";
    out.stream() << std::setw(53) << " " << "-- P. Grobas Illobre\n\n";
    out.stream() << " " << out.sticks << "\n\n";

    // Timing Info
    out.stream() << std::setw(42) << " " << "CPU Time:    "
                 << std::setw(3) << hours << " h "
                 << std::setw(2) << minutes << " min "
                 << std::setw(2) << seconds << " sec\n";

    out.stream() << std::setw(42) << " " << "Elapsed Time:"
                 << std::setw(3) << hours << " h "
                 << std::setw(2) << minutes << " min "
                 << std::setw(2) << seconds << " sec\n";

    out.stream() << "\n " << out.sticks << "\n\n";

    // Date/Time Stamp
    out.stream() << std::setw(4) << " "
                 << "Normal Termination of FQSolver program in date "
                 << std::put_time(tm_ptr, "%d/%m/%Y at %H:%M:%S") << "\n\n";

    out.stream() << " " << out.sticks << "\n";
}
//----------------------------------------------------------------------
