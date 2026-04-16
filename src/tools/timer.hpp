#ifndef TIMER_HPP
#define TIMER_HPP

#include "output.hpp"
#include <string>
#include <iostream>
#include <unordered_map>
#include <chrono>

//----------------------------------------------------------------------
/// @class Timer
/// @brief Handles performance timing for labeled code sections.
///
/// This class allows precise wall-time measurement of specific tasks or phases
/// within the program using high-resolution clocks. Each named timer can be
/// started, stopped, and summarized at the end of the program using `conclude()`.
class Timer {
public:
    /// @brief Initializes the default timers.
    ///
    /// Prepares internal timer map and clears any previous data.
    /// This should be called before starting any timers.
    void initialize();

    /// @brief Starts the timer for a specific task.
    ///
    /// Stores the current timestamp for the timer associated with the given name.
    /// Throws if the timer was not initialized via `initialize()`.
    ///
    /// @param name The label identifying the timer (e.g., "total", "step1").
    void start(const std::string& name);
    
    /// @brief Finishes the timer for a specific task.
    ///
    /// Stores the current timestamp as the end time for the named timer.
    /// Throws if the timer was not started first.
    ///
    /// @param name The label of the timer to finish.
    void finish(const std::string& name);

    /// @brief Finalizes and reports all timers.
    ///
    /// Calculates and prints the elapsed time for each timer to the output stream.
    /// Includes an end-of-run banner and timestamp.
    ///
    /// @param out Output object where the timer report will be written.
    void conclude(const Output& out);

private:
    /// @struct TimeData
    /// @brief Holds time-point data for a single timer.
    struct TimeData {
        std::chrono::high_resolution_clock::time_point start_time; ///< Start timestamp
        std::chrono::high_resolution_clock::time_point end_time;   ///< End timestamp
        double start_cpu_time = 0.0; ///< Start CPU time in seconds
        double end_cpu_time = 0.0;   ///< End CPU time in seconds
        bool started = false;   ///< True if the timer has been started
        bool finished = false;  ///< True if the timer has been finished
    };

    /// @brief Storage for all named timers.
    std::unordered_map<std::string, TimeData> timers;
};

#endif // TIMER_HPP
//----------------------------------------------------------------------
