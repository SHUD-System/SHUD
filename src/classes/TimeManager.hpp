#ifndef TIMEMANAGER_HPP
#define TIMEMANAGER_HPP

#include <chrono>
#include <string>

/**
 * TimeManager class - Manages all time-related functionality for the SHUD model
 * 
 * This class encapsulates time management including:
 * - Model base date (modelBaseDate) as the unified time reference
 * - Relative time variables (t0, t1, dt) for model calculations
 * - Absolute time points (timelocal, timeutc) for high-precision operations
 * - Julian day calculation
 * - Time formatting in 7 different formats
 */
class TimeManager {
public:
    TimeManager();
    ~TimeManager();
    
    // Initialization
    void setModelBaseDate(long yyyymmdd);
    bool isInitialized() const { return initialized; }
    
    // Time update
    void updateTime(double current_time_minutes);
    
    // Accessors (read-only)
    long getModelBaseDate() const { return modelBaseDate; }
    std::chrono::system_clock::time_point getBaseTime() const { return baseTime; }
    double getT0() const { return t0; }
    double getT1() const { return t1; }
    double getDT() const { return dt; }
    int getJulianDay() const { return jd; }
    std::chrono::system_clock::time_point getLocalTime() const { return timelocal; }
    std::chrono::system_clock::time_point getUTCTime() const { return timeutc; }
    
    // Unified formatting output
    std::string formatLocalTime(int type = 6, bool useSeparators = true) const;   // Default: YYYY-MM-DD hh:mm:ss
    std::string formatUTCTime(int type = 6, bool useSeparators = true) const;     // Default: YYYY-MM-DD hh:mm:ss
    std::string formatBaseTime(int type = 6, bool useSeparators = true) const;    // Default: YYYY-MM-DD hh:mm:ss
    
private:
    std::string formatTime(std::chrono::system_clock::time_point tp, int type, bool useSeparators) const;
    
    // Base time (two representations)
    long modelBaseDate;  // YYYYMMDD format, for file I/O
    std::chrono::system_clock::time_point baseTime;  // time_point format, for time calculations
    bool initialized;    // Protection flag: true after first initialization
    
    // Relative time (maintaining existing compatibility)
    double t0, t1, dt;  // minutes, relative to modelBaseDate
    
    // Absolute time (new functionality)
    std::chrono::system_clock::time_point timelocal;
    std::chrono::system_clock::time_point timeutc;
    
    // Derived time
    int jd;  // Julian day (1-366)
    
    // Helper methods
    int calculateJulianDay(long timestamp, double elapsed_days);
    std::chrono::system_clock::time_point convertYYYYMMDDtoTimePoint(long yyyymmdd);
    void calculateLocalTime(double elapsed_minutes);
    void calculateUTCTime(double elapsed_minutes);
};

#endif // TIMEMANAGER_HPP
