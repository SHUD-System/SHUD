#include "TimeManager.hpp"
#include <ctime>
#include <iomanip>
#include <sstream>
#include <cmath>

TimeManager::TimeManager() 
    : modelBaseDate(0), t0(0.0), t1(0.0), dt(0.0), jd(0), initialized(false) {
    // Initialize time points to epoch
    baseTime = std::chrono::system_clock::time_point();
    timelocal = baseTime;
    timeutc = baseTime;
}

TimeManager::~TimeManager() {
    // Nothing to clean up
}

void TimeManager::setModelBaseDate(long yyyymmdd) {
    // Protection: only allow setting once
    if (initialized) {
        fprintf(stderr, "WARNING: Attempt to modify modelBaseDate after initialization!\n");
        fprintf(stderr, "         Current value: %ld, Attempted value: %ld\n", modelBaseDate, yyyymmdd);
        fprintf(stderr, "         Ignoring the modification to maintain consistency.\n");
        return;
    }
    
    modelBaseDate = yyyymmdd;
    baseTime = convertYYYYMMDDtoTimePoint(yyyymmdd);
    
    // Initialize time variables
    t0 = 0.0;
    t1 = 0.0;
    dt = 0.0;
    timelocal = baseTime;
    timeutc = baseTime;
    
    // Calculate initial Julian day
    jd = calculateJulianDay(yyyymmdd, 0.0);
    
    // Mark as initialized
    initialized = true;
}

void TimeManager::updateTime(double current_time_minutes) {
    // Update relative time variables
    t0 = t1;
    t1 = current_time_minutes;
    dt = t1 - t0;
    
    // Calculate elapsed days for Julian day calculation
    double elapsed_days = t1 / 1440.0;  // 1440 minutes per day
    jd = calculateJulianDay(modelBaseDate, elapsed_days);
    
    // Calculate absolute times
    calculateLocalTime(t1);
    calculateUTCTime(t1);
}

std::string TimeManager::formatLocalTime(int type, bool useSeparators) const {
    return formatTime(timelocal, type, useSeparators);
}

std::string TimeManager::formatUTCTime(int type, bool useSeparators) const {
    return formatTime(timeutc, type, useSeparators);
}

std::string TimeManager::formatBaseTime(int type, bool useSeparators) const {
    return formatTime(baseTime, type, useSeparators);
}

std::string TimeManager::formatTime(std::chrono::system_clock::time_point tp, int type, bool useSeparators) const {
    // Convert time_point to time_t
    std::time_t tt = std::chrono::system_clock::to_time_t(tp);
    std::tm* tm = std::localtime(&tt);
    
    // Get microseconds
    auto duration = tp.time_since_epoch();
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(duration - seconds);
    
    std::ostringstream oss;
    oss << std::setfill('0');
    
    // Define separators based on useSeparators flag
    const char* dateSep = useSeparators ? "-" : "";
    const char* timeSep = useSeparators ? ":" : "";
    const char* dateTimeSep = useSeparators ? " " : "";
    const char* microSep = useSeparators ? "." : "";
    
    switch (type) {
        case 1:  // YYYY
            oss << std::setw(4) << (tm->tm_year + 1900);
            break;
            
        case 2:  // YYYY-MM or YYYYMM
            oss << std::setw(4) << (tm->tm_year + 1900) << dateSep
                << std::setw(2) << (tm->tm_mon + 1);
            break;
            
        case 3:  // YYYY-MM-DD or YYYYMMDD
            oss << std::setw(4) << (tm->tm_year + 1900) << dateSep
                << std::setw(2) << (tm->tm_mon + 1) << dateSep
                << std::setw(2) << tm->tm_mday;
            break;
            
        case 4:  // YYYY-MM-DD hh or YYYYMMDDhh
            oss << std::setw(4) << (tm->tm_year + 1900) << dateSep
                << std::setw(2) << (tm->tm_mon + 1) << dateSep
                << std::setw(2) << tm->tm_mday << dateTimeSep
                << std::setw(2) << tm->tm_hour;
            break;
            
        case 5:  // YYYY-MM-DD hh:mm or YYYYMMDDhhmm
            oss << std::setw(4) << (tm->tm_year + 1900) << dateSep
                << std::setw(2) << (tm->tm_mon + 1) << dateSep
                << std::setw(2) << tm->tm_mday << dateTimeSep
                << std::setw(2) << tm->tm_hour << timeSep
                << std::setw(2) << tm->tm_min;
            break;
            
        case 6:  // YYYY-MM-DD hh:mm:ss or YYYYMMDDhhmmss (default)
            oss << std::setw(4) << (tm->tm_year + 1900) << dateSep
                << std::setw(2) << (tm->tm_mon + 1) << dateSep
                << std::setw(2) << tm->tm_mday << dateTimeSep
                << std::setw(2) << tm->tm_hour << timeSep
                << std::setw(2) << tm->tm_min << timeSep
                << std::setw(2) << tm->tm_sec;
            break;
            
        case 7:  // YYYY-MM-DD hh:mm:ss.ssssss or YYYYMMDDhhmmssssssss
            oss << std::setw(4) << (tm->tm_year + 1900) << dateSep
                << std::setw(2) << (tm->tm_mon + 1) << dateSep
                << std::setw(2) << tm->tm_mday << dateTimeSep
                << std::setw(2) << tm->tm_hour << timeSep
                << std::setw(2) << tm->tm_min << timeSep
                << std::setw(2) << tm->tm_sec << microSep
                << std::setw(6) << microseconds.count();
            break;
            
        default:  // Default to type 6
            oss << std::setw(4) << (tm->tm_year + 1900) << dateSep
                << std::setw(2) << (tm->tm_mon + 1) << dateSep
                << std::setw(2) << tm->tm_mday << dateTimeSep
                << std::setw(2) << tm->tm_hour << timeSep
                << std::setw(2) << tm->tm_min << timeSep
                << std::setw(2) << tm->tm_sec;
            break;
    }
    
    return oss.str();
}

int TimeManager::calculateJulianDay(long timestamp, double elapsed_days) {
    // Extract year, month, day from YYYYMMDD
    int year = timestamp / 10000;
    int month = (timestamp / 100) % 100;
    int day = timestamp % 100;
    
    // Add elapsed days
    int total_days = day + static_cast<int>(std::floor(elapsed_days));
    
    // Days in each month (non-leap year)
    int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    
    // Check for leap year
    bool is_leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    if (is_leap) {
        days_in_month[1] = 29;
    }
    
    // Handle month overflow
    while (total_days > days_in_month[month - 1]) {
        total_days -= days_in_month[month - 1];
        month++;
        if (month > 12) {
            month = 1;
            year++;
            // Recalculate leap year for new year
            is_leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
            days_in_month[1] = is_leap ? 29 : 28;
        }
    }
    
    // Calculate Julian day (day of year)
    int julian_day = 0;
    for (int m = 1; m < month; m++) {
        julian_day += days_in_month[m - 1];
    }
    julian_day += total_days;
    
    return julian_day;
}

std::chrono::system_clock::time_point TimeManager::convertYYYYMMDDtoTimePoint(long yyyymmdd) {
    // Extract year, month, day
    int year = yyyymmdd / 10000;
    int month = (yyyymmdd / 100) % 100;
    int day = yyyymmdd % 100;
    
    // Create tm structure
    std::tm tm = {};
    tm.tm_year = year - 1900;  // Years since 1900
    tm.tm_mon = month - 1;     // Months since January (0-11)
    tm.tm_mday = day;
    tm.tm_hour = 0;
    tm.tm_min = 0;
    tm.tm_sec = 0;
    tm.tm_isdst = -1;  // Let mktime determine DST
    
    // Convert to time_t
    std::time_t tt = std::mktime(&tm);
    
    // Convert to time_point
    return std::chrono::system_clock::from_time_t(tt);
}

void TimeManager::calculateLocalTime(double elapsed_minutes) {
    // Convert minutes to microseconds for high precision
    long long microseconds = static_cast<long long>(elapsed_minutes * 60.0 * 1000000.0);
    
    // Add to base time
    timelocal = baseTime + std::chrono::microseconds(microseconds);
}

void TimeManager::calculateUTCTime(double elapsed_minutes) {
    // For now, UTC time is the same as local time
    // In the future, this could be adjusted for timezone differences
    long long microseconds = static_cast<long long>(elapsed_minutes * 60.0 * 1000000.0);
    timeutc = baseTime + std::chrono::microseconds(microseconds);
}
