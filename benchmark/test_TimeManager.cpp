#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../src/classes/TimeManager.hpp"
#include <cmath>
#include <chrono>

// Helper function: compare floating point numbers
bool almost_equal(double a, double b, double tol=1e-6) {
    return std::fabs(a - b) < tol;
}

TEST_CASE("TimeManager::setModelBaseDate initializes correctly", "[TimeManager]") {
    TimeManager tm;
    
    SECTION("Set base date 20230615") {
        tm.setModelBaseDate(20230615);
        
        REQUIRE(tm.getModelBaseDate() == 20230615);
        REQUIRE(tm.getT0() == 0.0);
        REQUIRE(tm.getT1() == 0.0);
        REQUIRE(tm.getDT() == 0.0);
        
        // Julian day for June 15 should be 166 (31+28+31+30+31+15)
        REQUIRE(tm.getJulianDay() == 166);
    }
    
    SECTION("Set base date 20240101 (leap year)") {
        tm.setModelBaseDate(20240101);
        
        REQUIRE(tm.getModelBaseDate() == 20240101);
        REQUIRE(tm.getJulianDay() == 1);
    }
    
    SECTION("Set base date 20231231") {
        tm.setModelBaseDate(20231231);
        
        REQUIRE(tm.getModelBaseDate() == 20231231);
        REQUIRE(tm.getJulianDay() == 365);
    }
    
    SECTION("Set base date 20240229 (leap year Feb 29)") {
        tm.setModelBaseDate(20240229);
        
        REQUIRE(tm.getModelBaseDate() == 20240229);
        // Julian day for Feb 29 in leap year should be 60 (31+29)
        REQUIRE(tm.getJulianDay() == 60);
    }
}

TEST_CASE("TimeManager::updateTime updates all time variables", "[TimeManager]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("Update to 1440 minutes (1 day)") {
        tm.updateTime(1440.0);
        
        REQUIRE(tm.getT0() == 0.0);
        REQUIRE(tm.getT1() == 1440.0);
        REQUIRE(tm.getDT() == 1440.0);
        
        // Julian day should be 167 (166 + 1)
        REQUIRE(tm.getJulianDay() == 167);
    }
    
    SECTION("Update twice to test t0 = previous t1") {
        tm.updateTime(1440.0);  // First update
        REQUIRE(tm.getT0() == 0.0);
        REQUIRE(tm.getT1() == 1440.0);
        
        tm.updateTime(2880.0);  // Second update
        REQUIRE(tm.getT0() == 1440.0);
        REQUIRE(tm.getT1() == 2880.0);
        REQUIRE(tm.getDT() == 1440.0);
    }
    
    SECTION("Update with fractional minutes") {
        tm.updateTime(1500.5);
        
        REQUIRE(almost_equal(tm.getT1(), 1500.5));
        REQUIRE(almost_equal(tm.getDT(), 1500.5));
    }
    
    SECTION("Multiple updates maintain consistency") {
        double times[] = {60.0, 120.0, 180.0, 240.0};
        double prev_t1 = 0.0;
        
        for (double t : times) {
            tm.updateTime(t);
            REQUIRE(tm.getT0() == prev_t1);
            REQUIRE(tm.getT1() == t);
            REQUIRE(tm.getDT() == (t - prev_t1));
            prev_t1 = t;
        }
    }
}

TEST_CASE("TimeManager::formatTime supports 7 formats", "[TimeManager]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("Format type 1: YYYY") {
        std::string result = tm.formatBaseTime(1);
        REQUIRE(result == "2023");
    }
    
    SECTION("Format type 2: YYYY-MM") {
        std::string result = tm.formatBaseTime(2);
        REQUIRE(result == "2023-06");
    }
    
    SECTION("Format type 3: YYYY-MM-DD") {
        std::string result = tm.formatBaseTime(3);
        REQUIRE(result == "2023-06-15");
    }
    
    SECTION("Format type 4: YYYY-MM-DD hh") {
        std::string result = tm.formatBaseTime(4);
        REQUIRE(result == "2023-06-15 00");
    }
    
    SECTION("Format type 5: YYYY-MM-DD hh:mm") {
        std::string result = tm.formatBaseTime(5);
        REQUIRE(result == "2023-06-15 00:00");
    }
    
    SECTION("Format type 6: YYYY-MM-DD hh:mm:ss (default)") {
        std::string result = tm.formatBaseTime(6);
        REQUIRE(result == "2023-06-15 00:00:00");
        
        // Test default parameter
        std::string result_default = tm.formatBaseTime();
        REQUIRE(result_default == "2023-06-15 00:00:00");
    }
    
    SECTION("Format type 7: YYYY-MM-DD hh:mm:ss.ssssss") {
        std::string result = tm.formatBaseTime(7);
        // Should start with the date and time
        REQUIRE(result.substr(0, 19) == "2023-06-15 00:00:00");
        // Should have microseconds (6 digits after the dot)
        REQUIRE(result.length() == 26);
        REQUIRE(result[19] == '.');
    }
}

TEST_CASE("TimeManager::formatLocalTime after time update", "[TimeManager]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("After 1 day (1440 minutes)") {
        tm.updateTime(1440.0);
        std::string result = tm.formatLocalTime(3);
        REQUIRE(result == "2023-06-16");
    }
    
    SECTION("After 1 hour (60 minutes)") {
        tm.updateTime(60.0);
        std::string result = tm.formatLocalTime(6);
        REQUIRE(result == "2023-06-15 01:00:00");
    }
    
    SECTION("After 90 minutes (1.5 hours)") {
        tm.updateTime(90.0);
        std::string result = tm.formatLocalTime(5);
        REQUIRE(result == "2023-06-15 01:30");
    }
    
    SECTION("After multiple days") {
        tm.updateTime(7200.0);  // 5 days
        std::string result = tm.formatLocalTime(3);
        REQUIRE(result == "2023-06-20");
    }
}

TEST_CASE("TimeManager::Julian day calculation across year boundary", "[TimeManager]") {
    TimeManager tm;
    
    SECTION("Dec 31 to Jan 1") {
        tm.setModelBaseDate(20231231);
        REQUIRE(tm.getJulianDay() == 365);
        
        tm.updateTime(1440.0);  // Add 1 day
        REQUIRE(tm.getJulianDay() == 1);  // Should wrap to day 1 of next year
    }
    
    SECTION("Leap year Dec 31 to Jan 1") {
        tm.setModelBaseDate(20241231);
        REQUIRE(tm.getJulianDay() == 366);
        
        tm.updateTime(1440.0);  // Add 1 day
        REQUIRE(tm.getJulianDay() == 1);
    }
}

TEST_CASE("TimeManager::Julian day calculation for leap years", "[TimeManager]") {
    TimeManager tm;
    
    SECTION("2024 is a leap year") {
        tm.setModelBaseDate(20240301);  // March 1
        // Jan(31) + Feb(29) + 1 = 61
        REQUIRE(tm.getJulianDay() == 61);
    }
    
    SECTION("2023 is not a leap year") {
        tm.setModelBaseDate(20230301);  // March 1
        // Jan(31) + Feb(28) + 1 = 60
        REQUIRE(tm.getJulianDay() == 60);
    }
    
    SECTION("2000 is a leap year (divisible by 400)") {
        tm.setModelBaseDate(20000301);
        REQUIRE(tm.getJulianDay() == 61);
    }
    
    SECTION("1900 is not a leap year (divisible by 100 but not 400)") {
        tm.setModelBaseDate(19000301);
        REQUIRE(tm.getJulianDay() == 60);
    }
}

TEST_CASE("TimeManager::Time consistency checks", "[TimeManager]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("baseTime and modelBaseDate are consistent") {
        auto base_time = tm.getBaseTime();
        std::time_t tt = std::chrono::system_clock::to_time_t(base_time);
        std::tm* local_tm = std::localtime(&tt);
        
        int year = local_tm->tm_year + 1900;
        int month = local_tm->tm_mon + 1;
        int day = local_tm->tm_mday;
        
        REQUIRE(year == 2023);
        REQUIRE(month == 6);
        REQUIRE(day == 15);
    }
    
    SECTION("timelocal updates correctly") {
        tm.updateTime(1440.0);  // 1 day
        
        auto local_time = tm.getLocalTime();
        std::time_t tt = std::chrono::system_clock::to_time_t(local_time);
        std::tm* local_tm = std::localtime(&tt);
        
        int day = local_tm->tm_mday;
        REQUIRE(day == 16);  // Should be June 16
    }
}

TEST_CASE("TimeManager::Edge cases", "[TimeManager]") {
    TimeManager tm;
    
    SECTION("Very large elapsed time") {
        tm.setModelBaseDate(20230101);
        tm.updateTime(525600.0);  // 365 days in minutes
        
        std::string result = tm.formatLocalTime(3);
        REQUIRE(result == "2024-01-01");
    }
    
    SECTION("Zero elapsed time") {
        tm.setModelBaseDate(20230615);
        tm.updateTime(0.0);
        
        REQUIRE(tm.getT1() == 0.0);
        REQUIRE(tm.getDT() == 0.0);
        std::string result = tm.formatLocalTime(3);
        REQUIRE(result == "2023-06-15");
    }
    
    SECTION("Small fractional time") {
        tm.setModelBaseDate(20230615);
        tm.updateTime(0.0166667);  // ~1 second
        
        REQUIRE(almost_equal(tm.getT1(), 0.0166667, 1e-5));
    }
}

// ============================================================================
// Tests for Task 5.1: Time Variable Synchronization
// Requirements 2.2, 2.3, 2.4, 2.5, 2.6
// ============================================================================

TEST_CASE("TimeManager::Single time step synchronization", "[TimeManager][Task5.1]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("All time variables update synchronously in single step") {
        // Initial state
        REQUIRE(tm.getT0() == 0.0);
        REQUIRE(tm.getT1() == 0.0);
        REQUIRE(tm.getDT() == 0.0);
        REQUIRE(tm.getJulianDay() == 166);
        
        // Update to 60 minutes (1 hour)
        tm.updateTime(60.0);
        
        // Verify all variables updated synchronously
        REQUIRE(tm.getT0() == 0.0);           // t0 = previous t1
        REQUIRE(tm.getT1() == 60.0);          // t1 = current time
        REQUIRE(tm.getDT() == 60.0);          // dt = t1 - t0
        REQUIRE(tm.getJulianDay() == 166);    // Still same day
        
        // Verify timelocal and timeutc are updated
        auto local_time = tm.getLocalTime();
        auto utc_time = tm.getUTCTime();
        REQUIRE(local_time == tm.getBaseTime() + std::chrono::minutes(60));
        REQUIRE(utc_time == tm.getBaseTime() + std::chrono::minutes(60));
    }
    
    SECTION("Time variables update correctly for 1 day step") {
        tm.updateTime(1440.0);  // 1 day = 1440 minutes
        
        REQUIRE(tm.getT0() == 0.0);
        REQUIRE(tm.getT1() == 1440.0);
        REQUIRE(tm.getDT() == 1440.0);
        REQUIRE(tm.getJulianDay() == 167);  // Next day
        
        // Verify formatted time
        std::string formatted = tm.formatLocalTime(3);
        REQUIRE(formatted == "2023-06-16");
    }
}

TEST_CASE("TimeManager::Multiple consecutive time steps", "[TimeManager][Task5.1]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("Three consecutive 30-minute steps") {
        // Step 1: 0 -> 30 minutes
        tm.updateTime(30.0);
        REQUIRE(tm.getT0() == 0.0);
        REQUIRE(tm.getT1() == 30.0);
        REQUIRE(tm.getDT() == 30.0);
        
        // Step 2: 30 -> 60 minutes
        tm.updateTime(60.0);
        REQUIRE(tm.getT0() == 30.0);   // t0 = previous t1
        REQUIRE(tm.getT1() == 60.0);
        REQUIRE(tm.getDT() == 30.0);   // dt = 60 - 30
        
        // Step 3: 60 -> 90 minutes
        tm.updateTime(90.0);
        REQUIRE(tm.getT0() == 60.0);   // t0 = previous t1
        REQUIRE(tm.getT1() == 90.0);
        REQUIRE(tm.getDT() == 30.0);   // dt = 90 - 60
        
        // Verify time formatting
        std::string formatted = tm.formatLocalTime(5);
        REQUIRE(formatted == "2023-06-15 01:30");
    }
    
    SECTION("Variable time steps") {
        double steps[] = {15.0, 45.0, 120.0, 240.0, 480.0};
        double expected_dt[] = {15.0, 30.0, 75.0, 120.0, 240.0};
        
        for (size_t i = 0; i < 5; i++) {
            double prev_t1 = tm.getT1();
            tm.updateTime(steps[i]);
            
            REQUIRE(tm.getT0() == prev_t1);
            REQUIRE(tm.getT1() == steps[i]);
            REQUIRE(almost_equal(tm.getDT(), expected_dt[i]));
        }
    }
}

TEST_CASE("TimeManager::Time variables across day boundaries", "[TimeManager][Task5.1]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("Step across midnight") {
        // Start at 23:00 (1380 minutes from midnight)
        tm.updateTime(1380.0);
        REQUIRE(tm.getJulianDay() == 166);
        std::string time1 = tm.formatLocalTime(5);
        REQUIRE(time1 == "2023-06-15 23:00");
        
        // Step to 01:00 next day (1500 minutes from original midnight)
        tm.updateTime(1500.0);
        REQUIRE(tm.getT0() == 1380.0);
        REQUIRE(tm.getT1() == 1500.0);
        REQUIRE(tm.getDT() == 120.0);  // 2 hours
        REQUIRE(tm.getJulianDay() == 167);  // Next day
        
        std::string time2 = tm.formatLocalTime(5);
        REQUIRE(time2 == "2023-06-16 01:00");
    }
    
    SECTION("Multiple days progression") {
        int jd_start = tm.getJulianDay();
        
        for (int day = 1; day <= 5; day++) {
            tm.updateTime(day * 1440.0);
            REQUIRE(tm.getJulianDay() == jd_start + day);
        }
    }
}

TEST_CASE("TimeManager::Julian day updates correctly", "[TimeManager][Task5.1]") {
    TimeManager tm;
    
    SECTION("Julian day increments with elapsed days") {
        tm.setModelBaseDate(20230615);
        int initial_jd = tm.getJulianDay();
        REQUIRE(initial_jd == 166);
        
        // Add 1 day
        tm.updateTime(1440.0);
        REQUIRE(tm.getJulianDay() == 167);
        
        // Add another day
        tm.updateTime(2880.0);
        REQUIRE(tm.getJulianDay() == 168);
        
        // Add 10 more days
        tm.updateTime(17280.0);  // 12 days total
        REQUIRE(tm.getJulianDay() == 178);
    }
    
    SECTION("Julian day wraps at year boundary") {
        tm.setModelBaseDate(20231230);
        REQUIRE(tm.getJulianDay() == 364);
        
        tm.updateTime(1440.0);  // Dec 31
        REQUIRE(tm.getJulianDay() == 365);
        
        tm.updateTime(2880.0);  // Jan 1, next year
        REQUIRE(tm.getJulianDay() == 1);
    }
}

TEST_CASE("TimeManager::timelocal and timeutc synchronization", "[TimeManager][Task5.1]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("timelocal updates with elapsed time") {
        auto base_time = tm.getBaseTime();
        
        tm.updateTime(120.0);  // 2 hours
        auto local_time = tm.getLocalTime();
        
        auto expected_time = base_time + std::chrono::minutes(120);
        auto diff = std::chrono::duration_cast<std::chrono::seconds>(
            local_time - expected_time
        ).count();
        
        REQUIRE(std::abs(diff) < 1);  // Within 1 second tolerance
    }
    
    SECTION("timelocal and timeutc are synchronized") {
        tm.updateTime(360.0);  // 6 hours
        
        auto local_time = tm.getLocalTime();
        auto utc_time = tm.getUTCTime();
        
        // For this test, we expect them to be the same
        // (actual timezone handling may differ in implementation)
        auto diff = std::chrono::duration_cast<std::chrono::seconds>(
            local_time - utc_time
        ).count();
        
        REQUIRE(std::abs(diff) < 1);
    }
}

TEST_CASE("TimeManager::All time variables remain consistent", "[TimeManager][Task5.1]") {
    TimeManager tm;
    tm.setModelBaseDate(20230101);
    
    SECTION("Consistency over long simulation") {
        // Simulate 30 days with varying time steps
        double current_time = 0.0;
        double time_steps[] = {30.0, 60.0, 45.0, 90.0, 120.0};
        
        for (int day = 0; day < 30; day++) {
            for (double step : time_steps) {
                double prev_t1 = tm.getT1();
                current_time += step;
                tm.updateTime(current_time);
                
                // Verify t0 = previous t1 (Requirement 2.3)
                REQUIRE(tm.getT0() == prev_t1);
                
                // Verify t1 = current_time
                REQUIRE(tm.getT1() == current_time);
                
                // Verify dt = t1 - t0 (Requirement 2.4)
                REQUIRE(almost_equal(tm.getDT(), current_time - prev_t1));
                
                // Verify Julian day is in valid range (Requirement 2.5)
                int jd = tm.getJulianDay();
                REQUIRE(jd >= 1);
                REQUIRE(jd <= 366);
                
                // Verify timelocal is consistent (Requirement 2.6)
                auto local_time = tm.getLocalTime();
                auto expected = tm.getBaseTime() + 
                    std::chrono::microseconds(static_cast<long long>(current_time * 60 * 1000000));
                auto diff = std::chrono::duration_cast<std::chrono::seconds>(
                    local_time - expected
                ).count();
                REQUIRE(std::abs(diff) < 1);
            }
        }
    }
}

TEST_CASE("TimeManager::Fractional minute precision", "[TimeManager][Task5.1]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("Sub-minute time steps maintain precision") {
        // Test with 0.5 minute (30 second) steps
        tm.updateTime(0.5);
        REQUIRE(almost_equal(tm.getT1(), 0.5));
        REQUIRE(almost_equal(tm.getDT(), 0.5));
        
        tm.updateTime(1.0);
        REQUIRE(almost_equal(tm.getT0(), 0.5));
        REQUIRE(almost_equal(tm.getT1(), 1.0));
        REQUIRE(almost_equal(tm.getDT(), 0.5));
        
        tm.updateTime(1.25);
        REQUIRE(almost_equal(tm.getT0(), 1.0));
        REQUIRE(almost_equal(tm.getT1(), 1.25));
        REQUIRE(almost_equal(tm.getDT(), 0.25));
    }
}
