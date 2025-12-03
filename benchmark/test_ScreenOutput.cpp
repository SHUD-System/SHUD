#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../src/classes/TimeManager.hpp"
#include <sstream>
#include <iostream>

// Test helper to capture stdout
class StdoutCapture {
private:
    std::streambuf* old_stdout;
    std::ostringstream captured;
    
public:
    StdoutCapture() {
        old_stdout = std::cout.rdbuf(captured.rdbuf());
    }
    
    ~StdoutCapture() {
        std::cout.rdbuf(old_stdout);
    }
    
    std::string getOutput() {
        return captured.str();
    }
};

TEST_CASE("TimeManager formatLocalTime produces correct format", "[screen_output]") {
    TimeManager tm;
    
    // Set a known base date: 2023-06-15
    tm.setModelBaseDate(20230615);
    
    SECTION("Format at base time (t=0)") {
        tm.updateTime(0.0);
        std::string formatted = tm.formatLocalTime();
        
        // Should be 2023-06-15 00:00:00
        REQUIRE(formatted.find("2023-06-15") != std::string::npos);
        REQUIRE(formatted.find("00:00:00") != std::string::npos);
    }
    
    SECTION("Format after 1 day (t=1440 minutes)") {
        tm.updateTime(1440.0);
        std::string formatted = tm.formatLocalTime();
        
        // Should be 2023-06-16 00:00:00
        REQUIRE(formatted.find("2023-06-16") != std::string::npos);
        REQUIRE(formatted.find("00:00:00") != std::string::npos);
    }
    
    SECTION("Format after 12 hours (t=720 minutes)") {
        tm.updateTime(720.0);
        std::string formatted = tm.formatLocalTime();
        
        // Should be 2023-06-15 12:00:00
        REQUIRE(formatted.find("2023-06-15") != std::string::npos);
        REQUIRE(formatted.find("12:00:00") != std::string::npos);
    }
    
    SECTION("Format with different precision levels") {
        tm.updateTime(0.0);
        
        // Type 3: YYYY-MM-DD
        std::string fmt3 = tm.formatLocalTime(3);
        REQUIRE(fmt3 == "2023-06-15");
        
        // Type 6: YYYY-MM-DD hh:mm:ss (default)
        std::string fmt6 = tm.formatLocalTime(6);
        REQUIRE(fmt6.find("2023-06-15") != std::string::npos);
        REQUIRE(fmt6.find("00:00:00") != std::string::npos);
    }
    
    SECTION("Format without separators (compact format)") {
        tm.updateTime(0.0);
        
        // Type 3 without separators: YYYYMMDD
        std::string fmt3_compact = tm.formatLocalTime(3, false);
        REQUIRE(fmt3_compact == "20230615");
        
        // Type 6 without separators: YYYYMMDDhhmmss
        std::string fmt6_compact = tm.formatLocalTime(6, false);
        REQUIRE(fmt6_compact == "20230615000000");
        
        // Type 5 without separators: YYYYMMDDhhmm
        std::string fmt5_compact = tm.formatLocalTime(5, false);
        REQUIRE(fmt5_compact == "202306150000");
    }
    
    SECTION("Format with separators (default behavior)") {
        tm.updateTime(720.0); // 12 hours
        
        // With separators (default)
        std::string with_sep = tm.formatLocalTime(6, true);
        REQUIRE(with_sep == "2023-06-15 12:00:00");
        
        // Without separators
        std::string without_sep = tm.formatLocalTime(6, false);
        REQUIRE(without_sep == "20230615120000");
    }
}

TEST_CASE("Screen output format verification", "[screen_output]") {
    SECTION("ScreenPrint output includes real time") {
        // This test verifies that the ScreenPrint function would include
        // real time in its output. Since we can't easily test the actual
        // function without a full Model_Data setup, we verify the format
        // that would be produced.
        
        TimeManager tm;
        tm.setModelBaseDate(20230615);
        tm.updateTime(1440.0); // 1 day
        
        std::string realTime = tm.formatLocalTime();
        
        // Verify the format matches what ScreenPrint expects
        REQUIRE(realTime.length() > 0);
        REQUIRE(realTime.find("2023-06-16") != std::string::npos);
        
        // Simulate the printf format used in ScreenPrint
        char buffer[256];
        snprintf(buffer, sizeof(buffer), 
                "%.2f day \t %s \t %.2f%% \t %.2f s \t %.2f s \t %ld \n",
                1440.0 / 1440.0, realTime.c_str(), 50.0, 1.23, 2.45, 1234L);
        
        std::string output(buffer);
        REQUIRE(output.find("1.00 day") != std::string::npos);
        REQUIRE(output.find("2023-06-16") != std::string::npos);
        REQUIRE(output.find("50.00%") != std::string::npos);
    }
    
    SECTION("ScreenPrintu output includes real time") {
        TimeManager tm;
        tm.setModelBaseDate(20230615);
        tm.updateTime(720.0); // 0.5 day
        
        std::string realTime = tm.formatLocalTime();
        
        // Simulate the printf format used in ScreenPrintu
        char buffer[256];
        snprintf(buffer, sizeof(buffer),
                "%6.2f d \t %s \t %5.2f%% \t %6.2f s \t %6ld %6ld %6ld %6ld %6ld\n",
                720.0 / 1440.0, realTime.c_str(), 25.0, 1.5, 10L, 20L, 30L, 40L, 50L);
        
        std::string output(buffer);
        REQUIRE(output.find("0.50 d") != std::string::npos);
        REQUIRE(output.find("2023-06-15") != std::string::npos);
        REQUIRE(output.find("12:00:00") != std::string::npos);
    }
    
    SECTION("modelSummary output format") {
        TimeManager tm;
        tm.setModelBaseDate(20230615);
        
        // Test start time
        tm.updateTime(0.0);
        std::string startTime = tm.formatLocalTime();
        
        // Test end time
        tm.updateTime(1440.0 * 365); // 365 days
        std::string endTime = tm.formatLocalTime();
        
        // Verify base date display
        long baseDate = tm.getModelBaseDate();
        REQUIRE(baseDate == 20230615);
        
        // Verify start time format
        REQUIRE(startTime.find("2023-06-15") != std::string::npos);
        REQUIRE(startTime.find("00:00:00") != std::string::npos);
        
        // Verify end time format (should be approximately 2024-06-14)
        REQUIRE(endTime.find("2024-06") != std::string::npos);
    }
}

TEST_CASE("Real time display integration", "[screen_output]") {
    SECTION("Time progression through simulation") {
        TimeManager tm;
        tm.setModelBaseDate(20000101); // Y2K
        
        // Simulate time progression
        std::vector<double> times = {0.0, 1440.0, 2880.0, 4320.0}; // 0, 1, 2, 3 days
        std::vector<std::string> expectedDates = {
            "2000-01-01",
            "2000-01-02", 
            "2000-01-03",
            "2000-01-04"
        };
        
        for (size_t i = 0; i < times.size(); i++) {
            tm.updateTime(times[i]);
            std::string formatted = tm.formatLocalTime();
            REQUIRE(formatted.find(expectedDates[i]) != std::string::npos);
        }
    }
}
