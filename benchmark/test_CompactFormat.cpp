#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../src/classes/TimeManager.hpp"

TEST_CASE("Compact format (no separators) functionality", "[compact_format]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("Type 1: YYYY") {
        tm.updateTime(0.0);
        REQUIRE(tm.formatLocalTime(1, true) == "2023");
        REQUIRE(tm.formatLocalTime(1, false) == "2023");  // No difference for type 1
    }
    
    SECTION("Type 2: YYYY-MM vs YYYYMM") {
        tm.updateTime(0.0);
        REQUIRE(tm.formatLocalTime(2, true) == "2023-06");
        REQUIRE(tm.formatLocalTime(2, false) == "202306");
    }
    
    SECTION("Type 3: YYYY-MM-DD vs YYYYMMDD") {
        tm.updateTime(0.0);
        REQUIRE(tm.formatLocalTime(3, true) == "2023-06-15");
        REQUIRE(tm.formatLocalTime(3, false) == "20230615");
    }
    
    SECTION("Type 4: YYYY-MM-DD hh vs YYYYMMDDhh") {
        tm.updateTime(720.0);  // 12 hours
        REQUIRE(tm.formatLocalTime(4, true) == "2023-06-15 12");
        REQUIRE(tm.formatLocalTime(4, false) == "2023061512");
    }
    
    SECTION("Type 5: YYYY-MM-DD hh:mm vs YYYYMMDDhhmm") {
        tm.updateTime(750.0);  // 12 hours 30 minutes
        REQUIRE(tm.formatLocalTime(5, true) == "2023-06-15 12:30");
        REQUIRE(tm.formatLocalTime(5, false) == "202306151230");
    }
    
    SECTION("Type 6: YYYY-MM-DD hh:mm:ss vs YYYYMMDDhhmmss") {
        tm.updateTime(765.0);  // 12 hours 45 minutes
        REQUIRE(tm.formatLocalTime(6, true) == "2023-06-15 12:45:00");
        REQUIRE(tm.formatLocalTime(6, false) == "20230615124500");
    }
    
    SECTION("Type 7: YYYY-MM-DD hh:mm:ss.ssssss vs YYYYMMDDhhmmssssssss") {
        tm.updateTime(0.0);
        std::string with_sep = tm.formatLocalTime(7, true);
        std::string without_sep = tm.formatLocalTime(7, false);
        
        // With separators should contain "-", ":", " ", "."
        REQUIRE(with_sep.find("-") != std::string::npos);
        REQUIRE(with_sep.find(":") != std::string::npos);
        REQUIRE(with_sep.find(" ") != std::string::npos);
        REQUIRE(with_sep.find(".") != std::string::npos);
        
        // Without separators should not contain any separators
        REQUIRE(without_sep.find("-") == std::string::npos);
        REQUIRE(without_sep.find(":") == std::string::npos);
        REQUIRE(without_sep.find(" ") == std::string::npos);
        REQUIRE(without_sep.find(".") == std::string::npos);
        
        // Should start with date
        REQUIRE(without_sep.substr(0, 8) == "20230615");
    }
}

TEST_CASE("Compact format across different dates", "[compact_format]") {
    TimeManager tm;
    
    SECTION("Year 2000") {
        tm.setModelBaseDate(20000101);
        tm.updateTime(0.0);
        REQUIRE(tm.formatLocalTime(3, false) == "20000101");
        REQUIRE(tm.formatLocalTime(6, false) == "20000101000000");
    }
    
    SECTION("Leap year date") {
        tm.setModelBaseDate(20200229);  // Feb 29, 2020 (leap year)
        tm.updateTime(0.0);
        REQUIRE(tm.formatLocalTime(3, false) == "20200229");
        REQUIRE(tm.formatLocalTime(6, false) == "20200229000000");
    }
    
    SECTION("End of year") {
        tm.setModelBaseDate(20231231);
        tm.updateTime(1439.0);  // 23:59
        REQUIRE(tm.formatLocalTime(6, false) == "20231231235900");
    }
}

TEST_CASE("Compact format for all format functions", "[compact_format]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    tm.updateTime(720.0);  // 12 hours
    
    SECTION("formatLocalTime") {
        REQUIRE(tm.formatLocalTime(6, false) == "20230615120000");
    }
    
    SECTION("formatUTCTime") {
        REQUIRE(tm.formatUTCTime(6, false) == "20230615120000");
    }
    
    SECTION("formatBaseTime") {
        REQUIRE(tm.formatBaseTime(6, false) == "20230615000000");
    }
}

TEST_CASE("Default behavior unchanged", "[compact_format]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    tm.updateTime(0.0);
    
    SECTION("Default parameter uses separators") {
        // Calling without second parameter should use separators (default)
        std::string default_format = tm.formatLocalTime(6);
        REQUIRE(default_format == "2023-06-15 00:00:00");
        
        // Explicitly passing true should give same result
        std::string explicit_true = tm.formatLocalTime(6, true);
        REQUIRE(explicit_true == default_format);
    }
}

TEST_CASE("Practical use cases", "[compact_format]") {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    
    SECTION("File naming with compact timestamp") {
        tm.updateTime(765.0);  // 12:45:00
        std::string timestamp = tm.formatLocalTime(6, false);
        std::string filename = "output_" + timestamp + ".dat";
        REQUIRE(filename == "output_20230615124500.dat");
    }
    
    SECTION("Database key with compact date") {
        tm.updateTime(0.0);
        std::string date_key = tm.formatLocalTime(3, false);
        REQUIRE(date_key == "20230615");
        REQUIRE(date_key.length() == 8);  // Exactly 8 characters
    }
    
    SECTION("Log entry with readable vs compact format") {
        tm.updateTime(1440.0);  // 1 day later
        
        // Readable format for display
        std::string readable = tm.formatLocalTime(6, true);
        REQUIRE(readable == "2023-06-16 00:00:00");
        
        // Compact format for storage/processing
        std::string compact = tm.formatLocalTime(6, false);
        REQUIRE(compact == "20230616000000");
        REQUIRE(compact.length() == 14);  // Exactly 14 characters
    }
}
