#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../src/classes/TimeManager.hpp"
#include <cstdio>

TEST_CASE("TimeManager::modelBaseDate protection", "[TimeManager][Protection]") {
    TimeManager tm;
    
    SECTION("First initialization succeeds") {
        REQUIRE_FALSE(tm.isInitialized());
        
        tm.setModelBaseDate(20230615);
        
        REQUIRE(tm.isInitialized());
        REQUIRE(tm.getModelBaseDate() == 20230615);
    }
    
    SECTION("Second initialization is blocked") {
        // First initialization
        tm.setModelBaseDate(20230615);
        REQUIRE(tm.getModelBaseDate() == 20230615);
        
        // Capture stderr to verify warning message
        FILE* original_stderr = stderr;
        FILE* temp_stderr = tmpfile();
        stderr = temp_stderr;
        
        // Attempt second initialization (should be blocked)
        tm.setModelBaseDate(20240101);
        
        // Restore stderr
        fflush(temp_stderr);
        stderr = original_stderr;
        fclose(temp_stderr);
        
        // Verify that modelBaseDate was NOT changed
        REQUIRE(tm.getModelBaseDate() == 20230615);
        REQUIRE_FALSE(tm.getModelBaseDate() == 20240101);
    }
    
    SECTION("Multiple attempts to modify are all blocked") {
        tm.setModelBaseDate(20230615);
        long original = tm.getModelBaseDate();
        
        // Try multiple times to change it
        tm.setModelBaseDate(20240101);
        tm.setModelBaseDate(20250101);
        tm.setModelBaseDate(19990101);
        
        // Should still be the original value
        REQUIRE(tm.getModelBaseDate() == original);
        REQUIRE(tm.getModelBaseDate() == 20230615);
    }
    
    SECTION("Protection persists across time updates") {
        tm.setModelBaseDate(20230615);
        
        // Update time
        tm.updateTime(1440.0);  // 1 day
        
        // Try to change base date after time update
        tm.setModelBaseDate(20240101);
        
        // Should still be original
        REQUIRE(tm.getModelBaseDate() == 20230615);
    }
}

TEST_CASE("TimeManager::isInitialized flag", "[TimeManager][Protection]") {
    TimeManager tm;
    
    SECTION("Initially not initialized") {
        REQUIRE_FALSE(tm.isInitialized());
    }
    
    SECTION("Becomes initialized after setModelBaseDate") {
        REQUIRE_FALSE(tm.isInitialized());
        
        tm.setModelBaseDate(20230615);
        
        REQUIRE(tm.isInitialized());
    }
    
    SECTION("Remains initialized after time updates") {
        tm.setModelBaseDate(20230615);
        REQUIRE(tm.isInitialized());
        
        tm.updateTime(1440.0);
        REQUIRE(tm.isInitialized());
        
        tm.updateTime(2880.0);
        REQUIRE(tm.isInitialized());
    }
}
