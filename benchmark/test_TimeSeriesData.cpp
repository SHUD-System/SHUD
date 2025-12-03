#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../src/classes/TimeSeriesData.hpp"
#include <fstream>
#include <cstdio>

// Helper function to create a temporary TSD file for testing
void createTestTSDFile(const char* filename, int nrows, int ncols, long startTime) {
    FILE* fp = fopen(filename, "w");
    if (fp) {
        fprintf(fp, "%d %d %ld\n", nrows, ncols, startTime);
        fprintf(fp, "Time_Day\tCol1\tCol2\n");
        for (int i = 0; i < nrows; i++) {
            fprintf(fp, "%d\t1.0\t2.0\n", i);
        }
        fclose(fp);
    }
}

TEST_CASE("TimeSeriesData::getStartTime returns correct timestamp", "[TimeSeriesData]") {
    _TimeSeriesData tsd;
    
    SECTION("Read StartTime from ccw LAI file") {
        tsd.fn = "input/ccw/ccw.tsd.lai";
        tsd.readDimensions();
        
        long startTime = tsd.getStartTime();
        REQUIRE(startTime == 20000101);
    }
    
    SECTION("Read StartTime from ccw MF file") {
        tsd.fn = "input/ccw/ccw.tsd.mf";
        tsd.readDimensions();
        
        long startTime = tsd.getStartTime();
        REQUIRE(startTime == 20000101);
    }
    
    SECTION("Read StartTime from heihe LAI file") {
        tsd.fn = "input/heihe/heihe.tsd.lai";
        tsd.readDimensions();
        
        long startTime = tsd.getStartTime();
        REQUIRE(startTime == 19890101);
    }
    
    SECTION("Read StartTime from heihe MF file") {
        tsd.fn = "input/heihe/heihe.tsd.mf";
        tsd.readDimensions();
        
        long startTime = tsd.getStartTime();
        REQUIRE(startTime == 19890101);
    }
}

TEST_CASE("TimeSeriesData::getStartTime with different timestamps", "[TimeSeriesData]") {
    _TimeSeriesData tsd;
    
    SECTION("Test with timestamp 20230615") {
        const char* testFile = "test_tsd_20230615.tmp";
        createTestTSDFile(testFile, 10, 3, 20230615);
        
        tsd.fn = testFile;
        tsd.readDimensions();
        
        long startTime = tsd.getStartTime();
        REQUIRE(startTime == 20230615);
        
        // Clean up
        remove(testFile);
    }
    
    SECTION("Test with timestamp 19900101") {
        const char* testFile = "test_tsd_19900101.tmp";
        createTestTSDFile(testFile, 5, 2, 19900101);
        
        tsd.fn = testFile;
        tsd.readDimensions();
        
        long startTime = tsd.getStartTime();
        REQUIRE(startTime == 19900101);
        
        // Clean up
        remove(testFile);
    }
    
    SECTION("Test with timestamp 20241231") {
        const char* testFile = "test_tsd_20241231.tmp";
        createTestTSDFile(testFile, 12, 4, 20241231);
        
        tsd.fn = testFile;
        tsd.readDimensions();
        
        long startTime = tsd.getStartTime();
        REQUIRE(startTime == 20241231);
        
        // Clean up
        remove(testFile);
    }
}

TEST_CASE("TimeSeriesData::getStartTime is const method", "[TimeSeriesData]") {
    _TimeSeriesData tsd;
    tsd.fn = "input/ccw/ccw.tsd.lai";
    tsd.readDimensions();
    
    SECTION("Can be called on const object") {
        const _TimeSeriesData& const_tsd = tsd;
        long startTime = const_tsd.getStartTime();
        REQUIRE(startTime == 20000101);
    }
}

TEST_CASE("TimeSeriesData::getStartTime consistency across multiple calls", "[TimeSeriesData]") {
    _TimeSeriesData tsd;
    tsd.fn = "input/ccw/ccw.tsd.lai";
    tsd.readDimensions();
    
    SECTION("Multiple calls return same value") {
        long startTime1 = tsd.getStartTime();
        long startTime2 = tsd.getStartTime();
        long startTime3 = tsd.getStartTime();
        
        REQUIRE(startTime1 == startTime2);
        REQUIRE(startTime2 == startTime3);
        REQUIRE(startTime1 == 20000101);
    }
}
