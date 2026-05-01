//
// Created by clz on 2024/10/23.
//

#include "gtest/gtest.h"

class TestCase : public ::testing::Test {

protected:
    void SetUp() override {
        filePath = TestDataPath;
    }
    void TearDown() override {
    }
public:
    std::string filePath ;
};

TEST_F(TestCase, CCW) {
    std::cout<<"ccw case test"<<std::endl;
}

