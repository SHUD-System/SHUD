#include "../src/classes/TimeManager.hpp"
#include <iostream>
#include <iomanip>

int main() {
    TimeManager tm;
    tm.setModelBaseDate(20230615);
    tm.updateTime(765.0);  // 12小时45分钟
    
    std::cout << "========================================" << std::endl;
    std::cout << "TimeManager 紧凑格式功能演示" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;
    
    std::cout << "基准日期: 20230615" << std::endl;
    std::cout << "经过时间: 765分钟 (12小时45分钟)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== 格式对照表 ===" << std::endl;
    std::cout << std::left << std::setw(8) << "Type" 
              << std::setw(30) << "带分隔符 (默认)" 
              << std::setw(30) << "紧凑格式 (无分隔符)" << std::endl;
    std::cout << std::string(68, '-') << std::endl;
    
    for (int type = 1; type <= 7; type++) {
        std::string with_sep = tm.formatLocalTime(type, true);
        std::string without_sep = tm.formatLocalTime(type, false);
        
        std::cout << std::left << std::setw(8) << type
                  << std::setw(30) << with_sep
                  << std::setw(30) << without_sep << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "=== 实际应用示例 ===" << std::endl;
    std::cout << std::endl;
    
    // 示例 1: 文件命名
    std::cout << "1. 文件命名:" << std::endl;
    std::string filename = "output_" + tm.formatLocalTime(6, false) + ".dat";
    std::cout << "   " << filename << std::endl;
    std::cout << std::endl;
    
    // 示例 2: 数据库键
    std::cout << "2. 数据库日期键:" << std::endl;
    std::string date_key = tm.formatLocalTime(3, false);
    std::cout << "   " << date_key << " (长度: " << date_key.length() << " 字符)" << std::endl;
    std::cout << std::endl;
    
    // 示例 3: 日志记录
    std::cout << "3. 日志记录:" << std::endl;
    std::cout << "   屏幕显示: [" << tm.formatLocalTime(6, true) << "] 模拟开始" << std::endl;
    std::cout << "   文件存储: " << tm.formatLocalTime(6, false) << ",START,0" << std::endl;
    std::cout << std::endl;
    
    // 示例 4: CSV 输出
    std::cout << "4. CSV 数据输出:" << std::endl;
    std::cout << "   紧凑格式: " << tm.formatLocalTime(6, false) << ",123.45,678.90" << std::endl;
    std::cout << "   可读格式: " << tm.formatLocalTime(6, true) << ",123.45,678.90" << std::endl;
    std::cout << std::endl;
    
    // 示例 5: 不同时间点
    std::cout << "5. 模拟时间进展:" << std::endl;
    std::cout << "   " << std::left << std::setw(20) << "相对时间" 
              << std::setw(25) << "可读格式" 
              << std::setw(20) << "紧凑格式" << std::endl;
    std::cout << "   " << std::string(65, '-') << std::endl;
    
    double times[] = {0.0, 360.0, 720.0, 1440.0, 2880.0};
    const char* labels[] = {"开始 (0分钟)", "6小时", "12小时", "1天", "2天"};
    
    for (int i = 0; i < 5; i++) {
        tm.updateTime(times[i]);
        std::cout << "   " << std::left << std::setw(20) << labels[i]
                  << std::setw(25) << tm.formatLocalTime(6, true)
                  << std::setw(20) << tm.formatLocalTime(6, false) << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "演示完成" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
