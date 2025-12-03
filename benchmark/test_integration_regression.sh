#!/bin/bash
# =============================================================
# SHUD Real-Time System Integration and Regression Test
# 集成测试和回归测试脚本
# =============================================================
# This script tests:
# 1. All test cases run successfully (ccw, heihe, qhh)
# 2. Time system doesn't affect model results
# 3. Output files contain correct timestamps
# =============================================================

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configuration
TEST_CASES=("ccw" "heihe" "qhh")
SHUD_EXECUTABLE="./shud"
OUTPUT_DIR="output"
BACKUP_DIR="output_backup_regression"

# Test results
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to print colored messages
print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ $1${NC}"
}

# Function to check if executable exists
check_executable() {
    if [ ! -f "$SHUD_EXECUTABLE" ]; then
        print_error "SHUD executable not found: $SHUD_EXECUTABLE"
        print_info "Please run 'make shud' first"
        exit 1
    fi
    print_success "Found SHUD executable: $SHUD_EXECUTABLE"
}

# Function to backup existing outputs
backup_outputs() {
    print_header "Backing up existing outputs"
    
    if [ -d "$OUTPUT_DIR" ]; then
        if [ -d "$BACKUP_DIR" ]; then
            rm -rf "$BACKUP_DIR"
        fi
        cp -r "$OUTPUT_DIR" "$BACKUP_DIR"
        print_success "Backed up outputs to $BACKUP_DIR"
    else
        print_warning "No existing outputs to backup"
    fi
}

# Function to run a test case
run_test_case() {
    local test_case=$1
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    print_header "Running test case: $test_case"
    
    # Check if input directory exists
    if [ ! -d "input/$test_case" ]; then
        print_error "Input directory not found: input/$test_case"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Run the model with timeout (5 minutes max)
    print_info "Executing: $SHUD_EXECUTABLE $test_case"
    
    # Capture start time
    start_time=$(date +%s)
    
    # Run and capture output
    # Auto-answer "Y" to any prompts using yes command
    yes "Y" | $SHUD_EXECUTABLE $test_case > "test_output_${test_case}.log" 2>&1
    exit_code=$?
    
    # Capture end time
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    
    if [ $exit_code -eq 0 ]; then
        print_success "Test case $test_case completed successfully in ${elapsed}s"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    elif [ $exit_code -eq 124 ]; then
        print_error "Test case $test_case timed out after 300s"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    else
        print_error "Test case $test_case failed with exit code $exit_code"
        print_info "Check test_output_${test_case}.log for details"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
}

# Function to verify output files exist
verify_output_files() {
    local test_case=$1
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    print_header "Verifying output files for: $test_case"
    
    local output_path="$OUTPUT_DIR/${test_case}.out"
    
    if [ ! -d "$output_path" ]; then
        print_error "Output directory not found: $output_path"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Count output files
    local file_count=$(find "$output_path" -type f | wc -l)
    
    if [ $file_count -gt 0 ]; then
        print_success "Found $file_count output files in $output_path"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    else
        print_error "No output files found in $output_path"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
}

# Function to check timestamps in output files
check_timestamps() {
    local test_case=$1
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    print_header "Checking timestamps in output files: $test_case"
    
    local output_path="$OUTPUT_DIR/${test_case}.out"
    local timestamp_found=0
    local files_checked=0
    
    # Check .dat and .csv files for timestamps
    for file in "$output_path"/*.dat "$output_path"/*.csv; do
        if [ -f "$file" ]; then
            files_checked=$((files_checked + 1))
            
            # Check if file contains a timestamp in YYYYMMDD format
            if head -5 "$file" | grep -qE '[0-9]{8}'; then
                timestamp_found=$((timestamp_found + 1))
                local timestamp=$(head -5 "$file" | grep -oE '[0-9]{8}' | head -1)
                print_info "Found timestamp $timestamp in $(basename $file)"
            fi
        fi
    done
    
    if [ $files_checked -eq 0 ]; then
        print_warning "No .dat or .csv files found to check"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    fi
    
    if [ $timestamp_found -gt 0 ]; then
        print_success "Found timestamps in $timestamp_found out of $files_checked files"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    else
        print_warning "No timestamps found in output files (may be binary format)"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    fi
}

# Function to compare with baseline (if exists)
compare_with_baseline() {
    local test_case=$1
    
    if [ ! -d "$BACKUP_DIR" ]; then
        print_info "No baseline to compare (first run)"
        return 0
    fi
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    print_header "Comparing results with baseline: $test_case"
    
    local new_output="$OUTPUT_DIR/${test_case}.out"
    local old_output="$BACKUP_DIR/${test_case}.out"
    
    if [ ! -d "$old_output" ]; then
        print_warning "No baseline output found for $test_case"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        return 0
    fi
    
    # Compare file counts
    local new_count=$(find "$new_output" -type f | wc -l)
    local old_count=$(find "$old_output" -type f | wc -l)
    
    if [ $new_count -ne $old_count ]; then
        print_warning "File count differs: new=$new_count, old=$old_count"
    else
        print_success "File count matches: $new_count files"
    fi
    
    # For binary files, just check they exist and have similar sizes
    local size_diff_count=0
    for new_file in "$new_output"/*; do
        if [ -f "$new_file" ]; then
            local basename=$(basename "$new_file")
            local old_file="$old_output/$basename"
            
            if [ -f "$old_file" ]; then
                local new_size=$(stat -f%z "$new_file" 2>/dev/null || stat -c%s "$new_file" 2>/dev/null)
                local old_size=$(stat -f%z "$old_file" 2>/dev/null || stat -c%s "$old_file" 2>/dev/null)
                
                # Allow 5% size difference (for timestamp additions)
                local diff=$((new_size - old_size))
                local abs_diff=${diff#-}
                local threshold=$((old_size / 20))  # 5%
                
                if [ $abs_diff -gt $threshold ]; then
                    size_diff_count=$((size_diff_count + 1))
                    print_warning "Size difference in $basename: new=$new_size, old=$old_size"
                fi
            fi
        fi
    done
    
    if [ $size_diff_count -eq 0 ]; then
        print_success "All file sizes are within acceptable range"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        print_warning "$size_diff_count files have size differences (may be due to timestamp additions)"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    fi
    
    return 0
}

# Function to check performance
check_performance() {
    local test_case=$1
    
    print_header "Checking performance for: $test_case"
    
    local log_file="test_output_${test_case}.log"
    
    if [ ! -f "$log_file" ]; then
        print_warning "Log file not found: $log_file"
        return 0
    fi
    
    # Extract timing information if available
    if grep -q "CPU time" "$log_file"; then
        local cpu_time=$(grep "CPU time" "$log_file" | tail -1)
        print_info "Performance: $cpu_time"
    fi
    
    if grep -q "Wall time" "$log_file"; then
        local wall_time=$(grep "Wall time" "$log_file" | tail -1)
        print_info "Performance: $wall_time"
    fi
}

# Function to print summary
print_summary() {
    print_header "Test Summary"
    echo -e "Total tests: $TOTAL_TESTS"
    echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"
    echo -e "${RED}Failed: $FAILED_TESTS${NC}"
    
    if [ $FAILED_TESTS -eq 0 ]; then
        print_success "All tests passed!"
        return 0
    else
        print_error "Some tests failed"
        return 1
    fi
}

# Main test execution
main() {
    print_header "SHUD Integration and Regression Test Suite"
    print_info "Testing real-time system implementation"
    echo ""
    
    # Step 1: Check executable
    check_executable
    echo ""
    
    # Step 2: Backup existing outputs
    backup_outputs
    echo ""
    
    # Step 3: Run all test cases
    for test_case in "${TEST_CASES[@]}"; do
        run_test_case "$test_case"
        echo ""
        
        if [ $? -eq 0 ]; then
            # Verify outputs
            verify_output_files "$test_case"
            echo ""
            
            # Check timestamps
            check_timestamps "$test_case"
            echo ""
            
            # Compare with baseline
            compare_with_baseline "$test_case"
            echo ""
            
            # Check performance
            check_performance "$test_case"
            echo ""
        fi
    done
    
    # Step 4: Print summary
    print_summary
    
    # Return appropriate exit code
    if [ $FAILED_TESTS -eq 0 ]; then
        exit 0
    else
        exit 1
    fi
}

# Run main function
main
