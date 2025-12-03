#!/bin/bash
# =============================================================
# SHUD Time Stamp Inconsistency Test
# 时间戳不一致场景测试
# =============================================================
# This script tests:
# 1. Detection of inconsistent timestamps
# 2. Warning message display
# 3. User confirmation mechanism (Y/n)
# =============================================================

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configuration
SHUD_EXECUTABLE="./shud"
TEST_DIR="input/test_timestamp"
BACKUP_DIR="input/ccw_backup"

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

# Function to setup test environment
setup_test_env() {
    print_header "Setting up test environment"
    
    # Backup original ccw directory
    if [ -d "input/ccw" ] && [ ! -d "$BACKUP_DIR" ]; then
        cp -r input/ccw "$BACKUP_DIR"
        print_success "Backed up original ccw directory"
    fi
    
    # Create test directory
    if [ -d "$TEST_DIR" ]; then
        rm -rf "$TEST_DIR"
    fi
    cp -r input/ccw "$TEST_DIR"
    print_success "Created test directory: $TEST_DIR"
}

# Function to cleanup test environment
cleanup_test_env() {
    print_header "Cleaning up test environment"
    
    if [ -d "$TEST_DIR" ]; then
        rm -rf "$TEST_DIR"
        print_success "Removed test directory"
    fi
    
    if [ -d "$BACKUP_DIR" ]; then
        rm -rf "$BACKUP_DIR"
        print_success "Removed backup directory"
    fi
}

# Function to modify timestamp in tsd.forc file
modify_tsd_forc_timestamp() {
    local new_timestamp=$1
    local file="$TEST_DIR/test_timestamp.tsd.forc"
    
    # Copy and modify the file
    cp "$TEST_DIR/ccw.tsd.forc" "$file"
    
    # Replace the timestamp (first line, second field)
    sed -i.bak "s/^1 [0-9]*/1 $new_timestamp/" "$file"
    rm -f "$file.bak"
    
    print_info "Modified tsd.forc timestamp to: $new_timestamp"
}

# Function to modify timestamp in forcing.csv file
modify_forcing_csv_timestamp() {
    local new_timestamp=$1
    local file="$TEST_DIR/forcing_modified.csv"
    
    # Copy and modify the file
    cp "$TEST_DIR/forcing.csv" "$file"
    
    # Replace the timestamp in the header (first line, third field)
    local first_line=$(head -1 "$file")
    local modified_line=$(echo "$first_line" | awk -v ts="$new_timestamp" '{$3=ts; print}')
    
    # Create temp file with modified header
    echo "$modified_line" > "$file.tmp"
    tail -n +2 "$file" >> "$file.tmp"
    mv "$file.tmp" "$file"
    
    print_info "Modified forcing.csv timestamp to: $new_timestamp"
}

# Function to create config file pointing to modified files
create_test_config() {
    local test_name=$1
    
    # Copy all config files
    for file in "$TEST_DIR"/ccw.*; do
        local basename=$(basename "$file")
        local newname=$(echo "$basename" | sed "s/ccw/$test_name/")
        cp "$file" "$TEST_DIR/$newname"
    done
    
    print_info "Created test config files for: $test_name"
}

# Test 1: Inconsistent timestamp - should show warning
test_inconsistent_timestamp_warning() {
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    print_header "Test 1: Inconsistent Timestamp Warning"
    
    setup_test_env
    
    # Modify tsd.forc to have timestamp 20000101
    modify_tsd_forc_timestamp "20000101"
    
    # Modify forcing.csv to have different timestamp 20050101
    modify_forcing_csv_timestamp "20050101"
    
    # Update tsd.forc to point to modified forcing file
    sed -i.bak "s/forcing.csv/forcing_modified.csv/" "$TEST_DIR/test_timestamp.tsd.forc"
    rm -f "$TEST_DIR/test_timestamp.tsd.forc.bak"
    
    # Create test config
    create_test_config "test_timestamp"
    
    # Run model and capture output (auto-answer Y)
    print_info "Running model with inconsistent timestamps..."
    echo "Y" | $SHUD_EXECUTABLE test_timestamp > test_inconsistency.log 2>&1
    
    # Check if warning was displayed
    if grep -q "Time Stamp Validation Warnings" test_inconsistency.log; then
        print_success "Warning message displayed correctly"
        
        # Check if it shows the inconsistency
        if grep -q "differs from base" test_inconsistency.log; then
            print_success "Inconsistency details shown"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            print_error "Inconsistency details not shown"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        print_error "Warning message not displayed"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    cleanup_test_env
}

# Test 2: User confirms to continue (Y)
test_user_confirms_continue() {
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    print_header "Test 2: User Confirms to Continue (Y)"
    
    setup_test_env
    
    # Use ccw which already has inconsistent timestamps
    print_info "Running model with auto-confirm (Y)..."
    echo "Y" | $SHUD_EXECUTABLE ccw > test_confirm_y.log 2>&1
    exit_code=$?
    
    # Check if model continued execution
    if [ $exit_code -eq 0 ]; then
        print_success "Model continued execution after user confirmed"
        
        # Check if output files were created
        if [ -d "output/ccw.out" ] && [ $(find output/ccw.out -type f | wc -l) -gt 0 ]; then
            print_success "Output files created successfully"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            print_error "Output files not created"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        print_error "Model did not continue execution"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
}

# Test 3: User declines to continue (n)
test_user_declines_continue() {
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    print_header "Test 3: User Declines to Continue (n)"
    
    # Clear previous outputs
    if [ -d "output/ccw.out" ]; then
        rm -rf output/ccw.out/*
    fi
    
    print_info "Running model with auto-decline (n)..."
    echo "n" | $SHUD_EXECUTABLE ccw > test_confirm_n.log 2>&1
    exit_code=$?
    
    # Check if model terminated
    if grep -q "Model execution terminated by user" test_confirm_n.log; then
        print_success "Model terminated correctly after user declined"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        print_error "Model did not terminate correctly"
        cat test_confirm_n.log
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
}

# Test 4: Check warning message format
test_warning_message_format() {
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    print_header "Test 4: Warning Message Format"
    
    print_info "Checking warning message format from previous test..."
    
    # Use the log from test 1
    if [ -f "test_inconsistency.log" ]; then
        # Check for required elements in warning message
        local has_header=0
        local has_base_date=0
        local has_count=0
        local has_file_list=0
        local has_prompt=0
        
        if grep -q "Time Stamp Validation Warnings" test_inconsistency.log; then
            has_header=1
        fi
        
        if grep -q "Model base date" test_inconsistency.log; then
            has_base_date=1
        fi
        
        if grep -q "Number of inconsistencies:" test_inconsistency.log; then
            has_count=1
        fi
        
        if grep -q "differs from base" test_inconsistency.log; then
            has_file_list=1
        fi
        
        if grep -q "Continue? \[Y\]/n:" test_inconsistency.log; then
            has_prompt=1
        fi
        
        local total=$((has_header + has_base_date + has_count + has_file_list + has_prompt))
        
        if [ $total -eq 5 ]; then
            print_success "Warning message has all required elements"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            print_error "Warning message missing elements (found $total/5)"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        print_error "Test log not found"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
}

# Test 5: Multiple inconsistent files
test_multiple_inconsistent_files() {
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    print_header "Test 5: Multiple Inconsistent Files"
    
    print_info "Running model with multiple inconsistent timestamps..."
    echo "Y" | $SHUD_EXECUTABLE ccw > test_multiple.log 2>&1
    
    # Check if warning shows number of inconsistencies
    if grep -q "Number of inconsistencies:" test_multiple.log; then
        local count=$(grep "Number of inconsistencies:" test_multiple.log | grep -oE '[0-9]+')
        print_success "Found $count inconsistencies reported"
        
        # Check if multiple files are listed
        if [ "$count" -gt 0 ]; then
            print_success "Multiple inconsistencies detected correctly"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            print_error "No inconsistencies detected"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        print_error "Inconsistency count not shown"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
}

# Function to print summary
print_summary() {
    print_header "Test Summary"
    echo -e "Total tests: $TOTAL_TESTS"
    echo -e "${GREEN}Passed: $PASSED_TESTS${NC}"
    echo -e "${RED}Failed: $FAILED_TESTS${NC}"
    
    if [ $FAILED_TESTS -eq 0 ]; then
        print_success "All timestamp inconsistency tests passed!"
        return 0
    else
        print_error "Some tests failed"
        return 1
    fi
}

# Main test execution
main() {
    print_header "SHUD Time Stamp Inconsistency Test Suite"
    print_info "Testing timestamp validation and user confirmation"
    echo ""
    
    # Run all tests
    test_inconsistent_timestamp_warning
    echo ""
    
    test_user_confirms_continue
    echo ""
    
    test_user_declines_continue
    echo ""
    
    test_warning_message_format
    echo ""
    
    test_multiple_inconsistent_files
    echo ""
    
    # Print summary
    print_summary
    
    # Cleanup
    cleanup_test_env
    
    # Return appropriate exit code
    if [ $FAILED_TESTS -eq 0 ]; then
        exit 0
    else
        exit 1
    fi
}

# Run main function
main
