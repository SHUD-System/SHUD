# =============================================================
# SHUD Hydrological Model Makefile
# SHUD 水文模型 Makefile
# Version: 2.0
# Author: Lele Shu (lele.shu@gmail.com)
# Last Update: 2024-07
# =============================================================
#
# 说明/Description：
# 1. 支持 Release/Debug/OpenMP 三种编译模式
#    Support Release/Debug/OpenMP build modes
# 2. 变量集中管理，便于移植和维护
#    Centralized variables for easy maintenance and portability
# 3. 伪目标声明，结构清晰
#    Use .PHONY for clear structure
# 4. 美化输出，便于阅读
#    Colorful output for better readability
# 5. 所有目标均采用分文件增量编译，所有.o文件集中于build/目录，且目录结构与src/镜像
#    All targets use per-file incremental build, all .o files in build/ mirroring src/
# =============================================================

# ========== 用户可配置部分 / User Configurable Section ==========
# SUNDIALS 安装路径（需自行修改）
# SUNDIALS install path (modify as needed)
SUNDIALS_DIR ?= $(HOME)/sundials
# OpenMP 路径（如需并行支持）
# OpenMP path (if parallel support is needed)
INC_OMP      ?= /usr/local/opt/libomp/include
LIB_OMP      ?= /usr/local/opt/libomp/lib
# 系统库路径
# System library path
LIB_SYS      ?= /usr/local/lib/
# 源码路径
# Source code directory
SRC_DIR      ?= src
# 构建输出路径
# Build output directory
BUILDDIR     ?= .
# 对象文件输出目录
OBJDIR       ?= build

# ========== 编译器与选项 / Compiler and Flags ==========
CXX          ?= g++
# 依赖关系跟踪选项 / Dependency tracking flags
DEPFLAGS     = -MMD -MP
CXXFLAGS     ?= -std=c++14 -O3 -g
CXXFLAGS_DBG ?= -std=c++14 -O0 -g -DDEBUG
CXXFLAGS_OMP ?= -std=c++14 -O3 -g -D_OPENMP_ON -fopenmp
LDFLAGS      ?= -lm -lsundials_cvode -lsundials_nvecserial
LDFLAGS_OMP  ?= -Xpreprocessor -fopenmp -lomp -lsundials_nvecopenmp

INCLUDES     = -I$(SUNDIALS_DIR)/include \
               -I$(INC_OMP) \
               -I$(SRC_DIR)/Model \
               -I$(SRC_DIR)/ModelData \
               -I$(SRC_DIR)/classes \
               -I$(SRC_DIR)/Equations
LIBRARIES    = -L$(LIB_OMP) -L$(SUNDIALS_DIR)/lib -L$(LIB_SYS)
RPATH        = -Wl,-rpath,$(SUNDIALS_DIR)/lib

# ========== 颜色美化 / Color Output ==========
C_RESET  = \033[0m
C_GREEN  = \033[32m
C_YELLOW = \033[33m
C_BLUE   = \033[34m
C_RED    = \033[31m

# ========== 源文件收集 / Source Files ==========
SRC_MAIN := $(SRC_DIR)/main.cpp
SRC_CPP := $(shell find $(SRC_DIR) -name '*.cpp')
ALL_CPP := $(SRC_CPP)

# ========== 生成镜像结构的.o文件路径 / Mirror structure for .o files ==========
OBJ_CPP      := $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR)/%.o, $(SRC_CPP))
ALL_OBJ      := $(OBJ_CPP)
OBJ_CPP_OMP  := $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR)/%.omp.o, $(SRC_CPP))
OMP_OBJ      := $(OBJ_CPP_OMP)
OBJ_CPP_DBG  := $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR)/%.dbg.o, $(SRC_CPP))
DBG_OBJ      := $(OBJ_CPP_DBG)
# 依赖关系文件 / Dependency files
DEP_CPP      := $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR)/%.d, $(SRC_CPP))
DEP_CPP_OMP  := $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR)/%.omp.d, $(SRC_CPP))
DEP_CPP_DBG  := $(patsubst $(SRC_DIR)/%.cpp, $(OBJDIR)/%.dbg.d, $(SRC_CPP))

# ========== 目标文件 / Targets ==========
TARGET        = $(BUILDDIR)/shud
TARGET_OMP    = $(BUILDDIR)/shud_omp
TARGET_DEBUG  = $(BUILDDIR)/shud_debug

# ========== 颜色美化 / Color Output ==========
C_RESET  = \033[0m
C_GREEN  = \033[32m
C_YELLOW = \033[33m
C_BLUE   = \033[34m
C_RED    = \033[31m

# ========== 默认目标 / Default Targets ==========
.PHONY: all help clean shud shud_omp shud_debug check cvode

all: clean shud shud_omp
	@echo
	@echo "${C_GREEN} All builds finished!${C_RESET}"
	@echo ""	

help:
	@echo "${C_BLUE}用法/Usage:${C_RESET}"
	@echo "  make all        - 编译 Release 和 OpenMP 版本 / Build Release and OpenMP versions"
	@echo "  make shud       - 编译 Release 版本 / Build Release version"
	@echo "  make shud_omp   - 编译 OpenMP 并行版本 / Build OpenMP parallel version"
	@echo "  make shud_debug - 编译 Debug 版本 / Build Debug version"
	@echo "  make clean      - 清理所有目标文件 / Clean all targets"
	@echo "  make check      - 检查 SUNDIALS 安装 / Check SUNDIALS installation"
	@echo "  make cvode      - 安装 SUNDIALS/CVODE 到 ~/sundials / Install SUNDIALS/CVODE to ~/sundials"
	@echo ""
	@echo ""

# ========== 自动创建build子目录 / Auto-create build subdirs ==========
$(OBJDIR):
	@mkdir -p $(OBJDIR)
$(OBJDIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJDIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@
$(OBJDIR)/%.omp.o: $(SRC_DIR)/%.cpp | $(OBJDIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS_OMP) $(DEPFLAGS) $(INCLUDES) -c $< -o $@
$(OBJDIR)/%.dbg.o: $(SRC_DIR)/%.cpp | $(OBJDIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS_DBG) $(DEPFLAGS) -fsanitize=address $(INCLUDES) -c $< -o $@

# ========== 包含依赖关系文件 / Include dependency files ==========
-include $(DEP_CPP)
-include $(DEP_CPP_OMP)
-include $(DEP_CPP_DBG)

# ========== 分文件编译规则 / Per-file Build Rules ==========
# 1. Release
$(TARGET): $(ALL_OBJ)
	@echo "${C_YELLOW}Linking $@ ...${C_RESET}"
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $@ $(ALL_OBJ) $(LDFLAGS)
	@echo "${C_GREEN}Build $@ success!${C_RESET}"

# 2. OpenMP
$(TARGET_OMP): $(OMP_OBJ)
	@echo "${C_YELLOW}Linking(OpenMP) $@ ...${C_RESET}"
	$(CXX) $(CXXFLAGS_OMP) $(INCLUDES) $(LIBRARIES) $(RPATH) -o $@ $(OMP_OBJ) $(LDFLAGS) $(LDFLAGS_OMP)
	@echo "${C_GREEN}Build $@ success!${C_RESET}"

# 3. Debug
$(TARGET_DEBUG): $(DBG_OBJ)
	@echo "${C_YELLOW}Linking(Debug) $@ ...${C_RESET}"
	$(CXX) $(CXXFLAGS_DBG) -fsanitize=address $(INCLUDES) $(LIBRARIES) $(RPATH) -o $@ $(DBG_OBJ) $(LDFLAGS)
	@echo "${C_GREEN} Build $@ success!${C_RESET}"

# ========== 依赖检查与安装 / Dependency Check & Install ==========
check:
	@echo "${C_BLUE}Checking SUNDIALS directory: ${C_RESET} $(SUNDIALS_DIR)"
	@ls $(SUNDIALS_DIR)
	@ls $(SUNDIALS_DIR)/lib
	@echo "${C_GREEN}SUNDIALS check finished.${C_RESET}"
	@echo ""

cvode:
	@echo "${C_BLUE} Installing SUNDIALS/CVODE...${C_RESET}"
	chmod +x configure
	./configure
	@echo "${C_GREEN}SUNDIALS/CVODE install finished.${C_RESET}"
	@echo ""

# ========== 清理 / Clean ==========
clean:
	@echo "${C_RED}Cleaning targets and objects...${C_RESET}"
	@rm -rf $(OBJDIR)
	@rm -f $(TARGET) $(TARGET_OMP) $(TARGET_DEBUG)
	@echo "${C_GREEN}Clean finished.${C_RESET}"





