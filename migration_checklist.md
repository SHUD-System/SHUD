# SHUD v2.1 迁移清单

**开始日期**: 2025-11-26  
**目标**: 从 SHUD_improved 选择性迁移改进到干净的 v2.1 分支  
**原则**: 小步提交，逐个测试，保持git历史清晰

---

## 📦 阶段1：准备工作

- [x] 创建备份: `SHUD_improved_backup_20251126.tar.gz` (792MB)
- [x] SHUD_improved创建git提交保存状态
- [x] 创建git tag: `backup/experimental-changes-20251126`
- [x] 创建决策文档: `REFACTORING_PLAN.md`
- [x] 创建本迁移清单: `migration_checklist.md`
- [ ] 创建开发分支: `develop/v2.1-clean`
- [ ] 确认当前版本可以编译运行

**完成时间**: 2025-11-26

---

## 📦 阶段2：核心功能迁移（第1-2周）

### 优先级1：独立新增模块

#### 1. SUNDIALS兼容层
- [ ] 复制 `sundials_includes.h` 到 `src/`
- [ ] 检查是否有文件引用此头文件
- [ ] 编译测试: `make clean && make shud`
- [ ] Git提交
- [ ] 更新日期: ___________

**来源**: `SHUD_improved/src/sundials_includes.h`  
**风险等级**: 🟢 低  
**依赖**: 无

---

#### 2. 水量平衡检查模块
- [ ] 复制 `WaterBalance.cpp` 到 `src/ModelData/`
- [ ] 复制 `WaterBalance.hpp` 到 `src/ModelData/`
- [ ] 更新 `Makefile` 添加编译目标
- [ ] 检查 `Model_Data.hpp` 中的 `CheckWaterBalance()` 函数声明
- [ ] 编译测试
- [ ] 功能测试（运行测试案例）
- [ ] Git提交
- [ ] 更新日期: ___________

**来源**: `SHUD_improved/src/ModelData/WaterBalance.{cpp,hpp}`  
**风险等级**: 🟢 低  
**依赖**: 可能需要 `Model_Data.hpp` 中的函数声明

---

#### 3. 入渗方程模块
- [ ] 复制 `eq_infiltration.cpp` 到 `src/Equations/`
- [ ] 复制 `eq_infiltration.hpp` 到 `src/Equations/`
- [ ] 更新 `Makefile` 添加编译目标
- [ ] 检查依赖关系（是否被其他模块调用）
- [ ] 编译测试
- [ ] 功能测试（对比结果）
- [ ] 回归测试（确保不影响现有功能）
- [ ] Git提交
- [ ] 更新日期: ___________

**来源**: `SHUD_improved/src/Equations/eq_infiltration.{cpp,hpp}`  
**风险等级**: 🟡 中（可能影响核心算法）  
**依赖**: 需要检查是否被 `MD_*.cpp` 调用

---

#### 4. 冻土深度模块
- [ ] 复制 `FrostDepth.cpp` 到 `src/classes/`
- [ ] 复制 `FrostDepth.hpp` 到 `src/classes/`
- [ ] 更新 `Makefile` 添加编译目标
- [ ] 检查是否需要在其他类中集成
- [ ] 编译测试
- [ ] 功能测试
- [ ] Git提交
- [ ] 更新日期: ___________

**来源**: `SHUD_improved/src/classes/FrostDepth.{cpp,hpp}`  
**风险等级**: 🟢 低（新功能，独立模块）  
**依赖**: 可能需要在 `Element` 或 `Model_Data` 中集成

---

### 优先级2：算法改进（第2周）

#### 5. functions.cpp 新增函数
- [ ] 对比两个版本的 `functions.cpp`
  ```bash
  git diff --no-index src/Equations/functions.cpp \
    SHUD_improved/src/Equations/functions.cpp > /tmp/functions_diff.txt
  ```
- [ ] 提取新增函数列表：
  - [ ] `getSecond()` - 时间测量函数
  - [ ] `Kcoeff_of_h()` - 导水率系数计算
  - [ ] `mean_Kcoeff_between_depths()` - 深度区间平均导水率
  - [ ] 其他改进（错误处理等）
- [ ] 逐个函数迁移（每个函数单独commit）
- [ ] 更新 `functions.hpp` 添加函数声明
- [ ] 编译测试
- [ ] 单元测试（如有）
- [ ] Git提交（每个函数一个commit）
- [ ] 更新日期: ___________

**来源**: `SHUD_improved/src/Equations/functions.cpp`  
**风险等级**: 🟡 中  
**依赖**: 需要检查函数调用关系

---

#### 6. 其他算法文件改进
- [ ] `cvode_config.cpp` - 对比并提取改进
- [ ] `funPlatform.cpp` - 对比并提取改进
- [ ] `is_sm_et.cpp` - 对比并提取改进
- [ ] `print.cpp` - 对比并提取改进
- [ ] 每个改进单独commit
- [ ] 更新日期: ___________

**方法**: 使用 `git diff --no-index` 对比，手动挑选改进

---

## 📦 阶段3：构建系统改进（第3周）

### 7. Makefile 重构
- [ ] 对比两个 Makefile
  ```bash
  diff -u Makefile SHUD_improved/Makefile > /tmp/makefile_diff.txt
  ```
- [ ] 提取改进点：
  - [ ] 增量编译支持（build/目录结构）
  - [ ] 彩色输出
  - [ ] 更清晰的目标结构
  - [ ] 注释和文档
- [ ] 逐步应用改进
- [ ] 测试三种编译模式：
  - [ ] Release: `make shud`
  - [ ] Debug: `make shud_debug`
  - [ ] OpenMP: `make shud_omp`
- [ ] Git提交
- [ ] 更新日期: ___________

**来源**: `SHUD_improved/Makefile`  
**风险等级**: 🟡 中  
**依赖**: 影响整个构建流程

---

## 📦 阶段4：文档与测试（第4周）

### 8. 文档迁移
- [ ] 复制 `wiki/` 目录
  - [ ] `Computation.md`
  - [ ] `WaterBalance.md`
  - [ ] `Output.md`
  - [ ] `equation.md`
  - [ ] `function.md`
  - [ ] `Conductivity.Rmd/html`
- [ ] 审查和更新文档内容
- [ ] Git提交
- [ ] 更新日期: ___________

---

### 9. 测试案例（可选）
- [ ] 评估是否需要新测试案例
  - [ ] `input/MC/` - MC案例
  - [ ] `input/WEM/` - WEM案例
- [ ] 如需要，复制并测试
- [ ] Git提交
- [ ] 更新日期: ___________

---

### 10. Python接口（可选）
- [ ] 评估是否需要Python接口
  - [ ] `python/pySHUD/`
  - [ ] `python/pytsd.py`
- [ ] 如需要，复制并测试
- [ ] Git提交
- [ ] 更新日期: ___________

---

### 11. 测试框架（可选）
- [ ] 评估是否需要测试框架
  - [ ] `Benchmark/`
  - [ ] `catch.hpp`
- [ ] 如需要，复制并测试
- [ ] Git提交
- [ ] 更新日期: ___________

---

## 📦 阶段5：验证与发布

### 12. 完整测试
- [ ] 编译所有目标
  - [ ] `make shud`
  - [ ] `make shud_debug`
  - [ ] `make shud_omp`
- [ ] 运行所有测试案例
  - [ ] ccw案例
  - [ ] heihe案例
  - [ ] qhh案例
  - [ ] 其他案例
- [ ] 性能对比测试（vs v2.0）
- [ ] 内存泄漏检查
- [ ] 更新日期: ___________

---

### 13. 文档更新
- [ ] 更新 `VersionUpdate.md`
  - [ ] 列出所有新功能
  - [ ] 列出所有bug修复
  - [ ] 列出所有改进
- [ ] 更新 `README.md`
- [ ] 创建 `CHANGELOG.md`
- [ ] Git提交
- [ ] 更新日期: ___________

---

### 14. 代码审查
- [ ] 审查所有新增代码
- [ ] 检查代码风格一致性
- [ ] 检查注释完整性
- [ ] 检查错误处理
- [ ] 更新日期: ___________

---

### 15. 发布准备
- [ ] 创建 release 分支
- [ ] 打 tag: `v2.1.0`
- [ ] 准备发布说明
- [ ] 更新网站文档
- [ ] 更新日期: ___________

---

## 🚫 不迁移的内容（需要重新实现或放弃）

### 变量名回退（不迁移）
- ❌ `Element.cpp` 中的 `zmax/zmin` 回退
- **决定**: 保持 `z_surf/z_bottom` 命名（更清晰）

### 删除的功能（需要评估）
- ⚠️ `MD_readin.cpp` 中删除的 `read_cfgout()` 函数
- **决定**: 先评估影响，确认不需要后再删除

### 大量数据文件改动（需要单独评估）
- ⚠️ `input/heihe/heihe.cfg.ic` (5006+行改动)
- **决定**: 单独评估，可能是测试数据更新

---

## 📊 进度统计

**总任务数**: 15个主要任务  
**已完成**: 5个（阶段1）  
**进行中**: 0个  
**待开始**: 10个  

**完成百分比**: 33%

---

## 📝 每日日志

### 2025-11-26
- ✅ 创建备份 `SHUD_improved_backup_20251126.tar.gz`
- ✅ SHUD_improved 创建git提交（commit: cf51553）
- ✅ 创建git tag: `backup/experimental-changes-20251126`
- ✅ 创建 `REFACTORING_PLAN.md`
- ✅ 创建 `migration_checklist.md`
- ⏳ 下一步: 创建开发分支 `develop/v2.1-clean`

---

### 2025-11-27
- [ ] 待记录...

---

## 🎯 里程碑

- [x] **M0**: 准备工作完成 - 2025-11-26
- [ ] **M1**: 核心模块迁移完成 - 目标: 2025-12-03
- [ ] **M2**: 算法改进迁移完成 - 目标: 2025-12-10
- [ ] **M3**: 构建系统改进完成 - 目标: 2025-12-17
- [ ] **M4**: 验证与文档完成 - 目标: 2025-12-24
- [ ] **M5**: v2.1正式发布 - 目标: 2025-12-31

---

## 📞 问题记录

### 待解决问题
1. 是否需要迁移Python接口？
2. 是否需要迁移新测试案例（MC, WEM）？
3. `read_cfgout()` 功能是否真的不需要？

### 已解决问题
- 无

---

**最后更新**: 2025-11-26  
**更新人**: Lele Shu
