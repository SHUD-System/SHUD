# SHUD 项目改造决策与执行计划

**决策日期**: 2025-11-26  
**决策人**: Lele Shu  
**文档版本**: 1.0

---

## 📊 背景分析

### 项目版本现状

| 版本 | 位置 | 文件数 | 状态 | 说明 |
|------|------|--------|------|------|
| v1.0 | `srcv1/` | 55个C++文件 | 稳定 | 2年前的基础版本 |
| v2.0 | `src/` | 57个C++文件 | 稳定 | 当前主版本，有清晰git历史 |
| v2.1实验版 | `SHUD_improved/` | 63个C++文件 | 实验中 | 前几个月的改进，未提交git |

### SHUD_improved 状态评估

**优点：**
- ✅ 可以编译运行（2025-11-03编译）
- ✅ 包含6个新功能模块
- ✅ 有算法改进和bug修复
- ✅ Makefile重构更清晰

**风险：**
- ❌ 75个文件未暂存修改
- ❌ 84个新文件未跟踪
- ❌ 47个文件有改动（+3763行，-3761行）
- ❌ **没有任何git提交记录**
- ❌ 变量名有回退（说明改动方向不一致）
- ❌ 删除了部分功能（`read_cfgout()`）
- ❌ 无法追溯bug来源
- ❌ 无法回滚错误改动

---

## 🎯 核心决策

### ✅ 决定：从当前版本（src/）重新开始，选择性迁移SHUD_improved的改进

### 决策理由

1. **风险控制**
   - 当前版本是稳定的基线，有清晰的git历史
   - SHUD_improved改动过大且无版本控制，风险不可控
   - 一旦出现问题无法回滚

2. **可维护性**
   - 渐进式迁移可以建立清晰的代码历史
   - 每个改动都有明确的commit记录
   - 便于未来的代码审查和问题追溯

3. **质量保证**
   - 每个功能迁移后都可以单独测试
   - 可以选择性地只迁移有价值的改进
   - 避免引入未知的bug

4. **时间成本**
   - 渐进式迁移：3-4周，质量有保证
   - 直接使用SHUD_improved：看似快，但后续可能花更多时间修bug
   - **长期来看，渐进式迁移更高效**

---

## 📋 执行计划

### 阶段1：准备工作（1天）

**目标**: 保护SHUD_improved的价值，创建清晰的迁移基线

**任务清单**:
- [x] 备份SHUD_improved目录（tar.gz格式）
- [x] 在SHUD_improved中创建git提交保存当前状态
- [x] 创建迁移清单文档
- [x] 创建本决策文档
- [ ] 在主项目创建开发分支 `develop/v2.1-clean`

**预期产出**:
- `SHUD_improved_backup_YYYYMMDD.tar.gz`
- `migration_checklist.md`
- `REFACTORING_PLAN.md`（本文档）
- Git tag: `backup/experimental-changes-20251126`

---

### 阶段2：核心功能迁移（第1-2周）

#### 优先级1：独立的新增模块（第1周）

**迁移顺序**（按依赖关系）:

1. **sundials_includes.h** - SUNDIALS兼容层
   - 位置: `src/sundials_includes.h`
   - 风险: 低（独立文件）
   - 测试: 编译测试

2. **WaterBalance.cpp/hpp** - 水量平衡检查模块
   - 位置: `src/ModelData/WaterBalance.{cpp,hpp}`
   - 风险: 低（独立模块）
   - 测试: 编译测试 + 功能测试

3. **eq_infiltration.cpp/hpp** - 入渗方程模块
   - 位置: `src/Equations/eq_infiltration.{cpp,hpp}`
   - 风险: 中（可能影响核心算法）
   - 测试: 编译测试 + 回归测试

4. **FrostDepth.cpp/hpp** - 冻土深度模块
   - 位置: `src/classes/FrostDepth.{cpp,hpp}`
   - 风险: 低（新功能）
   - 测试: 编译测试

**每个文件的迁移流程**:
```bash
1. 复制文件到目标位置
2. 更新Makefile（如需要）
3. 编译测试: make clean && make shud
4. 运行测试案例
5. Git提交（单独提交，清晰message）
6. 更新migration_checklist.md
```

#### 优先级2：算法改进（第2周）

**重点文件**:
- `src/Equations/functions.cpp/hpp` - 新增函数
  - 导水率计算函数
  - 时间测量函数
  - 错误处理改进

**迁移方法**:
```bash
# 使用diff对比，手动挑选改进
git diff --no-index src/Equations/functions.cpp \
  SHUD_improved/src/Equations/functions.cpp > functions_diff.txt

# 只提取明确的改进，不要整个文件替换
# 每个函数单独commit
```

---

### 阶段3：构建系统改进（第3周）

**任务**:
1. 对比两个Makefile
2. 提取有价值的改进：
   - 增量编译支持
   - 彩色输出
   - 更清晰的目标结构
3. 逐步应用，每次改进后测试编译

**测试重点**:
- macOS编译
- Linux编译（如有环境）
- Debug/Release/OpenMP三种模式

---

### 阶段4：验证与文档（第4周）

**任务清单**:
- [ ] 完整的回归测试（所有测试案例）
- [ ] 性能对比测试
- [ ] 更新 `VersionUpdate.md`
- [ ] 整理 `wiki/` 文档
- [ ] 准备发布说明
- [ ] 代码审查

**验收标准**:
- 所有测试案例通过
- 性能不低于v2.0
- 文档完整
- Git历史清晰

---

## 📦 SHUD_improved 中的可迁移资源

### 🟢 高价值资源（必须迁移）

1. **新增功能模块**（6个文件）
   - `sundials_includes.h` - SUNDIALS兼容层
   - `WaterBalance.cpp/hpp` - 水量平衡检查
   - `eq_infiltration.cpp/hpp` - 入渗方程
   - `FrostDepth.cpp/hpp` - 冻土深度

2. **算法改进**
   - `functions.cpp`: 导水率计算、时间测量、错误处理
   - Bug修复（从VersionUpdate.md）

3. **构建系统**
   - Makefile改进：增量编译、彩色输出

4. **文档系统**
   - `wiki/` 目录的文档

### 🟡 中等价值资源（选择性迁移）

1. **测试框架**
   - `Benchmark/` 目录
   - 测试脚本

2. **Python接口**
   - `python/pySHUD/`
   - `python/pytsd.py`

3. **新测试案例**
   - `input/MC/`
   - `input/WEM/`

### 🔴 需要谨慎处理的改动

1. **变量名回退**
   - `Element.cpp`: `z_surf/z_bottom` → `zmax/zmin`
   - **决定**: 保持新命名（更清晰），不回退

2. **删除的功能**
   - `MD_readin.cpp`: 删除了 `read_cfgout()`
   - **决定**: 先评估影响，确认不需要后再删除

3. **大量数据文件改动**
   - `input/heihe/heihe.cfg.ic` (5006+行)
   - **决定**: 单独评估，可能是测试数据

---

## 🔄 迁移流程规范

### 单个文件迁移标准流程

```bash
# 1. 复制文件
cp SHUD_improved/src/path/to/file.cpp src/path/to/

# 2. 检查依赖
grep -r "file.hpp" src/

# 3. 更新Makefile（如需要）
vim Makefile

# 4. 编译测试
make clean
make shud

# 5. 功能测试
./shud test_case

# 6. Git提交
git add src/path/to/file.{cpp,hpp}
git commit -m "feat: add [功能名称] from SHUD_improved

- 功能描述
- 来源: SHUD_improved/src/path/to/file.cpp
- 测试: [测试结果]
"

# 7. 更新清单
echo "- [x] file.cpp/hpp - $(date)" >> migration_checklist.md
```

### Commit Message 规范

```
类型: 简短描述（不超过50字符）

详细描述：
- 改动内容
- 来源: SHUD_improved/...
- 测试结果
- 相关issue（如有）

类型包括：
- feat: 新功能
- fix: Bug修复
- refactor: 重构
- docs: 文档
- test: 测试
- build: 构建系统
```

---

## ⚠️ 风险管理

### 已识别风险

| 风险 | 等级 | 应对措施 |
|------|------|----------|
| 迁移的代码有bug | 中 | 每个模块单独测试，保持小步提交 |
| 依赖关系复杂 | 中 | 按依赖顺序迁移，先独立模块 |
| 性能下降 | 低 | 每个阶段做性能测试 |
| 时间超期 | 低 | 优先迁移核心功能，次要功能可延后 |

### 回滚策略

- 每个commit都可以独立回滚
- 保留SHUD_improved备份，可随时参考
- 使用git tag标记重要节点

---

## 📊 进度跟踪

### 里程碑

- [ ] **M1**: 阶段1完成（1天）- 2025-11-26
- [ ] **M2**: 核心模块迁移完成（1周）- 2025-12-03
- [ ] **M3**: 算法改进迁移完成（2周）- 2025-12-10
- [ ] **M4**: 构建系统改进完成（3周）- 2025-12-17
- [ ] **M5**: 验证与文档完成（4周）- 2025-12-24
- [ ] **M6**: v2.1正式发布 - 2025-12-31

### 每日检查清单

- [ ] 今天迁移了哪些文件？
- [ ] 编译测试通过？
- [ ] 功能测试通过？
- [ ] Git提交完成？
- [ ] 文档更新？

---

## 📝 参考文档

- `SHUD_improved/CHANGES_SUMMARY.md` - 改动摘要
- `SHUD_improved/CHANGES_COMPARISON.md` - 详细对比
- `SHUD_improved/VersionUpdate.md` - 版本更新日志
- `migration_checklist.md` - 迁移清单（待创建）

---

## 🎯 成功标准

### 技术标准

- ✅ 所有测试案例通过
- ✅ 性能不低于v2.0
- ✅ 代码覆盖率不降低
- ✅ 无内存泄漏
- ✅ 支持macOS/Linux/Windows

### 质量标准

- ✅ Git历史清晰，每个commit可追溯
- ✅ 代码风格一致
- ✅ 文档完整
- ✅ 有完整的changelog

### 时间标准

- ✅ 4周内完成核心功能迁移
- ✅ 6周内完成全部工作并发布

---

## 📞 联系与支持

- **项目维护者**: Lele Shu (shulele@lzb.ac.cn)
- **项目网站**: https://www.shud.xyz
- **GitHub**: https://github.com/SHUD-System/SHUD

---

**最后更新**: 2025-11-26  
**下次审查**: 2025-12-03（第一周结束）
