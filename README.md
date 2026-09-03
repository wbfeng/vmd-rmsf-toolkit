# VMD + TCL 计算蛋白质 RMSF 工具包

> 面向分子动力学模拟后分析：基于 VMD 内置 TCL 解释器，从轨迹文件计算
> **Root Mean Square Fluctuation (RMSF, 均方根涨落)**，输出数据文件、
> 柔性可视化 PDB 与可视化曲线图。

---

## 1. 文件清单

| 文件 | 用途 |
|---|---|
| `rmsf_analysis.tcl` | VMD 主脚本。完成轨迹加载、叠合对齐、RMSF 计算、写出数据文件与 PDB |
| `plot_rmsf.py` | Python 绘图脚本，读取 dat 文件输出 PNG/PDF 曲线图（依赖 matplotlib） |
| `example_rmsf.residue.dat` | 示例数据（模拟 150 残基，用于验证绘图脚本） |
| `example_rmsf.residue.png` | 示例输出图（可作为预期结果参考） |

---

## 2. 工作流总览

```
┌──────────────────────────────────────────────┐
│  1. 手工读取轨迹  (PSF + DCD/XTC/TRR)         │
│  2. source rmsf_analysis.tcl                 │
│  3. ::rmsf::analyze -sel "..." -out prefix   │
│  4. python plot_rmsf.py prefix.residue.dat   │
└──────────────────────────────────────────────┘
   ↓
┌──────────────────────────────────────────────┐
│  输出: prefix.atom.dat   (逐原子 RMSF)        │
│  输出: prefix.residue.dat (逐残基 RMSF)        │
│  输出: prefix_rmsf.pdb  (B-factor = RMSF)     │
│  输出: prefix_average.pdb (平均结构)           │
│  输出: prefix.png / .pdf (曲线图)             │
└──────────────────────────────────────────────┘
```

---

## 3. 详细使用步骤

### Step 1 — 手工读取轨迹（两种方式任选）

**A. VMD 图形界面**
1. `File` → `New Molecule` → 浏览选择 PSF/PDB 作为结构文件
2. `File` → `Load Data Into Molecule` → 选择 DCD/XTC/TRR 轨迹
3. （可选）点击 `Load all frames at once` 一次性加载整段轨迹

**B. VMD TkConsole 命令行**
```tcl
mol new structure.psf type psf waitfor all
mol addfile traj.dcd type dcd first 0 last -1 step 1 waitfor all
mol addfile traj2.dcd type dcd first 0 last -1 step 1 waitfor all   ;# 多段轨迹追加
```

### Step 2 — 载入 TCL 脚本

在 VMD TkConsole 中执行：
```tcl
cd /d  ;# 切到脚本所在目录（若必要）
source rmsf_analysis.tcl
```
成功后控制台会显示：
```
==================================================================
 rmsf_analysis.tcl v1.0 载入成功
 用法示例:
   ::rmsf::analyze -sel "protein and name CA" -out rmsf_result
   ::rmsf::run structure.psf traj.dcd -first 0 -last -1
 可用选项: -sel -fitsel -first -last -step -molid -out
==================================================================
```

### Step 3 — 计算 RMSF

**方式一（推荐，单体系）** —— 显式两步：
```tcl
set molid [::rmsf::load_system structure.psf traj.dcd]
molinfo $molid set top 1
::rmsf::analyze \
    -sel     "protein and name CA"    \
    -fitsel  "protein and backbone"   \
    -first   0                        \
    -last    -1                       \
    -step    1                        \
    -out     rmsf_result
```

**方式二** —— 一键运行（自动加载+计算）：
```tcl
::rmsf::run structure.psf traj.dcd -sel "protein and name CA" -out rmsf_result
```

**方式三** —— 多段轨迹串联：
```tcl
::rmsf::run structure.psf {traj_1.dcd traj_2.dcd traj_3.dcd} -out rmsf_result
```

**方式四** —— 已加载分子直接计算（molid = top）：
```tcl
::rmsf::analyze -out rmsf_result
```

**常用选择语句示例**

| 场景 | `-sel` | `-fitsel` |
|---|---|---|
| 标准 Cα 残基 RMSF | `protein and name CA` | `protein and backbone` |
| 全原子逐残基 | `protein and noh` | `protein and backbone` |
| 配体口袋残基 | `protein and name CA and within 5 of resname LIG` | `protein and backbone` |
| 跨链对齐 | `segid PROA and name CA` | `segid PROA and backbone` |
| 全受体 | `protein` | `protein and backbone` |

### Step 4 — 绘制 RMSF 曲线

激活已安装 matplotlib 的 Python 环境（推荐虚拟环境）：
```bash
# 一次性安装（如未安装）
pip install matplotlib

# 绘制（自动识别 5列/6列 dat 格式）
python plot_rmsf.py rmsf_result.residue.dat
python plot_rmsf.py rmsf_result.residue.dat -o myfig --title "Receptor X" --sigma 1.5
python plot_rmsf.py rmsf_result.atom.dat  --no-highlight
```

脚本输出：
- `*_result.png` / `*_result.pdf`：曲线图
- 控制台摘要：均值/标准差/柔性残基数/最大RMSF残基

---

## 4. 输出文件含义

### `prefix.atom.dat`
```
# 逐原子 RMSF (单位: Angstrom)
# 列: resid resname chain segid atomname rmsf
     1   MET  A  PROA  N      0.8721
     1   MET  A  PROA  CA     0.7645
     ...
```

### `prefix.residue.dat`
```
# 逐残基 RMSF (单位: Angstrom)
# 列: resid resname chain segid rmsf
     1   MET  A  PROA  0.8123
     2   LYS  A  PROA  0.7451
     ...
```
残基 RMSF = 该残基所有被选原子的 RMSF 算术平均。

### `prefix_rmsf.pdb`
B-factor 列写入 RMSF 值（Å）。在 VMD 中：
```
Graphics → Representations → Coloring Method → Beta
```
即可在卡通/球棍/Cα 模式上以颜色梯度显示 RMSF 大小（蓝→白→红 表示 低→中→高 柔性）。

### `prefix_average.pdb`
对齐后的平均结构坐标（可直接叠合到原始结构做结构比较）。

---

## 5. 算法说明

1. **叠合对齐** —— 每一帧相对起始帧（`first`）以参考选择（`fitsel`）做
   最小二乘拟合：`measure fit`，消除平移/旋转；
2. **平均结构** —— 对齐后坐标按原子求算术平均：`measure avpos`；
3. **逐原子 RMSF** —— $\mathrm{RMSF}_i = \sqrt{ \frac{1}{N}\sum_{f} |\mathbf{r}_i(f) - \langle\mathbf{r}_i\rangle|^2 }$
4. **逐残基聚合** —— 同残基内原子 RMSF 算术平均；
5. **单位** —— 全部输出使用 Å（埃）。

---

## 6. 常见问题

**Q1. 报错 `0 atoms selected for selection '...'`**
- 检查分子是否已加载：`molinfo top` 看 molid
- 确认选择语法与结构类型匹配（如 PRB/全原子）
- 对 GROMACS 输出需先确认 `moltype`（一般 type auto 即可）

**Q2. 计算很慢**
- RMSF 计算复杂度 ≈ `O(选择原子数 × 帧数)`，CA 选择比全原子快 5-10 倍
- 增大 `-step`（如 2 或 5）可线性加速但损失时间分辨率
- 对于超大蛋白 (>1000 残基)，建议仅做 Cα 估算柔性

**Q3. 想做 RMSD 而不是 RMSF**
- 直接用 VMD 工具：`Extensions → Analysis → RMSD Trajectory Tool`
- 本脚本专用于 RMSF（依赖与平均结构的比较，RMSE 等价但概念上不同）

**Q4. 对齐修改了原始坐标是否安全？**
- 仅修改内存中的坐标，不写回原始 DCD 文件
- 如需保留原始未对齐坐标，计算前可先 `mol copy` 复制分子
- 或者输出对齐后的新轨迹：
  ```tcl
  animate write dcd aligned.dcd waitfor all sel [atomselect top "all"]
  ```

**Q5. RMSF 与 B-factor 转换**
- B-factor（晶体学温度因子）与 RMSF 关系：$B_i = \frac{8\pi^2}{3}\mathrm{RMSF}_i^2$
- 本脚本的 `*_rmsf.pdb` 已写入 RMSF 到 B-factor 列（不是 B），与文献中的
  "B-factor coloring" 不严格等价。如需输出真实 B 因子：
  ```tcl
  set B []
  foreach r $rmsf { lappend B [expr {8.0 * 3.14159265 * 3.14159265 / 3.0 * $r * $r}] }
  $sel set beta $B
  ```

**Q6. 多个链 / 复合物**
- `-sel "protein and chain X and name CA"` 按链筛选
- `-fitsel "protein and backbone"` 一般是全链 backbone
- 输出文件的 chain/segid 列保留全部标识，便于后处理拆分

---

## 7. 扩展建议

- **时序 RMSF (rolling RMSF)**：可在脚本基础上加入滑动窗口，每 100 帧输出一次局部 RMSF
- **残基相关性矩阵 (dynamic cross-correlation)**：用 `measure dipole` 等工具与本脚本联用
- **PCA / 本征投影**：用 `measure modebasis` 做主成分分析，结合 RMSF 找到柔性方向
- **X-ray B-factor 验证**：若同一蛋白有 PDB X-ray 结构，可与其 B 因子做皮尔森相关，验证模拟合理性

---

## 8. 参考与引用

- Humphrey W, Dalke A, Schulten K. *VMD — Visual Molecular Dynamics.* J Mol Graph 1996, 14:33-38
- Roe DR, Cheatham TE III. *PTRAJ and CPPTRAJ: Software for Processing and Analysis of Molecular Dynamics Trajectory Data.* J Chem Theory Comput 2013, 9:3084-3095
- 通用 RMSF 公式参考: Kearsley SK. *On the orthogonal transformation used for structural comparisons.* Acta Cryst A 1989, 45:208-210

---

_版本: 1.0 | 作者: RMSF analysis toolkit | 许可: 自由学术使用_
