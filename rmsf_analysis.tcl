# ======================================================================
# rmsf_analysis.tcl
# ------------------------------------------------------------------
# VMD TCL 脚本：计算蛋白质分子动力学模拟轨迹的 RMSF
#   (Root Mean Square Fluctuation, 均方根涨落)
#
# 算法流程：
#   (1) 逐帧叠合对齐（消除整体平移/旋转，默认以 backbone 为参考）
#   (2) 计算对齐后坐标的平均结构 (measure avpos)
#   (3) 逐原子累计平方偏差 -> RMSF = sqrt(<|r_i - <r_i>|^2>)
#   (4) 按残基聚合（同一残基内原子取平均）
#   (5) 输出数据文件 + B-factor 替换的 PDB（可直接在 VMD 中按 B-factor 着色）
#
# 使用方法（VMD TkConsole，两步）：
#   Step 1 - 手工读取轨迹（二选一）：
#       A. GUI: File -> New Molecule 载入 PSF/PDB，再 File -> Load Data
#              Into Molecule 载入 DCD/XTC/TRR
#       B. 命令:
#           mol new structure.psf type psf waitfor all
#           mol addfile traj.dcd type dcd first 0 last -1 step 1 waitfor all
#   Step 2 - 载入并运行脚本：
#       source rmsf_analysis.tcl
#       ::rmsf::analyze -sel "protein and name CA" -out rmsf_result
#       （或一键运行: ::rmsf::run structure.psf traj.dcd）
#
# 依赖: 仅需 VMD >= 1.9.2 自带 TCL 解释器
# 作者: RMSF analysis toolkit
# ======================================================================

namespace eval ::rmsf {
    variable version 1.0
}

# ------------------------------------------------------------------
# 辅助函数：载入体系（拓扑 + 单条或多条轨迹）
#   topology : PSF / PDB 文件路径
#   traj     : 轨迹文件路径或路径列表（DCD/XTC/TRR 均可，type 自动识别）
#   返回 molid
# ------------------------------------------------------------------
proc ::rmsf::load_system {topology traj} {
    if {![file exists $topology]} {
        error "拓扑文件不存在: $topology"
    }
    set molid [mol new $topology type auto waitfor all]
    puts "\[RMSF\] 载入拓扑: $topology (molid=$molid)"

    if {$traj ne ""} {
        if {[llength $traj] == 1} {
            set trajlist [list $traj]
        } else {
            set trajlist $traj
        }
        foreach f $trajlist {
            if {![file exists $f]} {
                error "轨迹文件不存在: $f"
            }
            mol addfile $f type auto first 0 last -1 step 1 waitfor all $molid
            puts "\[RMSF\] 载入轨迹: $f"
        }
    }
    set nf [molinfo $molid get numframes]
    puts "\[RMSF\] 共 $nf 帧已载入 (molid=$molid)"
    return $molid
}

# ------------------------------------------------------------------
# 核心函数：RMSF 计算
# 参数（-key value 形式）：
#   -sel    : 计算RMSF的原子选择，默认 "protein and name CA"（逐残基）
#             若要全原子: "protein and noh"
#   -fitsel : 叠合对齐参考原子，默认 "protein and backbone"
#   -first  : 起始帧（含），默认 0
#   -last   : 结束帧（含），-1 表示最后一帧
#   -step   : 帧间隔，默认 1
#   -molid  : 分子ID，"top" 表示当前顶层分子
#   -out    : 输出文件前缀，默认 "rmsf"
# 输出：
#   <out>.atom.dat     逐原子RMSF (resid resname chain segid name rmsf)
#   <out>.residue.dat  逐残基RMSF (resid resname chain segid rmsf)
#   <out>_rmsf.pdb     B-factor列写入RMSF的PDB（VMD按B-factor着色可视化）
#   <out>_average.pdb  平均结构坐标PDB
# 返回：逐原子RMSF列表
# ------------------------------------------------------------------
proc ::rmsf::analyze {args} {
    # ---------- 默认参数解析 ----------
    set opts(-sel)    "protein and name CA"
    set opts(-fitsel) "protein and backbone"
    set opts(-first)  0
    set opts(-last)   -1
    set opts(-step)   1
    set opts(-molid)  "top"
    set opts(-out)    "rmsf"

    foreach {key val} $args {
        if {![info exists opts($key)]} {
            error "未知选项: $key （可用: -sel -fitsel -first -last -step -molid -out）"
        }
        set opts($key) $val
    }

    set molid $opts(-molid)
    if {$molid eq "top"} { set molid [molinfo top] }

    set nf   [molinfo $molid get numframes]
    set f0   $opts(-first)
    set f1   $opts(-last)
    if {$f1 < 0 || $f1 >= $nf} { set f1 [expr {$nf - 1}] }
    if {$f0 < 0}    { set f0 0 }
    if {$f0 > $f1}  { error "起始帧($f0)大于结束帧($f1)" }
    set step $opts(-step)
    if {$step < 1} { set step 1 }

    # 实际参与统计的帧数
    set nframes 0
    for {set f $f0} {$f <= $f1} {incr f $step} { incr nframes }
    if {$nframes < 2} {
        error "参与计算的帧数不足 (仅 $nframes 帧)，RMSF无统计意义"
    }
    puts "\[RMSF\] 分子ID=$molid, 帧范围=($f0,$f1), 步长=$step, 共 $nframes 帧参与统计"

    # ---------- 原子选择 ----------
    set sel [atomselect $molid $opts(-sel)]
    if {[$sel num] == 0} {
        $sel delete
        error "选择 '$opts(-sel)' 未选中任何原子，请检查选择语句"
    }
    puts "\[RMSF\] RMSF计算选择: '$opts(-sel)' -> [$sel num] 个原子"

    # ==================================================================
    # Step 1: 逐帧叠合对齐（最小二乘拟合到起始帧，消除整体运动）
    # 注意：此操作会修改内存中的轨迹坐标（不写回磁盘文件），
    #       如需保留原始坐标请先另存轨迹或复制分子。
    # ==================================================================
    puts "\[RMSF\] Step 1/4: 逐帧叠合对齐 (参考选择: '$opts(-fitsel)') ..."
    set ref [atomselect $molid $opts(-fitsel) frame $f0]
    if {[$ref num] == 0} {
        $sel delete; $ref delete
        error "对齐参考选择 '$opts(-fitsel)' 未选中任何原子"
    }
    set iframe 0
    for {set f $f0} {$f <= $f1} {incr f $step} {
        set cur [atomselect $molid $opts(-fitsel) frame $f]
        set mat [measure fit $cur $ref]
        set allatoms [atomselect $molid all frame $f]
        $allatoms move $mat
        $cur delete
        $allatoms delete
        incr iframe
        if {[expr {$iframe % 50}] == 0} {
            puts "          已对齐 $iframe / $nframes 帧"
        }
    }
    $ref delete
    puts "\[RMSF\] 对齐完成"

    # ==================================================================
    # Step 2: 平均结构（对齐后坐标的逐原子算术平均）
    # ==================================================================
    puts "\[RMSF\] Step 2/4: 计算平均结构 ..."
    set avpos [measure avpos $sel first $f0 last $f1 step $step]

    # ==================================================================
    # Step 3: 逐原子累计平方偏差 -> RMSF
    #   RMSF_i = sqrt( (1/N) * sum_f |r_i(f) - <r_i>|^2 )
    # ==================================================================
    puts "\[RMSF\] Step 3/4: 累计平方偏差并计算RMSF ..."
    set natoms [$sel num]
    set sumsq {}
    for {set i 0} {$i < $natoms} {incr i} { lappend sumsq 0.0 }

    set iframe 0
    for {set f $f0} {$f <= $f1} {incr f $step} {
        $sel frame $f
        set coords [$sel get {x y z}]
        for {set i 0} {$i < $natoms} {incr i} {
            set c [lindex $coords $i]
            set a [lindex $avpos $i]
            set dx [expr {[lindex $c 0] - [lindex $a 0]}]
            set dy [expr {[lindex $c 1] - [lindex $a 1]}]
            set dz [expr {[lindex $c 2] - [lindex $a 2]}]
            lset sumsq $i [expr {[lindex $sumsq $i] + $dx*$dx + $dy*$dy + $dz*$dz}]
        }
        incr iframe
    }

    set rmsf {}
    foreach s $sumsq {
        lappend rmsf [expr {sqrt($s / double($nframes))}]
    }

    # 统计信息
    set rmin  1e30; set rmax -1e30; set rsum 0.0
    foreach r $rmsf {
        if {$r < $rmin} { set rmin $r }
        if {$r > $rmax} { set rmax $r }
        set rsum [expr {$rsum + $r}]
    }
    set ravg [expr {$rsum / double($natoms)}]
    puts "\[RMSF\] 逐原子RMSF统计 (单位: Angstrom):"
    puts [format "          min = %.4f | max = %.4f | mean = %.4f" $rmin $rmax $ravg]

    # ==================================================================
    # Step 4: 输出文件
    # ==================================================================
    puts "\[RMSF\] Step 4/4: 写出数据文件 ..."

    # --- 4a. 逐原子数据文件 ---
    set resids   [$sel get resid]
    set resnames [$sel get resname]
    set chains   [$sel get chain]
    set segids   [$sel get segid]
    set names    [$sel get name]

    set fp [open "${opts(-out)}.atom.dat" w]
    puts $fp "# 逐原子 RMSF (单位: Angstrom)"
    puts $fp "# 体系: [molinfo $molid get name] | 帧范围: $f0-$f1 步长 $step | 帧数: $nframes"
    puts $fp "# 列: resid resname chain segid atomname rmsf"
    for {set i 0} {$i < $natoms} {incr i} {
        puts $fp [format "%6d %4s %2s %6s %4s %10.4f" \
            [lindex $resids $i] [lindex $resnames $i] [lindex $chains $i] \
            [lindex $segids $i] [lindex $names $i] [lindex $rmsf $i]]
    }
    close $fp
    puts "\[RMSF\] 已写出: ${opts(-out)}.atom.dat"

    # --- 4b. 逐残基聚合（同残基内原子RMSF取平均） ---
    array unset racc
    array set racc {}
    set keyorder {}
    for {set i 0} {$i < $natoms} {incr i} {
        set key "[lindex $segids $i]:[lindex $chains $i]:[lindex $resids $i]"
        if {![info exists racc($key)]} {
            set racc($key) [list 0.0 0]
            lappend keyorder $key
        }
        set s [lindex $racc($key)]
        lset racc($key) 0 [expr {[lindex $s 0] + [lindex $rmsf $i]}]
        lset racc($key) 1 [expr {[lindex $s 1] + 1}]
    }

    set fp [open "${opts(-out)}.residue.dat" w]
    puts $fp "# 逐残基 RMSF (单位: Angstrom)"
    puts $fp "# 体系: [molinfo $molid get name] | 帧范围: $f0-$f1 步长 $step | 帧数: $nframes"
    puts $fp "# 列: resid resname chain segid rmsf"
    foreach key $keyorder {
        set s $racc($key)
        set mean [expr {[lindex $s 0] / double([lindex $s 1])}]
        set seg  [lindex [split $key ":"] 0]
        set chn  [lindex [split $key ":"] 1]
        set rid  [lindex [split $key ":"] 2]
        # 残基名取第一个匹配原子
        set idx [lsearch -exact $resids $rid]
        set rname [lindex $resnames $idx]
        puts $fp [format "%6d %4s %2s %6s %10.4f" $rid $rname $chn $seg $mean]
    }
    close $fp
    puts "\[RMSF\] 已写出: ${opts(-out)}.residue.dat"

    # --- 4c. B-factor 替换 PDB（用于柔性可视化着色） ---
    $sel set beta $rmsf
    animate write pdb "${opts(-out)}_rmsf.pdb" begiframe $f0 endiframe $f0 sel $sel
    puts "\[RMSF\] 已写出: ${opts(-out)}_rmsf.pdb (B-factor列 = RMSF值)"

    # --- 4d. 平均结构 PDB ---
    set avsel [atomselect $molid $opts(-sel)]
    $avsel set {x y z} $avpos
    $avsel frame $f0
    animate write pdb "${opts(-out)}_average.pdb" begiframe $f0 endiframe $f0 sel $avsel
    $avsel delete
    puts "\[RMSF\] 已写出: ${opts(-out)}_average.pdb (平均结构坐标)"

    # 提示可视化方法
    puts "\[RMSF\] 提示: 载入 ${opts(-out)}_rmsf.pdb 后,"
    puts "          Graphics->Representations->Coloring Method->Beta 即按RMSF着色"

    $sel delete
    return $rmsf
}

# ------------------------------------------------------------------
# 一键运行（载入体系 + 计算 + 输出）
#   ::rmsf::run structure.psf traj.dcd
#   ::rmsf::run complex.pdb {traj1.dcd traj2.dcd} -sel "protein and name CA"
# 额外参数透传给 ::rmsf::analyze
# ------------------------------------------------------------------
proc ::rmsf::run {topology traj args} {
    set molid [::rmsf::load_system $topology $traj]
    # 保证新载入的分子为 top
    molinfo $molid set top 1
    return [eval ::rmsf::analyze -molid $molid $args]
}

# ------------------------------------------------------------------
# 批量计算多体系（可选）：对每个 molid 分别输出
# ------------------------------------------------------------------
proc ::rmsf::analyze_all {args} {
    set result {}
    foreach molid [molinfo list] {
        if {[molinfo $molid get numframes] > 1} {
            lappend result $molid [eval ::rmsf::analyze -molid $molid $args]
        }
    }
    return $result
}

puts "=================================================================="
puts " rmsf_analysis.tcl v$::rmsf::version 载入成功"
puts " 用法示例:"
puts "   ::rmsf::analyze -sel \"protein and name CA\" -out rmsf_result"
puts "   ::rmsf::run structure.psf traj.dcd -first 0 -last -1"
puts " 可用选项: -sel -fitsel -first -last -step -molid -out"
puts "=================================================================="
