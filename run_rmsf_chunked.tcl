# ======================================================================
# run_rmsf_chunked.tcl —— 分块增量 RMSF 计算
# 适用于 32 位 VMD 无法一次性载入大轨迹（内存 >2GB 崩溃）的场景
# 数学等价性: RMSF_i = sqrt( E[|r_i|^2] - |E[r_i]|^2 )
#   单遍累积每原子的 Σr 与 Σ|r|^2，分块处理，内存峰值 = CHUNK 帧
# ======================================================================
cd {D:/keepworking/2026/beta-galactose/0827}

set TOPO   complex00000.pdb
set TRAJ   complex.xtc
set SELTXT "protein and name CA"
set FITTXT "protein and backbone"
set CHUNK  100
set OUTP   complex_rmsf

set sums    {}   ;# 每原子 {sx sy sz}
set sumsq   {}   ;# 每原子平方和
set natoms  0
set nframes 0
set a       0

while {1} {
    set m [mol new $TOPO type pdb waitfor all]
    set last [expr {$a + $CHUNK - 1}]
    if {[catch {
        mol addfile $TRAJ type xtc first $a last $last step 1 waitfor all $m
    } adderr]} {
        puts "CHUNK: 读取帧 $a-$last 出错: $adderr (视为轨迹结束)"
        mol delete $m
        break
    }
    set nf [molinfo $m get numframes]
    set added [expr {$nf - 1}]
    if {$added <= 0} {
        mol delete $m
        break
    }

    # 首块时初始化累积器
    if {$natoms == 0} {
        set s0 [atomselect $m $SELTXT]
        set natoms [$s0 num]
        $s0 delete
        if {$natoms == 0} { error "选择 '$SELTXT' 未选中原子" }
        for {set i 0} {$i < $natoms} {incr i} {
            lappend sums {0.0 0.0 0.0}
            lappend sumsq 0.0
        }
        puts "CHUNK: CA原子数 = $natoms"
    }

    # ---- 处理本块每一帧: 对齐 + 累积 ----
    set ref [atomselect $m $FITTXT frame 0]
    set sel [atomselect $m $SELTXT]
    for {set i 1} {$i < $nf} {incr i} {
        set cur   [atomselect $m $FITTXT frame $i]
        set mat   [measure fit $cur $ref]
        set allat [atomselect $m all frame $i]
        $allat move $mat
        $cur delete
        $allat delete

        $sel frame $i
        set coords [$sel get {x y z}]
        for {set j 0} {$j < $natoms} {incr j} {
            set c [lindex $coords $j]
            set s [lindex $sums $j]
            set cx [lindex $c 0]; set cy [lindex $c 1]; set cz [lindex $c 2]
            lset sums  $j [list [expr {[lindex $s 0] + $cx}] \
                                [expr {[lindex $s 1] + $cy}] \
                                [expr {[lindex $s 2] + $cz}]]
            lset sumsq $j [expr {[lindex $sumsq $j] + $cx*$cx + $cy*$cy + $cz*$cz}]
        }
        incr nframes
    }
    $ref delete
    $sel delete
    mol delete $m

    set a [expr {$a + $added}]
    puts "CHUNK: 已处理 $nframes 帧 (轨迹进度至第 [expr {$a - 1}] 帧)"
}

if {$nframes < 2} { error "有效帧数不足: $nframes" }
puts "TOTAL_FRAMES = $nframes"

# ---- 由累积量计算平均坐标与 RMSF ----
set rmsf {}
set avgpos {}
for {set j 0} {$j < $natoms} {incr j} {
    set s [lindex $sums $j]
    set mx [expr {[lindex $s 0] / double($nframes)}]
    set my [expr {[lindex $s 1] / double($nframes)}]
    set mz [expr {[lindex $s 2] / double($nframes)}]
    set var [expr {[lindex $sumsq $j] / double($nframes) - ($mx*$mx + $my*$my + $mz*$mz)}]
    if {$var < 0.0} { set var 0.0 }
    lappend rmsf   [expr {sqrt($var)}]
    lappend avgpos [list $mx $my $mz]
}

set rmin 1e30; set rmax -1e30; set rsum 0.0
foreach r $rmsf {
    if {$r < $rmin} { set rmin $r }
    if {$r > $rmax} { set rmax $r }
    set rsum [expr {$rsum + $r}]
}
set ravg [expr {$rsum / double($natoms)}]
puts [format "RMSF_STATS min=%.4f max=%.4f mean=%.4f" $rmin $rmax $ravg]

# ---- 输出文件 ----
set m [mol new $TOPO type pdb waitfor all]
set sel [atomselect $m $SELTXT frame 0]
set resids   [$sel get resid]
set resnames [$sel get resname]
set chains   [$sel get chain]
set segids   [$sel get segid]
set names    [$sel get name]

# 逐原子数据
set fp [open "${OUTP}.atom.dat" w]
puts $fp "# Per-atom RMSF (Angstrom), chunked incremental computation"
puts $fp "# System: complex | frames: $nframes | sel: $SELTXT"
puts $fp "# cols: resid resname chain segid atomname rmsf"
for {set i 0} {$i < $natoms} {incr i} {
    puts $fp [format "%6d %4s %2s %6s %4s %10.4f" \
        [lindex $resids $i] [lindex $resnames $i] [lindex $chains $i] \
        [lindex $segids $i] [lindex $names $i] [lindex $rmsf $i]]
}
close $fp
puts "WROTE ${OUTP}.atom.dat"

# 逐残基聚合
array unset racc; array set racc {}
set keyorder {}
for {set i 0} {$i < $natoms} {incr i} {
    set key "[lindex $segids $i]:[lindex $chains $i]:[lindex $resids $i]"
    if {![info exists racc($key)]} {
        set racc($key) [list 0.0 0]
        lappend keyorder $key
    }
    set s $racc($key)
    lset racc($key) 0 [expr {[lindex $s 0] + [lindex $rmsf $i]}]
    lset racc($key) 1 [expr {[lindex $s 1] + 1}]
}
set fp [open "${OUTP}.residue.dat" w]
puts $fp "# Per-residue RMSF (Angstrom), chunked incremental computation"
puts $fp "# System: complex | frames: $nframes | sel: $SELTXT"
puts $fp "# cols: resid resname chain segid rmsf"
foreach key $keyorder {
    set s $racc($key)
    set mean [expr {[lindex $s 0] / double([lindex $s 1])}]
    set parts [split $key ":"]
    set rid [lindex $parts 2]
    set idx [lsearch -exact $resids $rid]
    set rname [lindex $resnames $idx]
    puts $fp [format "%6d %4s %2s %6s %10.4f" $rid $rname [lindex $parts 1] [lindex $parts 0] $mean]
}
close $fp
puts "WROTE ${OUTP}.residue.dat"

# B-factor PDB
$sel set beta $rmsf
animate write pdb "${OUTP}_rmsf.pdb" begiframe 0 endiframe 0 sel $sel
puts "WROTE ${OUTP}_rmsf.pdb"

# 平均结构 PDB
set avsel [atomselect $m $SELTXT frame 0]
$avsel set {x y z} $avgpos
animate write pdb "${OUTP}_average.pdb" begiframe 0 endiframe 0 sel $avsel
$avsel delete
puts "WROTE ${OUTP}_average.pdb"

mol delete $m
puts "ALL_DONE_OK"
exit 0
