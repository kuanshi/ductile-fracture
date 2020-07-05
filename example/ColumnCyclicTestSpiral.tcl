# Declare a 2d model
model BasicBuilder -ndm 2

# Define constants
set pi [expr {2.0*asin(1.0)}];
set g 386.089;
set	mm_in	25.4;
set	mpa_ksi	6.895;

# Refer modeling subfiles or procedures
source DesignVariableSpiral.tcl
source GetGaussLobattoIP.tcl
source LoadingAlgorithm.tcl
source LoadingParameterSpiral.tcl

# Define nodes
node 1	0 0
node 2	0 $L
node 10001	0 0

# Define boundary condition
fix 10001	1 1 1
fix 1		0 0 0
fix 2		0 0 0
#equalDOF 10001 1 1
set fixtag 10001;

# Define material tags
set	mattag_steel	1;
set	mattag_coverconcrete	2;
set	mattag_coreconcrete	3;
set	mattag_barslip1	4;
set	mattag_barslip2	5;
set	mattag_barslip3	6
set	mattag_shear	7;

# Define section tags
set	sectag_flexure	1;
set	sectag_shear	2;
set	sectag_barslip	3;
set	sectag_fiber	4;
set	sectag_zerolength	5;

# Define shear section type
set ShearTag "NonLinear"

# Set up integration points
set LIP ""
set LIPR ""
set XIP ""
set IntegrationTag "GaussLobattol"
if {$IntegrationTag == "NewtonCotes"} {
	for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
		set LIP [lappend LIP [expr $L/$numIntgrPts]];	# NewtonCotes: uniformly distributed integration points
	}
} else {
	set tempIP [GetGaussLobattolLIP $numIntgrPts];
	set tempXIP [lindex $tempIP 0]; 
	set tempLIP [lindex $tempIP 1];
	for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
		set LIPR [lappend LIPR [expr 0.5*[lindex $tempLIP $IPTag-1]]];
		set LIP [lappend LIP [expr 0.5*[lindex $tempLIP $IPTag-1]*$L]];
		set XIP [lappend XIP [expr 0.5*[lindex $tempXIP $IPTag-1]+0.5]];
	}
}

# Create concrete mateirals
# cover concrete
set	eco	[expr 2.0*$fc/$Ec];
# core concrete
set	ds	[expr $D-2.0*$c-$dbt];
set	Ac	[expr 0.25*$pi*$ds*$ds];
set	rouCC	[expr $nsl*$Asl/$Ac];
set	Acc	[expr $Ac*(1-$rouCC)];
set	Ae	[expr 0.25*$pi*($ds-0.5*($s-$dbt))*($ds-0.5*($s-$dbt))];
set	ke	[expr $Ae/$Acc];
set	rouS	[expr 4.0*$Ast/$ds/$s];
set	fl	[expr 0.5*$rouS*$fyt];
set	flp	[expr $fl*$ke];
if {$LWTag == 0} {
	set	fcc	[expr $fc*(-1.254+2.254*sqrt(1+7.94*$flp/abs($fc))-2.0*$flp/abs($fc))];
	set	ecc	[expr $eco*(1.0+5.0*($fcc/$fc-1.0))];
	set	Gfc	[expr 1.7*2.0*(-$fc)*$mpa_ksi];
} else {
	set	fcc	[expr $fc*(-1.254+2.254*sqrt(1+7.94*$flp/abs($fc))-2.0*$flp/abs($fc))];
	set	ecc	[expr 0.75*$eco*(1.0+5.0*($fcc/$fc-1.0))];
	set	Gfc	[expr 0.75*1.7*2.0*(-$fc)*$mpa_ksi];
}
set	fpcu	[expr 0.2*$fcc];
puts "fcc = $fcc";
puts "ecc = $ecc";
for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
	set	ConcrCoverID	[expr $IPTag*10+1];
	set	ConcrCoreID		[expr $IPTag*10+2];
	set	SteelID		[expr $IPTag*10+3];
	set	tempLIP		[lindex $LIP $IPTag-1];
	set	ecu	[expr -($Gfc/0.6/(-$fcc*$mpa_ksi)/($tempLIP*$mm_in)-0.8*(-$fcc)/$Ec+(-$ecc))];
	uniaxialMaterial Concrete02 $ConcrCoverID $fc $eco 0.0 -0.004 0.4 $ft $Et;
	uniaxialMaterial Concrete02 $ConcrCoreID $fcc $ecc $fpcu $ecu 0.4 $ft $Et;
	# Create steel material
	if {$DFTag == 1} {
		uniaxialMaterial ReinforcingSteel [expr 1000+$SteelID] $fyl $ful $Esl $Esh $esh $esu; \
#			-GABuck [expr $s/$dbl] 2.0 0.4 0.5;
#			-DMBuck [expr $s/$dbl] 1.0;
		uniaxialMaterial DuctileFracture $SteelID [expr 1000+$SteelID] \
			-c_mono $c_mono -c_cycl $c_cycl -c_symm $c_symm \
			-E_s $Esl -esu $esu -k1 $k1 -k2 $k2 \
			-db $dbl -b1 $b1 -b2 $b2;
	} else {
		uniaxialMaterial ReinforcingSteel $SteelID $fyl $ful $Esl $Esh $esh $esu; \
#			-GABuck [expr $s/$dbl] 2.0 0.4 0.5;
#			-DMBuck [expr $s/$dbl] 1.0;
	}
}
# Create bar-slip material
set	alpha	0.4;
#set	Sy	[expr 0.1*pow($dbl/4000.0*$fyl*1000.0/sqrt(-$fc*1000.0)*(2.0*$alpha+1.0),1.0/$alpha)+0.013];
#set	Su_Sy	[expr 1.0+2.0*(1.0+$esu/($fyl/$Esl))*(($ful/$fyl)-1)];
#set	Su_Sy 40;
#set	Su	[expr $Su_Sy*$Sy];
set	bs	$bs;
set	R	$BR;
puts "Sy = $Sy";
puts "Su = $Su";
puts "bs = $bs";
puts "R = $R";
uniaxialMaterial Bond_SP01 $mattag_barslip1 $fyl $Sy $ful $Su $bs $R;
set	ecu_barslip	[expr -($Gfc/0.6/(-$fcc*$mpa_ksi)/(1.0*$mm_in)-0.8*(-$fcc)/$Ec+(-$ecc))];
#uniaxialMaterial Concrete02 $mattag_barslip2 $fcc $ecc [expr 4.0*$fpcu] [expr 100000.0*$ecu] 0.1 $ft [expr 0.05*$Ec];
uniaxialMaterial Concrete02 $mattag_barslip2 $fcc $ecc [expr $fcc*$CR] [expr $ecu_barslip*$ER] 0.4 $ft $Et;
uniaxialMaterial Concrete02 $mattag_barslip3 $fc $eco 0.0 -0.004 0.4 $ft $Et;

# Define section fiber size
set	numCirc	90;
set	numRad1	3;
set	numRad2	60;
set	intRad	[expr 0.5*$D-$c-$dbt];
set	extRad	[expr 0.5*$D];
# Create cross sections per IP point
# shear
set	Gc	[expr $Ec/2.0/(1.0+0.2)];
set	vn	[expr 3.0*sqrt(-$fc*1000.0)/1000.0+$rouS*$fyt];
if {$vn > [expr 8.0*sqrt(-$fc*1000.0)/1000.0]} {
	set vn [expr 8.0*sqrt(-$fc*1000.0)/1000.0];
}
set	s1p	[expr 0.002*sqrt(-$fc*1000.0)*$Ac*0.7];
set	e1p	[expr $s1p/$Gc/$Ac/0.7];
set	s2p	[expr 0.6*$vn*$Ac*0.7];
set	e2p	[expr $e1p+($s2p-$s1p)/0.4/$Gc/$Ac/0.7];
set	s3p	[expr $vn*$Ac*0.7];
set	e3p	[expr $e2p+0.4*$vn/0.1/$Gc];
set	s1n	[expr -$s1p];
set	e1n	[expr -$e1p];
set	s2n	[expr -$s2p];
set	e2n	[expr -$e2p];
set	s3n	[expr -$s3p];
set	e3n	[expr -$e3p];
set	pinchX	0.9;
set	pinchY	0.1;
set	damage1	0.0;
set	damage2	0.0;
set	beta	0.4;
if {$ShearTag == "Linear"} {
	uniaxialMaterial Elastic $mattag_shear [expr 0.1*$s1p/$e1p];
} else {
	uniaxialMaterial Hysteretic $mattag_shear $s1p $e1p $s2p $e2p $s3p $e3p $s1n $e1n $s2n $e2n $s3n $e3n $pinchX $pinchY $damage1 $damage2 $beta;
}
# flexure
for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
	section Fiber [expr $IPTag*100+1] {
	# cover concrete
	patch circ [expr $IPTag*10+1] $numCirc $numRad1 0.0 0.0 $intRad $extRad 0.0 360.0;
	# core concrete
	patch circ [expr $IPTag*10+2] $numCirc $numRad2 0.0 0.0 0.0 $intRad 0.0 360.0;
	# reinforcement
#	layer circ [expr $IPTag*10+3] $nsl $Asl 0.0 0.0 $intRad 0.0 [expr 360.0*(1-1/$nsl)];
	set	bartag	0;
	while {$bartag<$nsl} {
		if {$barlayout == 1} {
			set	yLoc	[expr $intRad*cos($bartag*2.0*$pi/$nsl)];
			set	zLoc	[expr $intRad*sin($bartag*2.0*$pi/$nsl)];
		}
		if {$barlayout == 2} {
			set	yLoc	[expr $intRad*cos($bartag*2.0*$pi/$nsl+$pi/$nsl)];
			set	zLoc	[expr $intRad*sin($bartag*2.0*$pi/$nsl+$pi/$nsl)];
		}
		fiber $yLoc $zLoc $Asl [expr $IPTag*10+3];
		incr bartag;
	}
	}
	section Aggregator [expr $IPTag*10+1] $mattag_shear Vy -section [expr $IPTag*100+1];
}
# bar-slip
section Fiber $sectag_barslip {
# cover concrete
patch circ $mattag_barslip3 $numCirc $numRad1 0.0 0.0 $intRad $extRad 0.0 360.0;
# core concrete
patch circ $mattag_barslip2 $numCirc $numRad2 0.0 0.0 0.0 $intRad 0.0 360.0;
# reinforcement
set	bartag	0;
while {$bartag<$nsl} {
	if {$barlayout == 1} {
		set	yLoc	[expr $intRad*cos($bartag*2.0*$pi/$nsl)];
		set	zLoc	[expr $intRad*sin($bartag*2.0*$pi/$nsl)];
	}
	if {$barlayout == 2} {
		set	yLoc	[expr $intRad*cos($bartag*2.0*$pi/$nsl+$pi/$nsl)];
		set	zLoc	[expr $intRad*sin($bartag*2.0*$pi/$nsl+$pi/$nsl)];
	}
	fiber $yLoc $zLoc $Asl $mattag_barslip1;
	incr bartag;	
}
}
section Aggregator $sectag_zerolength $mattag_shear Vy -section $sectag_barslip;

# Create fiber elements
set transfTag 1
#geomTransf PDelta $transfTag
geomTransf Corotational $transfTag
#geomTransf Linear $transfTag
set secTags ""
for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
	set secTags [lappend secTags [expr 10*$IPTag+1]]
}
puts $secTags
set integration "UserDefined $numIntgrPts $secTags $XIP $LIPR"
element forceBeamColumn 1 1 2 $transfTag $integration
element zeroLengthSection 10001 10001 1 $sectag_zerolength -orient 0 1 0 -1 0 0

# Define axial loads
pattern Plain 1 Linear {
	load 2 0.0 $P 0.0
}

# Update state
puts "Complete modeling."

# Conduct the gravity analysis
set numsteps_grav 10
set tol 1e-8
set maxiter 500
constraints Plain
numberer RCM
system BandGeneral
test RelativeEnergyIncr $tol $maxiter
algorithm Newton
integrator LoadControl [expr {1.0/$numsteps_grav}]
analysis Static
if {[analyze $numsteps_grav]} {
	puts "Application of gravity load failed"
} else {
	puts "Applied gravity loads."
}
loadConst -time 0.0
wipeAnalysis

# Set up recorder
recorder Node -file ./CyclicOutputSpiral/disp.out -time -node 2 -dof 1 disp
recorder Node -file ./CyclicOutputSpiral/force.out -time -node $fixtag -dof 1 reaction
recorder Element -file ./CyclicOutputSpiral/SteelBot.out -time -ele 1 section 1 fiber [expr -$D/2.0+$c+0.5*$dbl] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/SteelTop.out -time -ele 1 section 1 fiber [expr $D/2.0-$c-0.5*$dbl] 0.0 stressStrain

set	bartag	0;
while {$bartag<$nsl} {
	if {$barlayout == 1} {
		recorder Element -file ./CyclicOutputSpiral/SS$bartag.out -time -ele 1 section 1 fiber [expr $intRad*cos($bartag*2.0*$pi/$nsl)] [expr $intRad*sin($bartag*2.0*$pi/$nsl)] stressStrain;
		recorder Element -file ./CyclicOutputSpiral/FI$bartag.out -time -ele 1 section 1 fiber [expr $intRad*cos($bartag*2.0*$pi/$nsl)] [expr $intRad*sin($bartag*2.0*$pi/$nsl)] damage;
	}
	if {$barlayout == 2} {
		recorder Element -file ./CyclicOutputSpiral/SS$bartag.out -time -ele 1 section 1 fiber [expr $intRad*cos($bartag*2.0*$pi/$nsl+$pi/$nsl)] [expr $intRad*sin($bartag*2.0*$pi/$nsl+$pi/$nsl)] stressStrain;
		recorder Element -file ./CyclicOutputSpiral/FI$bartag.out -time -ele 1 section 1 fiber [expr $intRad*cos($bartag*2.0*$pi/$nsl+$pi/$nsl)] [expr $intRad*sin($bartag*2.0*$pi/$nsl+$pi/$nsl)] damage;
	}
	incr bartag;
}

recorder Element -file ./CyclicOutputSpiral/ConcrCoverBot.out -time -ele 1 section 1 fiber [expr -$D/2.0] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/ConcrCoverTop.out -time -ele 1 section 1 fiber [expr $D/2.0] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/ConcrCoreBot.out -time -ele 1 section 1 fiber [expr -$D/2.0+$c+1.5*$dbl] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/ConcrCoreTop.out -time -ele 1 section 1 fiber [expr $D/2.0-$c-1.5*$dbl] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/BarSlipForce.out -time -ele 10001 force
recorder Element -file ./CyclicOutputSpiral/BarSlipDisp.out -time -ele 10001 deformation
recorder Element -file ./CyclicOutputSpiral/BarSlipConcrCoverBot.out -time -ele 10001 section fiber [expr -$D/2.0] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/BarSlipConcrCoverTop.out -time -ele 10001 section fiber [expr $D/2.0] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/BarSlipConcrCoreBot.out -time -ele 10001 section fiber [expr -0.5*$D+$c+2.0*$dbl] 0.0 stressStrain
recorder Element -file ./CyclicOutputSpiral/BarSlipConcrCoreTop.out -time -ele 10001 section fiber [expr 0.5*$D-$c-2.0*$dbl] 0.0 stressStrain

# Conduct loading analysis
set LoadType "CyclicStep"
set Tol 1e-6
set numIter 800
RunStaticLoading 2 1 $LoadType $LoadHistory $Dincr $numIter $Tol