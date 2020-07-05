# Declare a 2d model
model BasicBuilder -ndm 2

# Define constants
set pi [expr {2.0*asin(1.0)}]
set g 386.089

# Refer modeling subfiles or procedures
source DesignVariableRectangular.tcl
source CreateConcreteMaterial.tcl
source CreateRCColumnSection.tcl
source GetGaussLobattoIP.tcl
source LoadingAlgorithm.tcl
source LoadingParameterRectangular.tcl

# Define nodes
node 1	0 0
node 2	0 $Lcol
node 10001	0 0

# Define boundary condition
fix 10001	1 1 1
fix 1		0 0 0
fix 2		0 0 0

# Define shear section type
set ShearTag "NonLinear"

# Set up integration points
set LIP ""
set LIPR ""
set XIP ""
set IntegrationTag "GaussLobattol"
if {$IntegrationTag == "NewtonCotes"} {
	for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
		set LIP [lappend LIP [expr $Lcol/$numIntgrPts]];	# NewtonCotes: uniformly distributed integration points
	}
} else {
	set tempIP [GetGaussLobattolLIP $numIntgrPts];
	set tempXIP [lindex $tempIP 0]; 
	set tempLIP [lindex $tempIP 1];
	for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
		set LIPR [lappend LIPR [expr 0.5*[lindex $tempLIP $IPTag-1]]];
		set LIP [lappend LIP [expr 0.5*[lindex $tempLIP $IPTag-1]*$Lcol]];
		set XIP [lappend XIP [expr 0.5*[lindex $tempXIP $IPTag-1]+0.5]];
	}
}

# Define material properties
set nfCoreY 1
set nfCoreZ 80
set nfCoverY 1
set nfCoverZ 2
for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
	set ConcrCoverID	[expr $IPTag*10+1]
	set ConcrCoreID	[expr $IPTag*10+2]
	set SteelID		[expr $IPTag*10+3]
	DefineRegularizedUnconfinedConcreteMaterial "Concrete02" 1 $ConcrCoverID $fc [lindex $LIP $IPTag-1]
	DefineRegularizedConfinedConcreteMaterial "Concrete02" 1 $ConcrCoreID $fc $nl $s [expr $b-2.0*$c-2.0*$db] $rou $fyt [lindex $LIP $IPTag-1]
#	set Lgage 8
	set Lgage [lindex $LIP $IPTag-1]
	puts "Lgage/LIP = [expr $Lgage/[lindex $LIP $IPTag-1]]"
	set Epyr_reg [expr [lindex $LIP $IPTag-1]/$Lgage*0.01]
	set esu_temp [expr ($eult-$fyl/$Es-$esh)*$Lgage/[lindex $LIP $IPTag-1]+$fyl/$Es+$esh]
	
	if {$DFTag == 1} {
		uniaxialMaterial ReinforcingSteel [expr $SteelID+10000] $fyl $ful \
			$Es [expr [lindex $LIP $IPTag-1]/$Lgage*$Esh] $esh $esu_temp \
			-MPCurveParams [expr 1/$MPR1] $MPR2 $MPR3;
#		uniaxialMaterial Steel02 [expr $SteelID+10000] 160.0 29000.0 0.001 25 0.925 0.15;
		uniaxialMaterial DuctileFracture $SteelID [expr $SteelID+10000] \
			-c_mono $c_mono -c_cycl $c_cycl -c_symm $c_symm \
			-E_s $Es -esu $esu_temp -k1 $k1 -k2 $k2 -db $db -b1 $b1 -b2 $b2;
	} else {
#		uniaxialMaterial Steel02 $SteelID 160.0 29000.0 0.001 25 0.925 0.15;
		uniaxialMaterial ReinforcingSteel $SteelID $fyl $ful \
			$Es [expr [lindex $LIP $IPTag-1]/$Lgage*$Esh] $esh $esu_temp \
			-MPCurveParams [expr 1/$MPR1] $MPR2 $MPR3;
	}
	CreateColumnSection [expr 10+$IPTag] $h $b $c $ConcrCoreID $ConcrCoverID $SteelID $db $Asli $nlb $nlt $Asli $nlm $ShearTag $fc $fyt $rou $nfCoreY $nfCoreZ $nfCoverY $nfCoverZ
}

# Set up bar-slip section
set ke	[expr ($nl-2.)/$nl*(1-$s/[expr $b-2.0*$c-2.0*$db])]
set fl	[expr $ke*$rou*$fyt]
set Kfc 	[expr -1.254+2.254*sqrt(1.0+7.94*$fl/(-$fc))-2*$fl/(-$fc)]
set BarslipSteelID	10001
set BarslipConcrCoverID	10002
set BarslipConcrCoreID	10003
set BarslipAlpha		0.4
set Barslipb		$bs
set BarslipR		$BR
set Su_Sy			[expr 1.0+2.0*(1.0+$eult/($fyl/$Es))*(($ful/$fyl)-1)]
#set Su_Sy			40;
#set Sy			0.013
#set Sy			[expr 0.013+0.1*pow(((2*$BarslipAlpha+1)*$fyl*1000.0/sqrt(-$fc*1000.0)*$db/4000.0),(1.0/$BarslipAlpha))]
#set Sy [expr $Sy*1.6];
#set Sy			[expr (0.34+2.54*pow((($db*25.4)/8437.0*($fyl*6.895)/sqrt(-$fc*6.895)*(2*$BarslipAlpha+1)),(1.0/$BarslipAlpha)))/25.4]
#set Sy			[expr pow($fyl,2.0)*$db/(8.0*$Es*21.0*sqrt(-$fc*1000.0)/1000.0)]
puts "Sy = $Sy"
puts "Su/Sy = $Su_Sy"
uniaxialMaterial Bond_SP01 $BarslipSteelID $fyl $Sy $ful $Su $Barslipb $BarslipR
#uniaxialMaterial Bond_SP01 $BarslipSteelID $fyl 0.055 $ful 0.30 $Barslipb $BarslipR
DefineRegularizedUnconfinedConcreteMaterial "Concrete02" 1 $BarslipConcrCoverID $fc -1
DefineRegularizedConfinedConcreteMaterial "Concrete02" 1 $BarslipConcrCoreID $fc $nl $s $b $rou $fyt -1
CreateColumnSection 10001 $h $b $c $BarslipConcrCoreID $BarslipConcrCoverID $BarslipSteelID $db $Asli $nlb $nlt $Asli $nlm $ShearTag $fc $fyt $rou $nfCoreY $nfCoreZ $nfCoverY $nfCoverZ

# Create fiber elements
set transfTag 1
#geomTransf PDelta $transfTag
geomTransf Corotational $transfTag
#geomTransf Linear $transfTag
set secTags ""
for {set IPTag 1} {$IPTag <= $numIntgrPts} {incr IPTag} {
	set secTags [lappend secTags [expr 10+$IPTag]]
}
puts $secTags
set integration "UserDefined $numIntgrPts $secTags $XIP $LIPR"
element forceBeamColumn 1 1 2 $transfTag $integration
element zeroLengthSection 10001 10001 1 10001 -orient 0 1 0 -1 0 0

## Define shear spring element
#set Ec [expr 57.0*sqrt(-$fc*1000.0)]
#set G [expr $Ec/2/(1+0.25)]
#set k_shear [expr $G*$Acol/$Lcol*5.0/6.0]
#uniaxialMaterial Elastic 20001 $k_shear
#element zeroLength 20001 10001 1 -mat 20001 -dir 1

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
recorder Node -file ./CyclicOutputRectangular/disp.out -time -node 2 -dof 1 disp
recorder Node -file ./CyclicOutputRectangular/force.out -time -node 10001 -dof 1 reaction
recorder Element -file ./CyclicOutputRectangular/SteelBot.out -time -ele 1 section 1 fiber [expr -$h/2.0+$c+0.5*$db] [expr $b/2.0-$c-0.5*$db] stressStrain
recorder Element -file ./CyclicOutputRectangular/SteelTop.out -time -ele 1 section 1 fiber [expr $h/2.0-$c-0.5*$db] [expr $b/2.0-$c-0.5*$db] stressStrain
recorder Element -file ./CyclicOutputRectangular/SteelBot2.out -time -ele 1 section 2 fiber [expr -$h/2.0+$c+0.5*$db] [expr $b/2.0-$c-0.5*$db] stressStrain
recorder Element -file ./CyclicOutputRectangular/SteelTop2.out -time -ele 1 section 2 fiber [expr $h/2.0-$c-0.5*$db] [expr $b/2.0-$c-0.5*$db] stressStrain

recorder Element -file ./CyclicOutputRectangular/FISteelBot.out -time -ele 1 section 1 fiber [expr -$h/2.0+$c+0.5*$db] [expr $b/2.0-$c-0.5*$db] damage
recorder Element -file ./CyclicOutputRectangular/FISteelTop.out -time -ele 1 section 1 fiber [expr $h/2.0-$c-0.5*$db] [expr $b/2.0-$c-0.5*$db] damage
recorder Element -file ./CyclicOutputRectangular/FISteelBot2.out -time -ele 1 section 2 fiber [expr -$h/2.0+$c+0.5*$db] [expr $b/2.0-$c-0.5*$db] damage
recorder Element -file ./CyclicOutputRectangular/FISteelTop2.out -time -ele 1 section 2 fiber [expr $h/2.0-$c-0.5*$db] [expr $b/2.0-$c-0.5*$db] damage

recorder Element -file ./CyclicOutputRectangular/ConcrCoverBot.out -time -ele 1 section 1 fiber [expr -$h/2.0] 0.0 stressStrain
recorder Element -file ./CyclicOutputRectangular/ConcrCoverTop.out -time -ele 1 section 1 fiber [expr $h/2.0] 0.0 stressStrain
recorder Element -file ./CyclicOutputRectangular/ConcrCoreBot.out -time -ele 1 section 1 fiber [expr -$h/2.0+$c+$db] 0.0 stressStrain
recorder Element -file ./CyclicOutputRectangular/ConcrCoreTop.out -time -ele 1 section 1 fiber [expr $h/2.0-$c-$db] 0.0 stressStrain
recorder Element -file ./CyclicOutputRectangular/BarSlipForce.out -time -ele 10001 force
recorder Element -file ./CyclicOutputRectangular/BarSlipDisp.out -time -ele 10001 deformation
#recorder Element -file ./CyclicOutput/ShearForce.out -time -ele 20001 force
#recorder Element -file ./CyclicOutput/ShearDisp.out -time -ele 20001 deformation

# Conduct loading analysis
set LoadType "CyclicStep"
set Tol 1e-6
set numIter 800
RunStaticLoading 2 1 $LoadType $LoadHistory $Dincr $numIter $Tol