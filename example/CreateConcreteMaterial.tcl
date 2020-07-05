####################################################################################################
# DefineConcreteMaterial: Defines a procedure which create the user-defined concrete material
#
# Authoer: Kuanshi Zhong
# Date: 11/2016
#
# Formal arguments
# fc: the conpressive strength (negative, ksi)
# nl: the number of longitudinal bars restrained by corners of hoops or legs of crossties
# s: the lateral steel spacing
# b: the column width
# rou: the lateral reinforcement ratio Ash/(b*s) in the loading direction
# fyt: the yield strength of lateral steel
#
# Notes
# Wall
# All materials defined below are regularized following Pugh et al. (2015)

# Define an unconfined concrete material
proc DefineRegularizedUnconfinedConcreteMaterial {type UnitTag matTag fc LIP} {
# Define unit converters
if {$UnitTag == 1} {
	set mm_in 25.4;
	set mpa_ksi 6.895;
	set in_mm [expr 1.0/$mm_in];
	set ksi_mpa [expr 1.0/$mpa_ksi];
} else {
	set mm_in 1;
	set mpa_ksi 1;
	set in_mm [expr 1.0/25.4];
	set ksi_mpa [expr 1.0/6.895];
}
# Concrete 01 Material: Zero Tensile Strength
if {$type == "Concrete01"} {
	set fpc	[expr $fc*$mpa_ksi*$ksi_mpa];
	set Ec	[expr 57.0*sqrt(-$fc*1000.0)];
	set epsc0	[expr 2.0*$fc/$Ec];
	set fpcu	[expr 0.2*$fpc];
	set Gfc	[expr 2.0*(-$fc)*$mpa_ksi];
	set epsU	-0.008;
	uniaxialMaterial Concrete01 $matTag $fpc $epsc0 $fpcu $epsU;
}
# Concrete 02 Material: Linear Tension Softening
if {$type == "Concrete02"} {
	set fpc	[expr $fc*$mpa_ksi*$ksi_mpa];
	set Ec	[expr 57.0*sqrt(-$fc*1000.0)];
	set epsc0	[expr 2.0*$fc/$Ec];
	set fpcu	[expr 0.2*$fpc];
	set Gfc	[expr 2.0*(-$fc)*$mpa_ksi];
#	set epsU	[expr -($Gfc/0.6/(-$fc*$mpa_ksi)/(abs($LIP)*$mm_in)-0.8*(-$fc)/$Ec+(-$epsc0))];
	set epsU	-0.008;
	set lambda	0.1;
	set ft	[expr 0.004*sqrt(-$fc*1000.0)];
	set Ets	[expr $Ec*0.05];
	uniaxialMaterial Concrete02 $matTag $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets;
}
# Concrete04 Material: Popovics Concrete Material
if {$type == "Concrete04"} {
	set fpc	[expr $fc*$mpa_ksi*$ksi_mpa];
	set Ec	[expr 57*sqrt(-$fc*1000)];
	set epsc0	[expr 2.0*$fc/$Ec];
	set Gfc	[expr 2*(-$fc)*$mpa_ksi];
	set epsU	[expr -($Gfc/0.6/(-$fc*$mpa_ksi)/($LIP*$mm_in)-0.8*(-$fc)/$Ec+(-$epsc0))];
	set ft	[expr 0.004*sqrt(-$fc*1000)];
	set et	0.002;
	set beta	0.0; # the residual stress ratio
	uniaxialMaterial Concrete04 $matTag $fpc $epsc0 $epsU $Ec $ft $et $beta;
}
puts "Unconfined Concrete Material Defined"
puts "fpc is $fpc";
puts "epsc0 is $epsc0";
puts "epsU is $epsU";
}

# Define a confined concrete material
proc DefineRegularizedConfinedConcreteMaterial {type UnitTag matTag fc nl s b rou fyt LIP} {
# Define unit converters
if {$UnitTag == 1} {
	set mm_in 25.4;
	set mpa_ksi 6.895;
	set in_mm [expr 1.0/$mm_in];
	set ksi_mpa [expr 1.0/$mpa_ksi];
} else {
	set mm_in 1;
	set mpa_ksi 1;
	set in_mm [expr 1.0/25.4];
	set ksi_mpa [expr 1.0/6.895];
}
# Concrete 01 Material: Zero Tensile Strength
if {$type == "Concrete01"} {
	set ke	[expr ($nl-2.)/$nl*(1-$s/$b)];
	set fl	[expr $ke*$rou*$fyt];
	set Kfc 	[expr -1.254+2.254*sqrt(1.0+7.94*$fl/(-$fc))-2*$fl/(-$fc)];
	set fpc	[expr $Kfc*$fc];
	set Keps	[expr 1+5*($Kfc-1)];
	set epsc0	[expr 2.0*$fc/$Ec];
	set epsc0	[expr $Keps*$epsc0];
	set fpcu	[expr 0.2*$fpc];
	set Gfc	[expr 1.7*2*(-$fc)*$mpa_ksi];
	set epsU	[expr -($Gfc/0.6/(-$fpc*$mpa_ksi)/($LIP*$mm_in)-0.8*(-$fpc)/$Ec+(-$epsc0))];

#	set flx	[expr $roux*$fyt];
#	set fly	[expr $rouy*$fyt];
#	set k2x	[expr 0.26*sqrt($bx/$s*$bx/$slx/$flx)];
#	if {$k2x > 1.0} {
#		set k2x 1.0;
#	}
#	set k2y	[expr 0.26*sqrt($by/$s*$by/$sly/$fly)];
#	if {$k2y > 1.0} {
#		set k2y 1.0;
#	}
#	set flex	[expr $k2x*$flx];
#	set fley	[expr $k2y*$fly];
#	set fle	[expr ($flex*$bx+$fley*$by)/($bx+$by)];
#	set k1	[expr 6.7*pow($fle,-0.17)];
#	set fpc	[expr $fc-$k1*$fle];
#	set Ec	[expr 57.*sqrt(-$fc*1000.)];
#	set epsc0	[expr 2.0*$fc/$Ec*(1+5.0*$k1*$fle/(-$fc))];
#	set roum	[expr ($roux*$bx+$rouy*$by)/($bx+$by)];
#	set eps85	[expr 260.0*$roum*$epsc0+0.0038];
#	set kunload	[expr 0.15*(-$fpc)/($eps85-$epsc0)];
#	set eps20	[expr $epsc0+0.8*(-$fpc)/$kunload];
#	set Gfc	[expr 1.7*2*(-$fc)*$mpa_ksi];

	if {$LIP > 0} {
		set fpcu	[expr 0.2*$fpc];
		set epsU	[expr -($Gfc/0.6/(-$fpc*$mpa_ksi)/($LIP*$mm_in)-0.8*(-$fpc)/$Ec+(-$epsc0))];
	} else {
		# for bar-slip section
		set fpcu	[expr 0.8*$fpc];
		set epsU	[expr 5.0*$epsc0];
	}

	uniaxialMaterial Concrete01 $matTag $fpc $epsc0 $fpcu $epsU;
}
# Concrete 02 Material: Linear Tension Softening
if {$type == "Concrete02"} {
	set Ec	[expr 57*sqrt(-$fc*1000)];
	set ke	[expr ($nl-2.)/$nl*(1-$s/$b)];
	set fl	[expr $ke*$rou*$fyt];
	set Kfc 	[expr -1.254+2.254*sqrt(1+7.94*$fl/(-$fc))-2*$fl/(-$fc)];
	set fpc	[expr $Kfc*$fc];
	set Keps	[expr 1+5*($Kfc-1)];
	set epsc0	[expr 2.0*$fc/$Ec];
	set epsc0	[expr $Keps*$epsc0];
	set fpcu	[expr 0.2*$fpc];
	set Gfc	[expr 1.7*2*(-$fc)*$mpa_ksi];
	set epsU	[expr -($Gfc/0.6/(-$fpc*$mpa_ksi)/($LIP*$mm_in)-0.8*(-$fpc)/$Ec+(-$epsc0))];
	set lambda	0.1;

#	set flx	[expr $roux*$fyt];
#	set fly	[expr $rouy*$fyt];
#	set k2x	[expr 0.26*sqrt($bx/$s*$bx/$slx/$flx)];
#	if {$k2x > 1.0} {
#		set k2x 1.0;
#	}
#	set k2y	[expr 0.26*sqrt($by/$s*$by/$sly/$fly)];
#	if {$k2y > 1.0} {
#		set k2y 1.0;
#	}
#	set flex	[expr $k2x*$flx];
#	set fley	[expr $k2y*$fly];
#	set fle	[expr ($flex*$bx+$fley*$by)/($bx+$by)];
#	set k1	[expr 6.7*pow($fle,-0.17)];
#	set fpc	[expr $fc-$k1*$fle];
#	set Ec	[expr 57.*sqrt(-$fc*1000.)];
#	set epsc0	[expr 2.0*$fc/$Ec*(1+5.0*$k1*$fle/(-$fc))];
#	set roum	[expr ($roux*$bx+$rouy*$by)/($bx+$by)];
#	set eps85	[expr 260.0*$roum*$epsc0+0.0038];
#	set kunload	[expr 0.15*(-$fpc)/($eps85-$epsc0)];
#	set eps20	[expr $epsc0+0.8*(-$fpc)/$kunload];
#	set Gfc	[expr 1.7*2*(-$fc)*$mpa_ksi];

	if {$LIP > 0} {
		set fpcu	[expr 0.2*$fpc];
		set epsU	[expr -($Gfc/0.6/(-$fpc*$mpa_ksi)/($LIP*$mm_in)-0.8*(-$fpc)/$Ec+(-$epsc0))];
	} else {
		set fpcu	[expr 0.8*$fpc];
#		set epsU	[expr $Kfc*10*$epsc0];
#		set epsU	[expr 5.0*$epsc0];
#		set fpcu	[expr 0.2*$fpc];
		set epsU	[expr -($Gfc/0.6/(-$fpc*$mpa_ksi)/(1.0*$mm_in)-0.8*(-$fpc)/$Ec+(-$epsc0))];
		puts "Barslip section concrete: epsU = $epsU";
	}
	set lambda	0.1;
	set ft	[expr 0.004*sqrt(-$fc*1000)];
	set Ets	[expr $Ec*0.05];
	uniaxialMaterial Concrete02 $matTag $fpc $epsc0 $fpcu $epsU $lambda $ft $Ets;
}
# Concrete04 Material: Popovics Concrete Material
if {$type == "Concrete04"} {
	set ke	[expr ($nl-2.)/$nl*(1-$s/$bx)];
	set fl	[expr $ke*$roux*$fyt];
	set Kfc 	[expr -1.254+2.254*sqrt(1+7.94*$fl/(-$fc))-2*$fl/(-$fc)];
	set fpc	[expr $Kfc*$fc];
	set Keps	[expr 1+5*($Kfc-1)];
	set Ec	[expr 57*sqrt(-$fc*1000)];
	set epsc0	[expr 2.0*$fc/$Ec];
	set epsc0	[expr $Keps*$epsc0];
	set Gfc	[expr 1.7*2*(-$fc)*$mpa_ksi];
	set epsU	[expr -($Gfc/0.6/(-$fpc*$mpa_ksi)/($LIP*$mm_in)-0.8*(-$fpc)/$Ec+(-$epsc0))];
	set ft	[expr 0.004*sqrt(-$fc*1000)];
	set et	0.002;
	set beta	0.0; # the residual tensile stress ratio
	uniaxialMaterial Concrete04 $matTag $fpc $epsc0 $epsU $Ec $ft $et $beta;
}
puts "Confined Concrete Material Defined"
puts "nl is $nl"
puts "s is $s"
puts "b is $b"
puts "ke is $ke"
puts "fl is $fl"
puts "Kfc is $Kfc"
#puts "Keps is $Keps"
#puts "flex is $flex";
#puts "fley is $fley";
puts "fpc is $fpc";
puts "epsc0 is $epsc0";
puts "epsU is $epsU";
}