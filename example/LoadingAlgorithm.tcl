#
# RunStaticLoading: this procedure will conduct static loading simulation
# Author: Kuanshi Zhong
# Email: kuanshi@stanford.edu
# Date: 11/2016
#
# Arguments:
# IDctrlNode: the loaded node tag
# IDctrlDOF: the loaded DOF tag
# LoadType: the loading protocol type (Mono,CyclicStep,CyclicPeak)
# LoadHistory: the loading history vector, for Mono and CyclicPeak, only the target amplitude
#Dincr: the spacing of incremental loading steps
# numIter: the num of iteration reducing the current step
# Tol: the tolerance used in convergence tests

proc RunStaticLoading {IDctrlNode IDctrlDOF LoadType LoadHistory Dincr numIter Tol} {

	# generate the loading profile
	if {$LoadType == "CyclicStep"} {
		set LoadProfile $LoadHistory;
	}
	if {$LoadType == "CyclicPeak"} {
		set currentD 0.0;
		set outFileID [open tmpDsteps.tcl w];
		puts $outFileID "set LoadProfile { ";
		puts $outFileID $currentD;
		foreach Dmax $LoadHistory {
			if {$Dmax < 0} {
				set dx [expr -$Dincr];
			} else {
				set dx [expr $Dincr];
			}
			set NstepsPeak [expr int(abs($Dmax/$Dincr))];
			for {set i 1} {$i <= $NstepsPeak} {incr i 1} {
				set currentD [expr $currentD+$dx];
				puts $outFileID $currentD;
			}
			for {set i 1} {$i <= $NstepsPeak} {incr i 1} {
				set currentD [expr $currentD-$dx];
				puts $outFileID $currentD;
			}
		}
		puts $outFileID " }";
		close $outFileID;
		source tmpDsteps.tcl
	}
	if {$LoadType == "Mono"} {
		set currentD 0.0;
		set outFileID [open tmpDsteps.tcl w];
		puts $outFileID "set LoadProfile { ";
		puts $outFileID $currentD;
		foreach Dmax $LoadHistory {
			if {$Dmax < 0} {
				set dx [expr -$Dincr];
			} else {
				set dx [expr $Dincr];
			}
			set NstepsPeak [expr int(abs($Dmax/$Dincr))];
			for {set i 1} {$i <= $NstepsPeak} {incr i 1} {
				set currentD [expr $currentD+$dx];
				puts $outFileID $currentD;
			}
		}
		puts $outFileID " };";
		close $outFileID;
		source tmpDsteps.tcl
	}
#	puts "Load Protocol:"
	

	# analyze
	set F_ref 1.0;
	pattern Plain 200 Linear {
		load $IDctrlNode $F_ref 0.0 0.0;
	}
	constraints Plain
	numberer RCM
	system BandGeneral
	test EnergyIncr $Tol $numIter 0;
	algorithm Newton
	analysis Static
	set currentDisp 0.0;
	set currentIter 1;
	foreach Disp $LoadProfile {
#		puts $Disp;
		set nextDisp $Disp;
		set incrD [expr $nextDisp-$currentDisp];
		integrator DisplacementControl $IDctrlNode $IDctrlDOF $incrD;
		# conduct one-step analysis
		set ok [analyze 1];
		# if the result doesn't converge
		if {$ok != 0} {
			if {$ok != 0} {
				puts "Trying Newton with Initial Tangent ..";
				test EnergyIncr [expr $Tol*1] [expr 10*$numIter] 0;
				algorithm Newton -initial;
				set ok [analyze 1];
				test EnergyIncr $Tol $numIter 0;
				algorithm Newton;
			}
			if {$ok != 0} {
				puts "Trying Modified Newton ..";
				test EnergyIncr [expr $Tol*10] [expr 10*$numIter] 0;
				algorithm ModifiedNewton;
				set ok [analyze 1];
				test EnergyIncr $Tol $numIter 0;
				algorithm Newton
			}
			if {$ok != 0} {
				puts "Trying Modified Newton with Initial Tangent ..";
				test EnergyIncr [expr $Tol*10] [expr 10*$numIter] 0;
				algorithm ModifiedNewton -initial;
				set ok [analyze 1];
				test EnergyIncr $Tol $numIter 0;
				algorithm Newton
			}
			if {$ok != 0} {
				puts "Trying Broyden ..";
				test EnergyIncr [expr $Tol*10.0] [expr 10*$numIter] 0;
				algorithm Broyden 20;
				set ok [analyze 1];
				test EnergyIncr $Tol $numIter 0;
				algorithm Newton
			}
			if {$ok != 0} {
				puts "Trying NewtonWithLineSearch ..";
				test EnergyIncr [expr $Tol*10.0] [expr 10*$numIter] 0;
				algorithm NewtonLineSearch 0.8 100;
				set ok [analyze 1];
				test EnergyIncr $Tol $numIter 0;
				algorithm Newton
			}
			if {$ok != 0} {
#				set putout [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT];
#				puts $putout;
				return -1;
			};
		}; # end if
		# update the current displacement
		set currentDisp $nextDisp;
	}
	if {$ok != 0 } {
		puts "PROBLEM"
	} else {
		puts "DONE"
	}
}
		