set b 18.0000; 
set h 18.0000; 
set Acol [expr $b*$h]; 
set c 1.0000; 
set fc -5.1600; 
set fyl 106.4000; 
set TY 1.1600; 
set ful [expr $TY*$fyl]; 
set Es 26000.0000; 
set esh 0.0160; 
set eult 0.0860; 
set steel_b [expr $fyl*($TY-1.0)/$Es/($eult-$fyl/$Es)]; 
set Epyr 0.0250; 
set Esh [expr $Epyr*$Es]; 
set fyt 84.6000; 
set nlt 3; 
set nlb 3; 
set nlm 1; 
set nl [expr $nlt+$nlb+2*$nlm]; 
set nt 3; 
set s 3.5000; 
set Ashi 0.2000; 
set Asli 0.4400; 
set Ash [expr $nt*$Ashi]; 
set rou [expr $Ash/$b/$s]; 
set SecTag 23; 
set db 7.480315e-01; 
set ub [expr 12*sqrt(-$fc*1000)/1000]; 
set theta_sy [expr $fyl/$Es*$fyl*$db/8.0/$ub/($h-2*$c)]; 
set theta_su [expr $db/8.0/$ub/($h-2*$c)*($fyl/$Es*$fyl+2*($fyl/$Es+$eult)*($ful-$fyl))]; 
set betabuckle 0.9929; 
set modenum 1.0000; 
set MPR1 3.0000; 
set MPR2 34.6897; 
set MPR3 6.3330; 
set c_mono 1.6621; 
set c_cycl 2.1471; 
set c_symm 1.0500; 
set k1 10.5288; 
set k2 5.9905; 
set b1 0.0892; 
set b2 0.0547; 
set bs 0.3000; 
set BR 1.0000; 
set Sy 0.0472; 
set Su 0.3800; 
set DFTag 1; 
