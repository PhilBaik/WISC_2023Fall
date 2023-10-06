%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Power Electronics ECE 412 Homework #2
%%%%%%% Writer: Hyeongmeen Baik
%%%%%%% Email: hyeongmeen.baik@wisc.edu
%%%%%%% Link: github.com/PhilBaik/WISC_2023Fall/blob/master/ECE412_PE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;

VF= 1
Vin = 100
D = 0.75
Vc = D*Vin - (1-D)*VF
Io = 10

tcr = 20e-9;
tvf = 30e-9;
tvr = 35e-9;
tcf = 15e-9;

Woff = 1/2*(Vin+VF)*Io*(tcf+tvr)
Won = 1/2*(Vin+VF)*Io*(tcr+tvf)
Rds_on = 0.18

Is = D*Io
Is_RMS = sqrt(D)*Io
is_peak = Io
Ps_cond = Rds_on*D*Io^2

Id = Io-D*Io
Id_RMS = sqrt(1-D)*Io
id_peak = Io
Pd_cond = VF*Io*(1-D)
