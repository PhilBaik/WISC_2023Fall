clc
clear
close all
J = 0.2
p = 8
V_R = 72
I_R = 85
N_R = 3700
wm_R = 3700/60*2*pi;
T_R = 14;

k_T = T_R/I_R %%%%% DC Machine Constant
% k_E = V_R/wm_R 

lambda_f_RMS = k_T*sqrt(2)/(p/2)/pi;

P_R = 5000

R_ph = 12.4e-3;
L_ph = 154e-6;


NN = 0.1:0.1:3750;
wm = NN/60*2*pi;
we = wm*p/2;

X_ph = wm*L_ph;

Pm = zeros(size(wm));
for i = 1:1:length(NN)
    if NN(1,i)<3000
        Pm(1,i) = (NN(1,i)/3000)^3*5000;
    else
        Pm(1,i) = 5000;
    end
end

Pe_ph = -Pm;
Te_ph = Pe_ph./wm;

I_AC_rms = Te_ph*(2/(3*p*lambda_f_RMS));
Z_ph = (R_ph + j*X_ph);

E_rms = we*lambda_f_RMS;

Sout = 3*E_rms.*I_AC_rms;

%%%% E_rms I_rms
V_in = E_rms - I_AC_rms.*Z_ph;
V_in_rms = abs(V_in);

%%%%
I2R_ph = 3*I_AC_rms.^2*R_ph;
Pin = 3*V_in_rms.*I_AC_rms;
I_DC = Pin/V_R;


%%
close all
figure(11)
plot(NN,Pm)
title('Power output')

figure(12)
plot(NN,Te_ph)
title('Torque output')

figure(13)
plot(NN,I_AC_rms)
title('I AC rms')

figure(14)
plot(NN,E_rms)
title('E rms')

figure(15)
plot(NN,Sout)
title('Sout')

figure(16)
plot(NN,V_in_rms)
title('Vin_{rms}')

figure(17)
plot(NN,I_DC)
title('I_{DC}')

%%

