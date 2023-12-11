clc
clear

J = 0.2

V_R = 72
I_R = 85
N_R = 3700
wm_R = 3700/60*2*pi
T_R = 14

k_T = T_R/I_R
k_E = V_R/wm_R
p = 4
P_R = 5000

NN = 0.1:0.1:3750;
ww = NN/60*2*pi;
Pm = zeros(size(ww));
for i = 1:1:length(NN)
    if NN(1,i)<3000
        Pm(1,i) = (NN(1,i)/3000)^3*5000;
    else
        Pm(1,i) = 5000;
    end
end
Pe = -Pm;
Te = Pe./ww;
I_AC_rms = Pm/k_T;


%%
figure(11)
plot(NN,Pm)
figure(12)
plot(NN,Te)
figure(13)
plot(NN,I_AC_rms)

