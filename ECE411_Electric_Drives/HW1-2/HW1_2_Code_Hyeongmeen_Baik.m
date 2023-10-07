%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Writer:  Hyeongmeen Baik 
%%%%%%%%%%%%%%%%%% Title:   HW1-2 Code
%%%%%%%%%%%%%%%%%% Date:    23-10-09
%%%%%%%%%%%%%%%%%% URL:     https://github.com/PhilBaik/WISC_2023Fall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;


%% Q1
clc;
clear;
close all;

hp = 30
n = 2500
%%%% Frame Size 288
Wf = 200;
Jm = 0.115;
Ra = 0.089;      %% Armature resistance 25C 
Ra_hot = Ra*1.2;
L = 0.00144;    %% Armature inducatnace
K = 0.850;      %% torque or voltage constant

%%%% a) 
%%%% Tau(w) = K*Va/Ra - K^2*w/Ra
%%%% w = (K*Va/Ra - Tau)*Ra/(K^2)

I1a_max = 25;
I1a_min = 0;
V1ain1 = 100;
V1ain2 = -100;
Tau1a_max = K*I1a_max;
Tau1a_min = K*I1a_min;
w1a_ur = (K*V1ain1/Ra_hot-Tau1a_min)*Ra_hot/(K^2)
w1a_lr = (K*V1ain1/Ra_hot-Tau1a_max)*Ra_hot/(K^2)
w1a_ul = (K*V1ain2/Ra_hot-Tau1a_min)*Ra_hot/(K^2)
w1a_ll = (K*V1ain2/Ra_hot-Tau1a_max)*Ra_hot/(K^2)

wnl1a1 = V1ain1/K
wnl1a2 = V1ain2/K
ww1a = 2*min(wnl1a1,wnl1a2):1:2*max(wnl1a1,wnl1a2);

Tau1a1 = K*V1ain1/Ra_hot-(K^2)*ww1a/Ra_hot;
Tau1a2 = K*V1ain2/Ra_hot-(K^2)*ww1a/Ra_hot;

Tau1a_max_plot = Tau1a_max*ones(size(ww1a));
Tau1a_min_plot = Tau1a_min*ones(size(ww1a));

Pmech1a_ul = w1a_ul*Tau1a_max;
Pmech1a_ur = w1a_ur*Tau1a_max;
Pmech1a_ll = w1a_ll*Tau1a_min;
Pmech1a_lr = w1a_lr*Tau1a_min;

figure(11)
plot(ww1a,Tau1a1,'DisplayName','Vin= 100 V')
title('Q1 a)')
xlabel('w [rad/s]')
ylabel('Tau [N*m]')
hold on;
plot(ww1a,Tau1a2,'DisplayName','Vin=-100 V')
xlim([min(w1a_lr,w1a_ll) max(w1a_ur,w1a_ul)])
plot(ww1a,Tau1a_max_plot,'DisplayName','Tau_{max}')
plot(ww1a,Tau1a_min_plot,'DisplayName','Tau_{min}')
ylim([Tau1a_min Tau1a_max])
grid on
legend('Location','best')

%%%% b) 
%%%% Tau(w) = K*Va/Ra - K^2*w/Ra
%%%% w = (K*Va/Ra - Tau)*Ra/(K^2)

I1b_max = 25;
I1b_min = -25;
V1bin1 = 100;
V1bin2 = -100;
Tau1b_max = K*I1b_max;
Tau1b_min = K*I1b_min;
w1b_ur = (K*V1bin1/Ra_hot-Tau1b_min)*Ra_hot/(K^2)
w1b_lr = (K*V1bin1/Ra_hot-Tau1b_max)*Ra_hot/(K^2)
w1b_ul = (K*V1bin2/Ra_hot-Tau1b_min)*Ra_hot/(K^2)
w1b_ll = (K*V1bin2/Ra_hot-Tau1b_max)*Ra_hot/(K^2)

wnl1b1 = V1bin1/K
wnl1b2 = V1bin2/K
ww1b = 2*min(wnl1b1,wnl1b2):1:2*max(wnl1b1,wnl1b2);

Tau1b1 = K*V1bin1/Ra_hot-(K^2)*ww1b/Ra_hot;
Tau1b2 = K*V1bin2/Ra_hot-(K^2)*ww1b/Ra_hot;

Tau1b_max_plot = Tau1b_max*ones(size(ww1b));
Tau1b_min_plot = Tau1b_min*ones(size(ww1b));

Pmech1b_ul = w1b_ul*Tau1b_max;
Pmech1b_ur = w1b_ur*Tau1b_max;
Pmech1b_ll = w1b_ll*Tau1b_min;
Pmech1b_lr = w1b_lr*Tau1b_min;

figure(12)
plot(ww1b,Tau1a1,'DisplayName','Vin= 100 V')
title('Q1 b)')
xlabel('w [rad/s]')
ylabel('Tau [N*m]')
hold on;
plot(ww1b,Tau1b2,'DisplayName','Vin=-100 V')
xlim([min(w1b_lr,w1b_ll) max(w1b_ur,w1b_ul)])
plot(ww1b,Tau1b_max_plot,'DisplayName','Tau_{max}')
plot(ww1b,Tau1b_min_plot,'DisplayName','Tau_{min}')
ylim([Tau1b_min Tau1b_max])
grid on
legend('Location','best')
%%
%%%% c)
%%%% Tau(w) = K*Va/Ra - K^2*w/Ra
Vin_max_1c = 100;
Vin_min_1c = -100;
nmax_1c = 3500;
nmin_1c = -3500;

wmax_1c = nmax_1c*0.10472;
wmin_1c = nmin_1c*0.10472;

KK = linspace(0,K,100000);

%%%% w = 3500rpm
Tau11_1c = KK*Vin_max_1c/Ra_hot -(KK.^2)*wmax_1c/Ra_hot;
Tau12_1c = KK*Vin_min_1c/Ra_hot -(KK.^2)*wmax_1c/Ra_hot;

Tau1max_1c = max(Tau11_1c);
Tau1min_1c = min(Tau12_1c);

Pout1_pos_1c = Tau1max_1c*wmax_1c
Pout1_neg_1c = Tau1min_1c*wmax_1c

%%%% w = -3500rpm
Tau21_1c = KK*Vin_max_1c/Ra_hot -(KK.^2)*wmin_1c/Ra_hot;
Tau22_1c = KK*Vin_min_1c/Ra_hot -(KK.^2)*wmin_1c/Ra_hot;

Tau2max_1c = max(Tau21_1c)
Tau2min_1c = min(Tau22_1c)

Pout2_pos_1c = Tau2max_1c*wmin_1c
Pout2_neg_1c = Tau2min_1c*wmin_1c

figure(13)
plot(KK,Tau11_1c,'DisplayName','Vin = 100, w = 3500rpm')
hold on;
plot(KK,Tau12_1c,'DisplayName','Vin = -100, w = 3500rpm')
grid on;
legend('Location','best')

figure(133)
plot(KK,Tau21_1c,'DisplayName','Vin = 100, w = -3500rpm')
hold on;
plot(KK,Tau22_1c,'DisplayName','Vin = -100, w = -3500rpm')
grid on;
legend('Location','best')
