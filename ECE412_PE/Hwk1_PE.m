clc
clear
close all

vp_hat = 110*sqrt(2)
fs = 60
Ts = 1/fs
Po = 13
Vo = 26
Io = Po/Vo
w = 2*pi*fs
vd = 1
vdce = 2*vd+Vo

%%%%% Q1
N_max = 28/110/sqrt(2)/0.9
N_min = 28/110/sqrt(2)
N = N_max

vs_hat = N*vp_hat

te = 1/w*sqrt(2*(1-28/vs_hat))
Rs = 4*te*(vs_hat-Vo-2*vd)/Ts/Io

C_min = Io*(Ts/2-2*te)/0.05/Vo

%%%% Q2
%%%%%% vs
Tss = 1e-4;
tt = 0:Tss:Ts;
vs = vs_hat*cos(w*tt)
vdce_plot = (2*vd + Vo)*ones(size(vs))

%%%%%% is
is_high = (vs-vdce)/Rs;
is_high(is_high<0) = 0; 
is_low = (vs+vdce)/Rs;
is_low(is_low>0) = 0; 

is = is_high + is_low
is_rect_high = is_high
is_rect_high(is_rect_high>0) = (vs_hat-vdce)/Rs
is_rect_low = is_low
is_rect_low(is_rect_low<0) = -(vs_hat-vdce)/Rs
is_rect = is_rect_high + is_rect_low

%%%%%% id1
id = abs(is);
id_rect = abs(is_rect);

id1 = is;
id1_rect = is_rect;
id1(id1<0) = 0;
id1_rect(id1_rect<0) = 0;

%%%%%% vd1
vd1 = 13+vs/2-26;
vd1(is>0)=1;
vd1(is<0)=-27;

%%%%%% Capacitor
ic_rect = id_rect - Io;

%%%%%% Load
io = 26*ones(size(tt));

%%%%%%%%% plotting
figure(1)
plot(tt,vs,'DisplayName','v_s')
title('v_s vs time')
ylabel('voltage [V]')
xlabel('time [s]')
xlim([0 Ts])
grid on;
hold on;
plot(tt,vdce_plot,'DisplayName','v_{dce}')
plot(tt,-vdce_plot,'DisplayName','-v_{dce}')
scatter(te,vdce,'black','DisplayName','t_e')
legend('Location','best')

figure(2)
plot(tt,is,'DisplayName','i_s')
title('i_s vs time')
ylabel('current [A]')
xlabel('time [s]')
xlim([0 Ts])
grid on;
hold on;
plot(tt,is_rect,'DisplayName','i{s.approximation}')
legend('Location','best')

figure(3)
yyaxis left;
plot(tt,id1,'DisplayName','i_{d1}')
hold on;
plot(tt,id1_rect,'DisplayName','i_{d1.approx}')
ylabel('current [A]')
title('diode vs time')

yyaxis right;
ylabel('voltage [V]')
plot(tt,vd1,'DisplayName','v_{d1}')
xlabel('time [s]')
xlim([0 Ts])
grid on;
legend('Location','best')

figure(4)
plot(tt,ic_rect,'DisplayName','i_{C}')
yyaxis left;
title('Capacitor vs time')
ylabel('Current [A]')
hold on;
xlabel('time [s]')
xlim([0 Ts])
yyaxis right;
ylabel('Voltage [V]')
ylim([20 30])
grid on;
legend('Location','best')

figure(5)
plot(tt,io,'DisplayName','i_{o}')
yyaxis left;
title('Load vs time')
ylabel('Current [A]')
hold on;
xlabel('time [s]')
xlim([0 Ts])
yyaxis right;
ylabel('Voltage [V]')
ylim([20 30])
grid on;
legend('Location','best')
