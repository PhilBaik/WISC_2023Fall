%%%%%%%%%%%% ECE412 HW 5
%%%%%%%%%%%% Hyeongmeen Baik

clc
clear

Lm = 150e-6
d = 0.4
Fs = 75e3
Ts = 1/Fs
Vin = 100
Vo = 10
n = 15/100
Io = 20


syms M(D)
M(D) = D*n/(1-D)
subs(M(D),D,0.4)

del_im = Vin*d*Ts/Lm
del_is = del_im/n

is_min = Io/(1-d)-del_is/2
is_mid = is_min + del_is/2
is_max = is_min + del_is

vp_low  = -1*Vo/n

im_max = n*(is_min+del_is)
im_min = n*is_min
Im = (im_max+im_min)/2
Ip = Im*d
Im/n*(1-d)