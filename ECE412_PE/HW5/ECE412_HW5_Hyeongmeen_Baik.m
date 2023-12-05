%%%%%%%%%%%% ECE412 HW 5
%%%%%%%%%%%% Hyeongmeen Baik

hc = 30e-3
hw = 20e-3
d = 3.183e-3
w0 = 12e-3
w1 = 12e-3
w2 = 6e-3
ur = 1000
u0 = 4*pi*(10e-7);
A1 = w0*d;
A2 = (hc-hw)/2*d;
A3 = w2*d;
A4 = A2;
l1 = (hc+hw)/2
l2 = (w0+w1+w2+w1)/2
l3 = l1
l4 = l2


R1 = l1/u0/ur/A1
R2 = l2/u0/ur/A2
R3 = l3/u0/ur/A3
R4 = l4/u0/ur/A4

Rtot = R1+R2+R3+R4

i1 = 0.2
i2 = 0.15
N1 = 150
N2 = 300

flux = (i1*N1+i2*N2)/Rtot

fluxlinkage1 = flux*N1
fluxlinkage2 = flux*N2

B1= flux/A1
B3 = flux/A3


%% q3
f = 150
Lm = N1*N2/Rtot

V1_ph = 1
R = 5

Ip_ph = 4*1/R
Im_ph = V1_ph/(j*2*pi*f*Lm)

I1_ph = Ip_ph + Im_ph;
I1 = abs(I1_ph);