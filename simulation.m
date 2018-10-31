clc, clear
close all

global mF ml l AL cDL g rho p q r fT fHx fHy fHz

mF = 0.93;
ml = 0.2;
l = 1.0;
AL = 0.0142;
cDL = 1.2;
g = 0.9810;
rho = 1.1840;

p = 0;
q = 0;
r = 0.001;
fT = (mF + ml) * g;
fHx = 0;
fHy = 0;
fHz = ml * g;

tic
x0 = [1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0];
T_val = 1;
[t, tmpx] = ode45(@QP_Model, [0 T_val], x0);
toc
