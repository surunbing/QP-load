clc, clear
close all

global g;
global mq ml;
global L;
global j11 j12 j13 j21 j22 j23 j31 j32 j33;

g = 1.6;
mq = 200;
ml = 300;
L = 8;
j11 = 100;
j12 = 0;
j13 = 0;
j21 = 0;
j22 = 100;
j23 = 0;
j31 = 0;
j32 = 0;
j33 = 100;

%% ³õÖµ

x0 = [0; 10; 10];
v0 = [0; 0; 0];
q0 = [0; 0; -1];
omega0 = [0; 0; 0];
R0 = [1; 0; 0; 0; 1; 0; 0; 0; 1];
OMEGA0 = [0; 0; 0];
tic
x_start = [x0; v0; q0; omega0; R0; OMEGA0; 0; 0; 1; (mq + ml) * g + 10];

T_val = 1;
[t, tmpx] = ode45(@QP_Model, [0 T_val], x_start);
toc
