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

Np = 10;
Nc = 10;
Nu = 4;
Nx = 24;
NOD = 16;

ACADOvariables.x = zeros(Np + 1, Nx);
ACADOvariables.u = zeros(Nc, Nu);
ACADOvariables.od = zeros(Np + 1, NOD);
ACADOvariables.y = zeros(Np, 7);
ACADOvariables.yN = zeros(1, 3);
ACADOvariables.x0 = zeros(1, Nx);

%% ³õÖµ
for i = 1 : Nc
%    ACADOvariables.u(i, 1) = 800;
%     ACADOvariables.u(i, 2) = 800;
%     ACADOvariables.u(i, 3) = 800;
    ACADOvariables.u(i, 4) = 800;
end

x0 = [0; 10; 10];
v0 = [1; 0; 0];
q0 = [0; 0; -1];
omega0 = [0; 0; 0];
R0 = [1; 0; 0; 0; 1; 0; 0; 0; 1];
OMEGA0 = [0; 0; 0];
tic
x_start = [x0; v0; q0; omega0; R0; OMEGA0; 0; 0; 0; (mq + ml) * g];

ACADOvariables.x0 = [x0; v0; q0; omega0; R0; OMEGA0]';
for i = 1:Np + 1
    ACADOvariables.od(i, 1) = x0(1);
    ACADOvariables.od(i, 2) = x0(2);
    ACADOvariables.od(i, 3) = x0(3) + 10;
    ACADOvariables.od(i, 4) = g;
    ACADOvariables.od(i, 5) = mq;
    ACADOvariables.od(i, 6) = ml;
    ACADOvariables.od(i, 7) = L;
    ACADOvariables.od(i, 8) = j11;
    ACADOvariables.od(i, 9) = j12;
    ACADOvariables.od(i, 10) = j13;
    ACADOvariables.od(i, 11) = j21;
    ACADOvariables.od(i, 12) = j22;
    ACADOvariables.od(i, 13) = j23;
    ACADOvariables.od(i, 14) = j31;
    ACADOvariables.od(i, 15) = j32;
    ACADOvariables.od(i, 16) = j33;
end
tic
data = acado_solver(ACADOvariables);
toc

T_val = 1;
[t, tmpx] = ode45(@QP_Model, [0 T_val], x_start);
toc
