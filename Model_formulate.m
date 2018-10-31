clc,clear
close all

syms phi theta psi

Rx = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
Ry = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
Rz = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];

R = (Rx * Ry * Rz)';

syms p q r vx vy vz
omega = [p q r];
v = [vx vy vz];
ans = -cross(omega, v);

syms g 
RT = R';
ans = g * RT * [0; 0; 1];

syms mF fT fHx fHy fHz
fH = [fHz; fHy; fHz];
ans = (RT * fH - fT * [0; 0; 1]) / mF;

J = [cos(theta), sin(phi) * sin(theta), cos(phi) * sin(theta); 0, cos(phi) * cos(theta), -sin(phi) * cos(theta); 0, sin(phi), cos(phi)] / cos(theta);
ans = J * [p; q; r];

syms a b c x y z
ans = cross([a, b, c], [x, y, z]);
