clc, clear
close all

syms wx wy wz;
syms xl yl zl;
syms dxl dyldzl;
syms qx qy qz;
syms wx wy wz;
syms r11 r12 r13 r21 r22 r23 r31 r32 r33;
syms wwx wwy wwz;


syms mx my mz;
syms f;

syms xr yr zr;
syms g;
syms mq ml;
syms L;
syms j11 j12 j13 j21 j22 j23 j31 j32 j33;

omega = [wx; wy; wz];
q = [qx; qy; qz];
ans = cross(omega, q);

syms dqx dqy dqz;
dq = [dqx; dqy; dqz];

R = [r11 r12 r13; r21 r22 r23; r31 r32 r33];
J = [j11, j12, j13; j21, j22, j23; j31, j32, j33];
e3 = [0; 0; 1];

OMEGA = [0 -wwz wwy; wwz 0 -wwx; -wwy wwx 0];

ans = (dot(q, f * R * e3) - mq * L * dot(dq, dq)) * q / (mq + ml) - g * e3;

ans = -cross(q, f * R * e3) / (mq * L);

ans = R * OMEGA;

ans = cross(omega, J * omega);

syms dwwx dwwy dwwz;
DOMEGA = [dwwx; dwwy; dwwz];
M = [mx; my; mz];
ans = J * DOMEGA + cross(omega, J * omega) - M;

