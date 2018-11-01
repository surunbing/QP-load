function dx = bycicle_kinematics(t, x)

global g;
global mq ml;
global L;
global j11 j12 j13 j21 j22 j23 j31 j32 j33;

xl = x(1);
yl = x(2);
zl = x(3);
dxl = x(4);
dyl = x(5);
dzl = x(6);
qx = x(7);
qy = x(8);
qz = x(9);
wx = x(10);
wy = x(11);
wz = x(12);
r11 = x(13);
r12 = x(14);
r13 = x(15);
r21 = x(16);
r22 = x(17);
r23 = x(18);
r31 = x(19);
r32 = x(20);
r33 = x(21);
wwx = x(22);
wwy = x(23);
wwz = x(24);

mx = x(25);
my = x(26);
mz = x(27);
F = x(28);

dqx = qz*wy - qy*wz;
dqy = qx*wz - qz*wx;
dqz = qy*wx - qx*wy;
qfre3mqlqq = F * r13 * qx + F * r23 * qy + F * r33 * qz - L * mq * (dqx * dqx + dqy * dqy + dqz * dqz);
            
dx(1) = dxl; 
dx(2) = dyl; 
dx(3) = dzl; 
dx(4) = qfre3mqlqq * qx / (mq + ml);
dx(5) = qfre3mqlqq * qy / (mq + ml);
dx(6) = qfre3mqlqq * qz / (mq + ml) - g;
dx(7) = dqx;
dx(8) = dqy;
dx(9) = dqz;
dx(10) =  -(F * qy * r33 - F * qz * r23) / (L * mq);
dx(11) = (F * qx * r33 - F * qz * r13) / (L * mq);
dx(12) = -(F * qx * r23 - F * qy * r13)/(L * mq);
dx(13) = r12 * wwz - r13 * wwy;
dx(14) = r13 * wwx - r11 * wwz;
dx(15) = r11 * wwy - r12 * wwx;
dx(16) = r22 * wwz - r23 * wwy;
dx(17) = r23 * wwx - r21 * wwz;
dx(18) = r21 * wwy - r22 * wwx;
dx(19) = r32 * wwz - r33 * wwy;
dx(20) = r33 * wwx - r31 * wwz;
dx(21) = r31 * wwy - r32 * wwx;
beq = [mx - wy*(j31*wx + j32*wy + j33*wz) + wz*(j21*wx + j22*wy + j23*wz);
     my + wx*(j31*wx + j32*wy + j33*wz) - wz*(j11*wx + j12*wy + j13*wz);
     mz - wx*(j21*wx + j22*wy + j23*wz) + wy*(j11*wx + j12*wy + j13*wz)];
J = [j11, j12, j13; j21, j22, j23; j31, j32, j33];
res = J ^ (-1) * beq;
dx(22) = res(1);
dx(23) = res(2);
dx(24) = res(3);
dx(25) = 0;
dx(26) = 0;
dx(27) = 0;
dx(28) = 0;

dx = dx';
