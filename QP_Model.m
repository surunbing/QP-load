function dx = bycicle_kinematics(t, x)

global mF ml l AL cDL g rho p q r fT fHx fHy fHz
xx = x(1);
yy = x(2);
zz = x(3);
phi = x(7);
theta = x(8);
psi = x(9);
vx = x(4);
vy = x(5);
vz = x(6);
vxl = x(12);
vyl = x(13);
xl = x(10);
yl = x(11);

dot_x = cos(psi) * cos(theta) * vx + (cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)) * vy + (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) * vz;
dot_y = cos(theta) * sin(psi) * vx + (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) * vy + (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)) * vz;
dot_z = -sin(theta) * vx + cos(theta) * sin(phi) * vy + cos(phi) * cos(theta) * vz;
ddot_vx = r * vy - q * vz - g * sin(theta) + (fHx * cos(psi) * cos(theta) - fHz * sin(theta) + fHy * cos(theta) * sin(psi)) / mF;
ddot_vy = p * vz - r * vx + g * cos(theta) * sin(phi) +  (fHy * (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) - fHx * (cos(phi) * sin(psi) - cos(psi) * sin(phi) * sin(theta)) + fHz * cos(theta) * sin(phi)) / mF;
ddot_vz = q * vx - p * vy + g * cos(phi) * cos(theta) -(fT - fHx * (sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)) + fHy * (cos(psi) * sin(phi) - cos(phi) * sin(psi) * sin(theta)) - fHz * cos(phi) * cos(theta))/mF;
ddot_x = cos(psi) * cos(theta) * ddot_vx + (cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)) * ddot_vy + (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) * ddot_vz;
ddot_y = cos(theta) * sin(psi) * ddot_vx + (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) * ddot_vy + (cos(phi) * sin(psi) * sin(theta) - cos(psi) * sin(phi)) * ddot_vz;
ddot_z = -sin(theta) * ddot_vx + cos(theta) * sin(phi) * ddot_vy + cos(phi) * cos(theta) * ddot_vz;


delta_x = xx - xl;
delta_y = yy - yl;
zl = zz + sqrt(l * l - delta_x * delta_x - delta_y * delta_y);
delta_z = zz - zl;
dot_delta_x = dot_x - vxl;
dot_delta_y = dot_y - vyl;
H = dot_delta_x * delta_x + dot_delta_y * delta_y;
dot_delta_z = H / delta_z;
dot_zl = dot_delta_z + dot_z;
vl_square = vxl * vxl + vyl * vyl + dot_zl * dot_zl;
fDL_m = -0.5 * cDL * rho * AL * vl_square;
fDLx = fDL_m * vxl;
fDLy = fDL_m * vyl;
fDLz = fDL_m * dot_zl;
rcx = xl - xx;
rcy = yl - yy;
rcz = zl - zz;
% 方程的另一边
diff_r = dot_delta_x * dot_delta_x + dot_delta_y * dot_delta_y + dot_delta_z * dot_delta_z + ddot_x * delta_x + ddot_y * delta_y + ddot_z * ddot_z;
fp = fDLx;
fq = fDLy;
fr = fDLz; + ml * g;
b_eq = [rcz * fq - rcy * fr;
        rcz * fp + rcx * fr;
        rcy * fp - rcx * fq;
        diff_r];
A_eq = [0, rcz * ml, -rcy * ml;
        -rcz * ml, 0, rcx * ml;
        rcy * ml, -rcx * ml, 0;
        delta_x, delta_y, delta_z];
res = linsolve(A_eq, b_eq);
alx = res(1);
aly = res(2);



dx(1,1) = dot_x;
dx(2,1) = dot_y;
dx(3,1) = dot_z;
dx(4,1) = ddot_vx; %r * vy - q * vz - g * sin(theta) + (fHz * cos(psi) * cos(theta) - fHz * cos(theta) * sin(psi)) / mF;
dx(5,1) = ddot_vy; %p * vz - r * vx + g * cos(theta) * sin(phi) +  (fHy * (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) - fHz * (cos(phi) * sin(psi) - cos(psi) * sin(phi) * sin(theta)) + fHz * cos(theta) * sin(phi)) / mF;
dx(6,1) = ddot_vz; %q * vx - p * vy + g * cos(phi) * cos(theta) -(fT + fHy * (cos(psi) * sin(phi) - cos(phi) * sin(psi) * sin(theta)) - fHz * (sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)) - fHz * cos(phi) * cos(theta)) / mF;
dx(7,1) = p + (r * cos(phi) * sin(theta)) / cos(theta) + (q * sin(phi) * sin(theta)) / cos(theta);
dx(8,1) = q * cos(phi) - r * sin(phi);
dx(9,1) = (r * cos(phi)) / cos(theta) + (q * sin(phi)) / cos(theta);
dx(10,1) = vxl;
dx(11,1) = vyl;
dx(12,1) = alx;
dx(13,1) = aly;