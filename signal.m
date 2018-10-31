clc, clear
close all

N = 2000;
T = 1 / 80;
N_E = 4 / T;
x_r = zeros(N, 1);
y_r = zeros(N, 1);
z_r = zeros(N, 1);
for i = 1 : N
    if i < N_E
        omega = i / N_E * 0.4;
        t = i * T;
        x_r(i) = 2 * sin(omega * t);
        y_r(i) = 2 * sin(omega * t) * cos(omega * t);
        z_r(i) = -1.5;  
    else
        omega = 0.4;
        t = i * T;
        x_r(i) = 2 * sin(omega * t);
        y_r(i) = 2 * sin(omega * t) * cos(omega * t);
        z_r(i) = -1.5; 
    end
end
figurename("X_Y");
plot(x_r, y_r, '-o');
grid on
axis equal