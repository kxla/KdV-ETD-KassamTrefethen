% KT_EDTRK4_conv.m
%
% Convergence study of exponential time differencing method for solving
% the KdV equation given by Kassam & Trefethen (2005) using ETDRK4 scheme
%
% u_t = (1/6) * epsilon * u_xxx - (F-1) * u_x + (3/2) * alpha * u * u_x

clc
clear
close all

h_values = [1/4 10e-2 10e-3 10e-4 10e-5];
N_values = [256 512 1024 2048 4096 8192];
%u_exact = zeros(size(h_values,2),size(N_values,2),N);
%u_numerical = zeros(size(h_values,2),size(N_values,2),N);
infnorm = zeros(size(h_values,2),size(N_values,2));
twonorm = zeros(size(h_values,2),size(N_values,2));
dx = zeros(size(h_values,2),size(N_values,2));
dt = zeros(size(h_values,2),size(N_values,2));


for index1 = 1:size(h_values,2)
    for index2 = 1:size(N_values,2)
        h = h_values(index1);
        N = N_values(index2);
        %u_exact = KT_ETDRK4(h,N);           % Exact solution
        [u_numerical, u_exact] = KT_ETDRK4(h,N);   % Numerical solution
        infnorm(index1,index2) = norm(u_exact - u_numerical,Inf);
        twonorm(index1,index2) = norm(u_exact - u_numerical)/sqrt(N);
        dx(index1,index2) = N;
        dt(index1, index2) = h;
        clc
        % close all
        clear u_exact u_numerical
    end
end


figure
surf(dx,dt,infnorm)
title('Infinity norm')
xlabel('dx')
ylabel('dt')
zlabel('Error')
figure
surf(dx,dt,twonorm)
title('Two norm')
xlabel('dx')
ylabel('dt')
zlabel('Error')