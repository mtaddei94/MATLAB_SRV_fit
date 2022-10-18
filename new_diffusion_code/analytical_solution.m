
% This is the code to solve the analytical solution of carrier density
% 1D diffusion from https://www.nature.com/articles/nenergy2016207#Sec12 --
% see sI page 6

D = 0.04; %units of cm^2/s
tau = 1e-7; %units of s
srv = 100; %units of cm/s 
x = -250e-7:1e-7:250e-7; %units of cm, the lower and upper limit are the first and last value and the range is in the middle
t = 0:1e-9:1e-6; %linspace(1*(10^-6), 5*(10^-6), 501) %units of s

sol = @(x,t)((1/2)*exp(-((x^2)/(4*D*t))))*(@(x, D, t)erf_func(-(x/(2*(sqrt(D*t))))) - @(x, D, t)erf_func((x/(2*(sqrt(D*t))))) + @(x, S, D, t)2*erf_func(sqrt(t/D) + ((2*x)/(2*sqrt(D*t)))))*@(t, tau)exp(-t/tau);


mesh(x, t, sol) % to get 3D graphic
% 
% figure
% plot(x, sol(1,:), 'LineWidth', 5); % to get a single line profile of the
% change in carrier density in space
% 
% figure
% plot(x, sol(10,:), 'LineWidth', 5);
% 
% figure
% plot(x, sol(1000,:), 'LineWidth', 5);
