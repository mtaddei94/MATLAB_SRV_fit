

diff_coeff = 0.004; %units of cm^2/s
mono_recomb_coeff = 1e6; %units of s^-1
bi_recomb_coeff = 4e-11; %units of cm/s
init_carrier_density = 1e13; %units of cm^-1
srv_1 = 1000; %units of cm/s 
srv_2 = 1000; %units of cm/s 
film_thickness = 1000; %in nm




% distance array: x values are a linear sequence (= linspace in python)
xmesh_step = 1; % distance step in nm. It is going to be in the middle of the min and max of the array (min:step:max) shpwn below
x = 0:(xmesh_step*1e-7):(film_thickness*1e-7); %units of cm, the lower and upper limit are the first and last value and the range is in the middle

%time array
t_step = 1e-10; % time step in s
t_simulate = 0:t_step:1e-7; %linspace(1*(10^-6), 5*(10^-6), 501) %units of s

sol= pdepe(0,@(x,t,u,DuDx)mono_bi_recomb_pde(x,t,u,DuDx, mono_recomb_coeff,bi_recomb_coeff,diff_coeff),@(x)ic_exponential(x, init_carrier_density, absorbance, film_thickness),@(xl,ul,xr,ur,t)SRV_fixed_bc(xl,ul,xr,ur,t, srv_1, srv_2),x,t_simulate);

t_size = size(t_simulate); % find size of time array
n_integrate = zeros(1,t_size(2)); %create empty array to store x integrated solutions for each time point
    
for j = 1:t_size(2);
    n_int = trapz(x, sol(j,:)); %integrate over x for time point j
    n_integrate(1,j) = n_int; %store in array        
end

figure
plot(t_simulate,n_integrate/max(n_integrate))

%Calculate spatially integrated PL decay from solution to time dependent carrier density

pl_x_integrated = pl_function_mono_bi(t_simulate,n_integrate, mono_recomb_coeff, bi_recomb_coeff);

pl_x_integrated_norm = pl_x_integrated / max(pl_x_integrated);


simulated_noise = sqrt(poissrnd(pl_x_integrated * t_step));


noisy_pl = simulated_noise;


x0 = [0,1];
lb = [10,10];
ub = [1000,1000];
x_fit = lsqcurvefit(pl_function_mono_bi, x0, t_simulate, pl_x_integrated, lb, ub);

plot(t_simulate,x,'ko',xdata,pl_function_mono_bi(x_fit,pl_x_integrated),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')


r= optimvar('r', 3, "LowerBound",10, "UpperBound", 1000);
