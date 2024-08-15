function self_synchronous(theta_input)
% Author: Collin Victor. Last modified on 2023-09-12
% Solves the 2D Navier-Stokes equations in stream-function form
% Using spectral methods, IMEX Euler with integrating factor, 2/3
% dealiasing

%Remove plots
% close all;
% clear all; %#ok<CLALL>

%Add helpers directory to path, this contains useful utility functions
addpath(genpath('helpers'));

%Add routines directory to path, this contains functions used for various
% routines, such as saving and restarting.



% Restart calls go here if needed
% [p, settings, psi_hat, vars] = restartRoutine();

% If not restarting a run, start as normal.
if ~exist('p','var') %If p exists, you are running from a restart.
    %Initialize settings.
    % If empty use default values, otherwise specify field and value in pairs e.g.
    % settings = Settings();
    settings = Settings('save_data', true,'save_interval', 10, ...
        'save_tag', "SS_NSE_noise=5_theta_tests",'plot_spectrum_ensembles', false,...
        'noise_enable', true, 'noise_level', 1e-5, 'noise_type','white',...
        'plot_skip', 100, 'plot_ref', false, 'plot_error', true, 'plot_energy', true,...
        'final_time', 25);
    % settings = Settings('plot_skip',10, 'noise_enable',true, 'save_data', false, 'save_interval', 0.15);
    % settings = Settings('zero_initial_data', true,'plot_enable', true, 'plot_ref', true, 'plot_spectrum', true, 'plot_energy', true, 'plot_error', false, 'waitbar', false, 'final_time', 1, 'plot_skip', 1);


    %Initialize parameters.
    % If empty use default values, otherwise specify field and value in pairs e.g.
    %p = Parameters(settings, 'G', 1e5,'resolution', 2^9)
    p = Parameters(settings, 'resolution', 2^9, 'trunc', 50);


    %Initialize data assimilation variables
    vars = [];
    vars = repelem(DA_obs(p),1,1);


    psi_hat = p.psi_hat_initial;

    en = generate_energy(p, psi_hat);
    p.energy(1) = en;

    for i = 1: length(vars)
        en = generate_energy(p, vars(i).v_hat);
        vars(i).energy(1) = en;
    end


end



if(settings.waitbar)
    wb = waitbar(0,'Initializing.');
end


psi = ifftn(psi_hat); %#ok<NASGU>

size_vars = length(vars);


ev = evolution_parameters(p);



if(~isempty(vars))
    end_flag = zeros(1,size_vars);
end






theta_1 = theta_input; % 0 1/4 1/2 3/4 1
theta_2 = 1 - theta_1;


vars(1).mu = theta_2;

psi_hat_2 = psi_hat;
variance = settings.noise_level;
noise = generate_noisy_obs(p, variance);

psi_hat_2 = psi_hat_2 + noise;
% psi_hat_2 = 0.*psi_hat_2;
% psi_hat_2(vars(1).trunc_index) = psi_hat(vars(1).trunc_index);


vars(1).v_hat = psi_hat_2;
psi_hat_1 = psi_hat;




%% Plots
if(settings.plot_enable)
    plots = MasterPlotter(settings,p, psi_hat, vars);
end




if(settings.save_data)
saveInitialize(settings, p, plots, psi_hat, vars);
end

%% Main Program Loop
for ti = p.ti+1:p.num_timesteps
    p.ti = ti;
    p.time_current = p.t(ti);


    if(settings.waitbar)
        waitbar(ti/p.num_timesteps,wb,'Running the main loop.');
    end

    if (settings.plot_enable && mod(ti-1, settings.plot_skip)==0)
       plots = plots.update(p, psi_hat, vars);
    end

    [energy, k1, nonlin_1] = compute_rhs_hat_ss(p.time_current, psi_hat_1, p);
    [energy_2, k2, nonlin_2] = compute_rhs_hat_ss(p.time_current, psi_hat_2, p);

    psi_hat_1 = ev.E2.*psi_hat_1 + ev.dt*ev.E2.*k1 + ev.dt*ev.E2.*(theta_1*p.trunc_array.*( nonlin_1 - nonlin_2));
    psi_hat_2 = ev.E2.*psi_hat_2 + ev.dt*ev.E2.*k2 + ev.dt*ev.E2.*(theta_2*p.trunc_array.*( nonlin_2 - nonlin_1));





    vars(1).v_hat = psi_hat_2;
    vars(1).v_obs = psi_hat_2;
    vars(1).u_obs = psi_hat_1;





  
    if settings.save_data && (mod(p.time_current, settings.save_interval) < p.dt && mod(p.time_current, settings.save_interval) >= 0)
        % [p.time_current, settings.save_interval, p.dt]
        saveRoutine(settings, p, plots, psi_hat, vars);
    

    end
    

    for i = 1:size_vars
        errorcheck = 1; %#ok<NASGU>
        if(end_flag(i) == 0)
            % if()
            total_diff = abs(vars(i).v_hat-psi_hat_1);
            low_diff = norm(total_diff(vars(i).trunc_index))/sqrt(p.Nx*p.Ny);
            high_diff = norm(total_diff(vars(i).trunc_index_comp))/sqrt(p.Nx*p.Ny);
            low_diff = max(low_diff, p.eps);

            vars(i).error_low(ti) = low_diff;
            vars(i).error_high(ti) = high_diff;
            vars(i).error(ti) = low_diff + high_diff;
        end
        slope = (log(vars(i).error(ti)) - log(vars(i).error(ti-1)))/p.dt;
        vars(1).energy(ti) = slope;



    end

    
    if(~size_vars == 0)
        if(sum(end_flag) == size_vars)
            break;
        end
    end

end

if settings.waitbar
        waitbar(1,wb,'Finishing details.');
end

if settings.save_data
        saveFinalize(settings, p, plots, psi_hat, vars);
end


if(settings.waitbar)
    close(wb)
end


close all;

end

function [en, k1, nonlin] = compute_rhs_hat_ss(time, psi_hat, p)
% Dealias by setting the "middle" Fourier coefficients to zero.
% This is where Matlab stores the highest frequencies.
psi_hat(p.dealias_modes) = 0;
% psi_2(p.dealias_modes) = 0;

% enstrophy_sq = u_x_hat*u_x_hat'/p.Nsq;

% psi_x_hat = bsxfun(@times,p.ikx,psi_hat);
% psi_y_hat = bsxfun(@times,p.iky,psi_hat);
psi_x_hat = p.IKX.*psi_hat;
psi_y_hat = p.IKY.*psi_hat;


en = (norm(psi_x_hat)^2 + norm(psi_y_hat)^2)/p.Nx/p.Ny;
%energy = 0;

% % Compute the nonlinear terms in physical space using dealiased versions.
% psi_x      = ifftn(              psi_x_hat,'symmetric'); % psi_x = -v
% psi_y      = ifftn(              psi_y_hat,'symmetric'); % psi_y =  u
% psi_x_lap  = ifftn(p.k_lap.*psi_x_hat,'symmetric');
% psi_y_lap  = ifftn(p.k_lap.*psi_y_hat,'symmetric');
% % rhs_hat = p.helm_inv.*p.k_lap_inv.*fft2(-psi_y.*psi_x_lap + psi_x.*psi_y_lap + p.f_hat);
% rhs_hat = p.k_lap_inv.*(fft2(-psi_y.*psi_x_lap + psi_x.*psi_y_lap) + p.f_hat);

% Basdevant formula (4 fft/ifft operations, 6 products). See:
%   P. Emami and J. C. Bowman,"On the Global Attractor of 2D Incompressible Turbulence with Random Forcing",
%   Journal of Differential Equations. 264, 4036-4066 (2018).
u1 =  ifftn(psi_y_hat,'symmetric');
u2 = -ifftn(psi_x_hat,'symmetric');
% u1       =  ifftn(bsxfun(@times,p.iky,psi_hat),'symmetric'); % psi_y =  u
% u2       = -ifftn(bsxfun(@times,p.ikx,psi_hat),'symmetric'); % psi_x = -v
nonlin = -p.k_lap_inv.*(p.k1sq_m_k2sq.*fftn(u1.*u2) + p.k1k2.*fftn(u2.^2 - u1.^2));
k1 = -nonlin + p.k_lap_inv.*p.f_hat;
% k1 = p.k_lap_inv.*(p.k1sq_m_k2sq.*fftn(u1_1.*u2_1) + p.k1k2.*fftn(u2_1.^2 - u1_1.^2) + p.f_hat);

% nonlin = p.k1sq_m_k2sq.*fftn(u1_1.*u2_1) + p.k1k2.*fftn(u2_1.^2 - u1_1.^2);
% nonlin_term_2 = p.k1sq_m_k2sq.*fftn(u1_2.*u2_2) + p.k1k2.*fftn(u2_2.^2 - u1_2.^2);
% k1 = p.k_lap_inv.*(nonlin_term_1 - nonlin_term_2);
% k2 = p.k_lap_inv.*(nonlin_term_2 - nonlin_term_1);



% rhs_hat = p.helm_inv.*p.k_lap_inv.*(p.k1sq_m_k2sq.*fftn(u1.*u2) + p.k1k2.*fftn(u2.^2 - u1.^2) + p.f_hat);

% if p.compute_norms == 1
%     % Compute expensive norms while we have the data.
%     p.compute_norms = 0;
%     u_Linf = sqrt(max(u1(:).^2 + (u2(:)).^2));
% else
%     u_Linf = 0;
% end
end

function ev = evolution_parameters(p)

%Generate all of the parameters necessary for evolving solution via
%numerical scheme.
E  = exp(p.nu*p.k_lap.*p.dt/2);
E2 = exp(p.nu*p.k_lap.*p.dt);
E3 = exp(p.nu*p.k_lap.*2*p.dt);
E4 = exp(p.nu*p.k_lap.*3*p.dt);
% E_test = E2;
% E_test(p.trunc_array == 1) = 1;
E_test = exp(1.5409e-4*p.k_lap.*p.dt);

% E2_alpha = exp(p.nu_alpha*p.k_lap.^alpha.*p.dt);

%TODO: Implement a variety of these based on the numerical scheme chosen
%using a switch statement a p.evolution_scheme
ev = struct;
ev.type = p.evolution_scheme;
ev.E = E;
% ev.E2_alpha = E2_alpha;
ev.E2 = E2;
ev.E3 = E3;
ev.E4 = E4;
ev.dt = p.dt;
ev.E_test = E_test;
end

function noisy_obs = generate_noisy_obs(p, variance)
noise_v1 = generate_non_div_free_noise(p, variance);
noise_v2 = generate_non_div_free_noise(p, variance);
% [v1,v2] = psi_converter(psi_hat,p);
noisy_v1 = noise_v1;
noisy_v2 = noise_v2;

noisy_obs = velocity_converter(noisy_v1,noisy_v2,p);

end


function noise = generate_non_div_free_noise(p, variance)
% Define the size of your 2D domain
N = p.Nx;

% Generate the noise in Fourier space
noise_fourier_space = variance*(randn(N, N/2+1) + 1i * randn(N, N/2+1));

% Set the 0 wave mode to 0
noise_fourier_space(1, 1) = 0 + 0i;

% Create the complex conjugate mirrored version of your noise
noise_fourier_space_full = zeros(N, N);
noise_fourier_space_full(:, 1:N/2+1) = noise_fourier_space;
% Flipud function flips the array up down, and fliplr function flips it left to right

% conj function takes the complex conjugate
noise_fourier_space_full(:, N/2+2:N) = rot90(conj(noise_fourier_space(:, 2:N/2)),2);

noise = noise_fourier_space_full;

end
