function simulation_runner()
% Author: Collin Victor. Last modified on 2023-09-12
% Solves the 2D Navier-Stokes equations in stream-function form
% Using spectral methods, IMEX Euler with integrating factor, 2/3
% dealiasing

%Remove plots
close all;
clear all; %#ok<CLALL>

%Add helpers directory to path, this contains useful utility functions
addpath(genpath('helpers'));

%Add routines directory to path, this contains functions used for various
% routines, such as saving and restarting.



% Restart calls go here if needed
% [p, settings, psi_hat, vars] = restartRoutine();

tol = 0.1;


% If not restarting a run, start as normal.
if ~exist('p','var') %If p exists, you are running from a restart.
    %Initialize settings.
    % If empty use default values, otherwise specify field and value in pairs e.g.
    % settings = Settings();
    settings = Settings('save_data', true,'save_interval', 5, 'save_tag', "adaptive_alt",'plot_spectrum_ensembles', false, 'noise_enable', true, 'noise_level', 0.1, 'noise_type','white', 'plot_skip', 10, 'plot_ref', false, 'plot_error', true, 'plot_energy', false, 'final_time', 1);
    % settings = Settings('plot_skip',10, 'noise_enable',true, 'save_data', false, 'save_interval', 0.15);
    % settings = Settings('zero_initial_data', true,'plot_enable', true, 'plot_ref', true, 'plot_spectrum', true, 'plot_energy', true, 'plot_error', false, 'waitbar', false, 'final_time', 1, 'plot_skip', 1);


    %Initialize parameters.
    % If empty use default values, otherwise specify field and value in pairs e.g.
    %p = Parameters(settings, 'G', 1e5,'resolution', 2^9)
    p = Parameters(settings, 'resolution', 2^10, 'trunc', 100);



    %Initialize data assimilation variables
    vars = [];
    vars = repelem(DA_obs(p),1,1);

    % Data assimilation trials should be initialized similar to Settings
    % and Parameters, i.e. with optional inputs inserted in pairs. Only p is
    % a required input.
    vars(1) = DA_obs(p,"IC_type", "Projection", "assimilate_type", "AOT", "mu", 10000); %Default DA_obs object is made
    % vars(1) = DA_obs(p,"IC_type", "Zero", "assimilate_type", "AOT", "mu", 1 ); %Default DA_obs object is made

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


psi_hat_alt = zeros(size(psi_hat));
% psi_hat_alt(vars(1).trunc_index) = psi_hat(vars(1).trunc_index);

psi = ifftn(psi_hat); %#ok<NASGU>

size_vars = length(vars);


ev = evolution_parameters(p);

noisy_obs = psi_hat;
if(settings.noise_enable)
    noisy_obs = generate_noisy_obs(psi_hat, p, settings);
end

for i = 1: size_vars
    vars(i) = vars(i).initialize(noisy_obs);

    if vars(i).assimilate_type == "EnKF"
        M = vars(i).ensemble_count;
           ens_mean = zeros(p.Nx, p.Ny);

        for j = 1:M 
            noisier_obs = generate_noisy_obs(noisy_obs, p, settings);
            noisier_obs(vars(i).trunc_index_comp) = 0;
            vars(i).ensemble(j).forecast = noisier_obs;
            vars(i).ensemble(j).analysis = noisier_obs;
            ens_mean = ens_mean + noisier_obs;

        end
        ens_mean = ens_mean./M;
        vars(i).v_hat = ens_mean;
    end
end



if(~isempty(vars))
    end_flag = zeros(1,size_vars);
end



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


    psi_hat_old = psi_hat;
    [en, psi_hat] = evolve(p, ev, psi_hat);


     % psi_hat_ = psi_hat;
    [~, psi_hat_alt] = evolve(p, ev, psi_hat_alt);


    if(settings.noise_enable)
        noisy_obs_old = noisy_obs;
        noisy_obs = generate_noisy_obs(psi_hat, p, settings);
    end

    % energy(ti) = en;
    p.energy(ti) = en;

    % tic;
    for i = 1:size_vars
            if(settings.noise_enable)
                vars(i) = vars(i).observe(noisy_obs_old, noisy_obs);
            else
                vars(i).observe(psi_hat_old, psi_hat);
            end
    end

    for i =1:size_vars

            %Forward Euler
            [en,k1] = compute_rhs_hat(p.time_current, vars(i).v_hat, p);
            v_pred = ev.E2.*vars(i).v_hat + ev.dt*ev.E2.*k1;

            vars(i).energy(ti) = en;

            vars(i) = vars(i).assimilate(ev, v_pred, p.time_current);
    end

    % toc;

    if settings.save_data && (mod(p.time_current, settings.save_interval) < p.dt && mod(p.time_current, settings.save_interval) >= 0)
        saveRoutine(settings, p, plots, psi_hat, vars);
    

    end

    for i = 1:size_vars
        if(end_flag(i) == 0)
            psi_hat_noise = psi_hat;
            psi_hat_noise(vars(1).trunc_index) = noisy_obs(vars(1).trunc_index);
            % if()
            total_diff = abs(vars(i).v_hat-psi_hat);

            obs_diff = abs(psi_hat_noise - vars(i).v_hat);
            low_obs_diff = norm(obs_diff(vars(i).trunc_index))/sqrt(p.Nx*p.Ny);

            low_diff = norm(total_diff(vars(i).trunc_index))/sqrt(p.Nx*p.Ny);
            high_diff = norm(total_diff(vars(i).trunc_index_comp))/sqrt(p.Nx*p.Ny);
            low_diff = max(low_diff, p.eps);

            vars(i).error_low(ti) = low_diff;
            vars(i).error_high(ti) = high_diff;
            vars(i).error(ti) = low_diff + high_diff;
            vars(i).error_observed(ti) = low_obs_diff;
        end
        if vars(i).assimilate_subtype == "Adaptive"
            if ti > 7 && vars(i).switch_times_mu(end) < p.t(ti-5)
                prev_error = vars(i).error_observed(ti-5);
                current_error = vars(i).error_observed(ti);
                diff = log(current_error) - log(prev_error);
                if diff > tol
                    vars(i).mu = vars(i).mu/10;
%                     last_ti = ti;
                    vars(i).mu_values = [vars(i).mu_values; vars(i).mu];
                    vars(i).switch_times_mu = [vars(i).switch_times_mu; p.t(ti)];



                end

            end


        end

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
function noisy_obs = generate_noisy_obs(psi_hat,p, settings)
noise_v1 = generate_non_div_free_noise(p, settings);
noise_v2 = generate_non_div_free_noise(p, settings);
[v1,v2] = psi_converter(psi_hat,p);
noisy_v1 = v1 + noise_v1;
noisy_v2 = v2 + noise_v2;

noisy_obs = velocity_converter(noisy_v1,noisy_v2,p);

end

function noise = generate_non_div_free_noise(p, settings)
% Define the size of your 2D domain
N = p.Nx;

% Generate the noise in Fourier space
noise_fourier_space = settings.noise_level*(randn(N, N/2+1) + 1i * randn(N, N/2+1));

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
