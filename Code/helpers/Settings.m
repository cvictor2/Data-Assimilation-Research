classdef Settings
    properties
        noise_enable = false; % Add noise to simulations?
        noise_level = 0.1; % Variance of noise
        noise_type = 'white'; % White or pink

        plot_enable = true; % Plot something
        plot_vars = true; % Plot DA 2x2 plots for each trial
        plot_ref = true; % Plot reference solution
        plot_energy = true; % Plot energy
        plot_skip = 10; % Timestep interval to skip
        plot_spectrum = true; % Plot spectrum in separate window
        plot_error = true; % Plot errors in separate window
        plot_pause = false;
        plot_error_debug = false;
        plot_spectrum_ensembles = false;


        movie = false; % Save plots into gif file

        seed

        restart_extend = 10;
        
        zero_forcing = false;

        initial_random = false; % Should initial state be randomized
        initial_time = 0; % Time to ramp up to initially

        zero_initial_data = false; %Start with 0 initial data.

        trunc = 50; % Number of Fourier modes to truncate to
        final_time = 100; % Final time to run out to

        save_data = true; % Save anything?
        save_interval = 5; % time between saving files
        save_tag = ''; % Naming convention for filenames
        % save_time = datetime("now",'Format','yyyy_mm_dd_HH_MM_SS');
        currentDateTime = datetime('now');
        save_time = "";

        save_var = true; % Save variables
        save_ref = true; % Save reference solution?
        save_plots = true; % Save plotted figures?


        waitbar = false;

        save_name %name of folder inside runs that data will be saved in.

        note = ''; % Additional string used for save_tag and other misc. uses
    end

    % properties (Dependent)
    %     seed = randi([1,1e6]);
    % end

    methods
        % Constructor with input parser
        function obj = Settings(varargin)
            %Note input parser lets user input whatever

            % Create a 'random' initial seed (this will be overwritten
            % based on the input settings).
            rng('shuffle');
            obj.seed = randi([1,1e6]);

            if isempty(obj.save_tag)
                default_name = obj.save_time;
            else
                default_name = strcat(obj.save_time, "_", obj.save_tag);
            end
            % Create an input parser object
            p = inputParser;

            % Define default values for the properties
            addParameter(p, 'noise_enable', obj.noise_enable);
            addParameter(p, 'noise_level', obj.noise_level);
            addParameter(p, 'noise_type', obj.noise_type);
            addParameter(p, 'plot_enable', obj.plot_enable);
            addParameter(p, 'plot_ref', obj.plot_ref);
            addParameter(p, 'plot_vars', obj.plot_vars);
            addParameter(p, 'plot_energy', obj.plot_energy);            
            addParameter(p, 'plot_error', obj.plot_error);
            addParameter(p, 'plot_spectrum', obj.plot_spectrum);
            addParameter(p, 'plot_spectrum_ensembles', obj.plot_spectrum_ensembles);
            addParameter(p, 'plot_skip', obj.plot_skip);
            addParameter(p, 'plot_error_debug', obj.plot_error_debug);
            addParameter(p, 'movie', obj.movie);
            addParameter(p, 'initial_random', obj.initial_random);
            addParameter(p, 'initial_time', obj.initial_time);
            addParameter(p, 'trunc', obj.trunc);
            addParameter(p, 'final_time', obj.final_time);
            addParameter(p, 'save_data', obj.save_data);
            addParameter(p, 'save_interval', obj.save_interval);
            addParameter(p, 'save_tag', obj.save_tag);
            addParameter(p, 'save_var', obj.save_var);
            addParameter(p, 'save_ref', obj.save_ref);
            addParameter(p, 'save_plots', obj.save_plots);
            addParameter(p, 'note', obj.note);
            addParameter(p, 'waitbar', obj.waitbar);
            addParameter(p, 'seed', obj.seed);
            addParameter(p, 'save_name', string(obj.save_name));
            addParameter(p, 'restart_extend',obj.restart_extend);
            addParameter(p, 'zero_initial_data',obj.zero_initial_data);
            addParameter(p, 'plot_pause',obj.plot_pause);
            addParameter(p, 'zero_forcing', obj.zero_forcing);



            % Parse the inputs
            parse(p, varargin{:});

            % Assign the parsed values to the properties
            fields = fieldnames(p.Results);
            for i = 1:length(fields)
                obj.(fields{i}) = p.Results.(fields{i});
            end
            %Seed is definitely set, so now set the rng.
            % rng(obj.seed);
            
            obj.save_time = datestr(obj.currentDateTime, 'yyyy_mm_dd_HH_MM_SS');


            %Initialize file name for saving, by defualt its just the
            %timestamp.
            if isempty(obj.save_name)
                if isempty(obj.save_tag)
                    obj.save_name = string(obj.save_time);
                else
                    obj.save_name = strcat(obj.save_time, "_", obj.save_tag);
                end
            end

        end

        function obj = set.seed(obj, value)
            %Setter for seed. This is called anytime seed is changed.
            %This is necessary as anytime the seed is changed rng must be
            %called.
            obj.seed = value;
            rng(value); % Set the rng every time seed is changed
        end


    end
end