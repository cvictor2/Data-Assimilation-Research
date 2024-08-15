close all;

% Specify the folder where the files live.
myFolder = './runs'; % Replace with your folder path

plotTitle = 'Error Decay Rates for Mutual Nudging';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end
load('/Users/victor/Documents/MATLAB/Synchronization/Code/runs/2024_08_03_14_46_29_different_IC/parameters.mat')
p = parameters;
filePattern = fullfile(myFolder, '*zero_convergence/DA/*.mat'); % Change to whatever pattern you need.

% Get a list of all files in the folder with the desired file name pattern.
% filePattern = fullfile(myFolder, '*_Mutual_Nudging_NSE_*percent_final'); % Change to whatever pattern you need.
% filePattern = fullfile(myFolder, '*AOT_aot.mat'); % Change to whatever pattern you need.
theFiles_AOT = dir(filePattern);

% % Get a list of all files in the folder with the desired file name pattern.
% filePattern = fullfile(myFolder, '*SS_NSE_noise=15_*_tests'); % Change to whatever pattern you need.
% % filePattern = fullfile(myFolder, '*Synchronization_aot.mat'); % Change to whatever pattern you need.
% theFiles_Synch = dir(filePattern);

% filePattern = fullfile(myFolder, '*adaptive*_aot.mat'); % Change to whatever pattern you need.
% theFiles_Hybrid = dir(filePattern);

theFiles = [theFiles_AOT];
% theFiles = [theFiles_AOT; theFiles_Synch; theFiles_Hybrid];
% theFiles = [theFiles_AOT; theFiles_Synch];

plot_fig1 = figure;
hold on;
plot_fig2 = figure;
hold on;
plot_fig3 = figure;
hold on; % Holds the plot so all data can be plotted on same graph


legendLabels = {};
% data = [];

% Define lower threshold for plotting
tol = 1e-16;

% Define marker spacing and base offset
markerSpacing = 500;
baseOffset = 400;

mu_files = containers.Map('KeyType','double','ValueType','any');
params_files = containers.Map('KeyType','double','ValueType','any');
settings_files = containers.Map('KeyType','double','ValueType','any');
ens_files = containers.Map('KeyType', 'double', 'ValueType','any');

for k = 1 : length(theFiles)

    baseFileName = theFiles(k).name;
    % baseFileName = "aot_data.mat";
    % baseFileName = "DA/DA1.mat";
    % baseFileName = "enkf_data.mat";
    fullFileName = fullfile(theFiles(k).folder + "/" +theFiles(k).name)%, baseFileName);
    % fullFileName_params = fullfile(theFiles(k).folder + "/" +theFiles(k).name, "parameters.mat");
    % fullFileName_settings = fullfile(theFiles(k).folder + "/" +theFiles(k).name, "settings.mat");
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as calling load() to load the data.
    data = load(fullFileName);
    % data = load(fullFileName_settings);



    % elseif data.aot.assimilate_subtype ~= []
 
    mu = data.var.mu;
    if data.var.assimilate_type == "Synchronization"
        mu = inf;
    end
    % mu = data.settings.noise_level;
    % mu = data.aot.mu;
    % ens_num = data.enkf.num_ensembles;
    % ens_num = data..inf_variance;
    % end

    if isKey(mu_files,mu)
        mu_files(mu) = [mu_files(mu), {fullFileName}];
        % params_files(mu) = [params_files(mu), {fullFileName_params}];
        % settings_files(mu) = [settings_files(mu), {fullFileName_settings}];
    else
        mu_files(mu) = {fullFileName};
        % params_files(mu) = {fullFileName_params};
        % settings_files(mu) = {fullFileName_settings};
    end


    % if isKey(ens_files,ens_num)
    %     ens_files(ens_num) = [ens_files(ens_num), {fullFileName}];
    %     params_files(ens_num) = [params_files(ens_num), {fullFileName_params}];
    % else
    %     ens_files(ens_num) = {fullFileName};
    %     params_files(ens_num) = {fullFileName_params};
    % end




    %    missing_mu = ~ismember(mu, mu_values);
    %     if missing_mu
    %         mu_values = [mu_values, mu];
    %     end
end


ens_values = cell2mat(keys(ens_files));

ens_values = sort(ens_values,'descend');

% 
mu_values = cell2mat(keys(mu_files));

mu_values = sort(mu_values,'ascend');


target_ens = ens_values;
% % Specify specific mu values for plotting
% % target_mus = [1,1000,1000000,1000000000,3];
target_mus = mu_values;
% target_mus = [0,-5,-15,-25,-30,-35,-40,-45,-50];
% target_mus = [0,5,25,50];
% target_mus = [0 10 20 30 40 50];  %45 49 49.9
% % target_mus = [1, 10 25, 50, 75, 100, 125, 150, 175, 199];
% target_mus(target_mus < 1) = [];
% % target_mus = [-inf, target_mus];
% % target_mus = [-inf];
target_mus =  sort(target_mus,'descend');
colors = jet(length(target_mus)); % Generates unique colors

% Define the names for the 4 archetypes
% names = {'\mu = 1', '\mu = 1000', '\mu = 1000000', '\mu = 1000000000','Synchornization'};



% Set up coloring for all figures for each mu value

% Initialize the handles array for the 4 archetypes
h1 = gobjects(length(target_ens),1);
h2 = gobjects(length(target_ens),1);
h3 = gobjects(length(target_ens),1);

mu_values = target_mus;
% ens_values = target_ens;

% figure(plot_fig1);
% hold on;
for i = 1: length(mu_values)
    if  abs(mu_values(i)) ~=inf
        if mod(mu_values(i),1) ~= 0

            % name = sprintf('$\\theta_2 = %.2f$', mu_values(i)); % Replace 'name' with your property name

            name = sprintf('$\\mu = %.0e$', mu_values(i)); % Replace 'name' with your property name

        else
            name = sprintf('$\\mu = %.0f$', mu_values(i)); % Replace 'name' with your property name
            % name = sprintf('$\\theta_2 = %.0f$', mu_values(i)); % Replace 'name' with your property name


        end
    elseif mu_values(i) == inf
        name = "Synchronization";
    % elseif mu_values(i) == -inf
    %     name = 'Adaptive $\mu$';
    % else
    %     name = sprintf('$\\mu = %.1f$', mu_values(i)); % Replace 'name' with your property name
    end
    figure(plot_fig1);
    h1(i) = plot(NaN, NaN, 'Color', colors(i,:), 'DisplayName', name, 'LineWidth',2);
    figure(plot_fig2);
    h2(i) = plot(NaN, NaN, 'Color', colors(i,:), 'DisplayName', name,'LineWidth',2);
    figure(plot_fig3);
    h3(i) = plot(NaN, NaN, 'Color', colors(i,:), 'DisplayName', name,'LineWidth',2);

end

figure(plot_fig1);
legend(h1,'Location', 'northeast');

figure(plot_fig2);
legend(h2,'Location', 'northeast');

figure(plot_fig3);
legend(h3,'Location', 'northeast');





% k = 0;
% i = 1;
% for i = 1 : length(mu_values)
% mu = 0;
% files_for_mu = mu_files(mu);
% colors = jet(length(files_for_mu));
for i = 1 : length(mu_values)
    % mu_values(i)
    % files_for_ens = ens_files(ens_values(i));
    files_for_mu = mu_files(mu_values(i));
    % files_for_params = params_files(mu_values(i));
    color = colors(i,:);
    k = 1;
    for j = 1:length(files_for_mu)
    % for j = k:k
        file_name = files_for_mu{j};
        % param_name = files_for_params{j};

        %     k = mod(k,length(theFiles)) + 1;
        %     missing_mu = true;
        %     while(missing_mu)
        %         baseFileName = theFiles(k).name;
        %         fullFileName = fullfile(theFiles(k).folder, baseFileName);
        %     fprintf(1, 'Now reading %s\n', fullFileName);
        % Now do whatever you want with this file name,
        % such as calling load() to load the data.
        data = load(file_name);
        % data_params = load(param_name);

        % if (data.aot.T == 25)
            % && data.aot.IC_type == "Projection")



            % if isfield(data, 'var') && isprop(data.aot, 'error') % checking if the 'var' and 'error' fields exist


            try 
                if data.aot.assimilate_subtype == "Adaptive"
                    data.aot.assimlate_type ="Adaptive";
                end

            catch ME

            end

                % Compute offset for this series
                offset = mod(j, baseOffset)+1;

                % Compute indices of points to be marked
                markedIndices = offset:markerSpacing:length(data.var.error);
                
                % data.aot.error_low(1)
                % if data.aot.error_high(1) < 1e-10
                %     data.aot.error_high(1) = NaN;
                % % end
                % 
                % % if data.aot.error(1) < 1e-10
                %     data.aot.error(1) = NaN;
                % % end
                % 
                % % if data.aot.error_low(1) < 1e-10
                %     data.aot.error_low(1) = NaN;
                % end



                % data.aot.error_high = max(tol,data.aot.error_high);
                % data.aot.error_low = max(tol,data.aot.error_low);
                % data.aot.error = max(tol,data.aot.error);

                % 
                % figure(plot_fig1);
                % dt = p.dt;
                % rate_low = (log(data.var.error_low(2:end)) - log(data.var.error_low(1:end-1)))/dt;
                % rate_high= (log(data.var.error_high(2:end)) - log(data.var.error_high(1:end-1)))/dt;
                % rate = (log(data.var.error(2:end)) - log(data.var.error(1:end-1)))/dt;
                % rate(1) = nan;
                % rate_low(1) = nan;
                % rate_low(data.var.error_low(2:end)==eps) = NaN;
                % rate_high(1) = nan;
                % 
                % p1 = semilogy(p.t,-rate_low,'LineWidth',2,'Color',color); % plot error value
                % 
                % figure(plot_fig2);
                % 
                % p2 = plot(p.t,rate_high,'LineWidth',2,'Color',color); % plot error value
                % 
                % figure(plot_fig3);
                % 
                % p3 = plot(p.t,rate,'LineWidth',2,'Color',color); % plot error value


                figure(plot_fig1);
                dt = p.dt;

                data.var.error_low(data.var.error_low == eps) = NaN;
                data.var.error_low(data.var.error_low <= eps) = eps;
                data.var.error_high(data.var.error_high <= eps) = eps;
                data.var.error(data.var.error <= eps) = eps;

                p1 = semilogy(p.t,(data.var.error_low(2:end)),'LineWidth',2,'Color',color); % plot error value

                figure(plot_fig2);

                p2 = semilogy(p.t,(data.var.error_high(2:end)),'LineWidth',2,'Color',color);
                figure(plot_fig3);

                p3 = semilogy(p.t,(data.var.error(2:end)),'LineWidth',2,'Color',color);



                % Set 'MarkerIndices' property of the plot
                p1.MarkerIndices = markedIndices;
                p2.MarkerIndices = markedIndices;
                p3.MarkerIndices = markedIndices;

                % Set 'Marker' property of the plot
                p1.Marker = 'o';
                p2.Marker = 'o';
                p3.Marker = 'o';

                hold on;
                % switch data.aot.assimilate_type
                %     case "AOT"
                %         if fix(data.aot.mu)|| data.aot.mu == 0
                %             legendLabels{i} = sprintf('$\\mu = %d$', data.aot.mu); % Replace 'name' with your property name
                %         else
                %             legendLabels{i} = sprintf('$\\mu = %.2f$', data.aot.mu); % Replace 'name' with your property name
                %         end
                %         % Set 'Marker' property of the plot
                %         p1.Marker = 'o';
                %     case "Synchronization"
                %         legendLabels{i} = 'Synchronization';
                %         % Set 'Marker' property of the plot
                %         p1.Marker = '+';
                %         p2.Marker = '+';
                %         p3.Marker = '+';
                %     case "Adaptive"
                %         legendLabels{i} = 'Adaptive';
                %         % Set 'Marker' property of the plot
                %         p1.Marker = '+';
                %         p2.Marker = '+';
                %         p3.Marker = '+';
                % end
            % else
                % fprintf('The file %s does not contain var.error field\n', baseFileName);
            % end
            missing_mu = false;

        % end
        k = mod(k,length(theFiles)) + 1;

    end
end

figure(plot_fig1);
% set(gca, 'YScale', 'log');
% set(gca, 'YScale', 'linear');
% legend(legendLabels,'Interpreter', 'latex');
legend(h1,'Interpreter', 'latex','FontSize', 30,'Location', 'northeast');
% title(plotTitle,'FontSize', 40);
xlabel('Time','FontSize', 40);
% ylabel('Decay rate constant for observed error','FontSize', 40)
ylabel('L^2 Error on observed modes','FontSize', 40)

% axis[];
hold off; % Releases hold on the plot

% Get current y-axis limits
% axis([0 25 1e-18 1e+2]);
% axis([0 1 1e-16 1e+2]);

% Set new y-axis limits, preserving the max and setting the min to 1e-15
% ylim([1e-18, 1e+2]);
% ylim([-17, 1])
% ylim([-4, 1])


ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2; % Set the width of the axes lines, which also affects the grid lines
  % Force major ticks on 1e-n marks
                ax.YScale = 'log';
                ax.YTick = 10.^(-30:2);
                % ax.XTick = 0:1:time_adj(end);
                ax.YGrid = 'on';
                ax.XGrid = 'off';
                ax.MinorGridLineStyle = 'none'; % Optional: Turn off minor grid lines
                ax.YMinorTick = 'off';
                ax.GridLineWidth = 2;

figure(plot_fig2);
% set(gca, 'YScale', 'log');
legend(h2,'Interpreter', 'latex','FontSize', 30,'Location', 'northeast');
% title(plotTitle,'FontSize', 40);
xlabel('Time','FontSize', 40);
% ylabel('Decay rate constant for unobserved error','FontSize', 40)
ylabel('L^2 Error on unobserved modes','FontSize', 40)
% axis[];
hold off; % Releases hold on the plot

% axis([0 25 1e-18 1e+2]);
% axis([0 1 1e-6 1e+2]);

% Get current y-axis limits
y_limits = ylim;

% Set new y-axis limits, preserving the max and setting the min to 1e-15
% ylim([1e-18, 1e+2]);
% ylim([-4, 1])
% ylim([-4, 1])


ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2; % Set the width of the axes lines, which also affects the grid lines
  % Force major ticks on 1e-n marks
                ax.YScale = 'log';
                ax.YTick = 10.^(-30:2);
                ax.YGrid = 'on';
                ax.XGrid = 'off';
                ax.MinorGridLineStyle = 'none'; % Optional: Turn off minor grid lines
                ax.YMinorTick = 'off';
                ax.GridLineWidth = 2;

figure(plot_fig3);
% set(gca, 'YScale', 'log');
legend(h3,'Interpreter', 'latex','FontSize', 30,'Location', 'northeast');
% title(plotTitle,'FontSize', 40);
xlabel('Time','FontSize', 40);
% ylabel('Decay rate constant for L^2 error','FontSize', 40)
ylabel('L^2 Error','FontSize', 40)
% axis([0 25 1e-18 1e+2]);
% axis([0 1 1e-6 1e+2]);
% axis[];
hold off; % Releases hold on the plot
% Get current y-axis limits
% y_limits = ylim;
% ylim([-17, 1])
% ylim([-4, 1])


% Set new y-axis limits, preserving the max and setting the min to 1e-15
% ylim([1e-18, 1e+2]);

ax = gca;
ax.FontSize = 20;
ax.LineWidth = 2; % Set the width of the axes lines, which also affects the grid lines

  % Force major ticks on 1e-n marks
                ax.YScale = 'log';
                ax.YTick = 10.^(-30:2);
                % ax.XTick = 0:1:time_adj(end);
                ax.YGrid = 'on';
                ax.XGrid = 'off';
                ax.MinorGridLineStyle = 'none'; % Optional: Turn off minor grid lines
                ax.YMinorTick = 'off';
                ax.GridLineWidth = 2;
                % 

% Get the screen size
    screenSize = get(0, 'ScreenSize');

    % Set the figure position to the screen size to maximize it
    set(plot_fig1, 'Position', [screenSize(1), screenSize(2), screenSize(3), screenSize(4)]);
    set(plot_fig2, 'Position', [screenSize(1), screenSize(2), screenSize(3), screenSize(4)]);
    set(plot_fig3, 'Position', [screenSize(1), screenSize(2), screenSize(3), screenSize(4)]);

% % % 
saveas(plot_fig1, 'error_primary_zero_low.jpg');
saveas(plot_fig2, 'error_primary_zero_high.jpg');
saveas(plot_fig3, 'error_primary_zero_total.jpg');

saveas(plot_fig1, 'error_primary_zero_low.fig');
saveas(plot_fig2, 'error_primary_zero_high.fig');
saveas(plot_fig3, 'error_primary_zero_total.fig');