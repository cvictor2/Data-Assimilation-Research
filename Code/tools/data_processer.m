close all;

% Specify the folder where the files live.
myFolder = 'C:\Users\Collin\Desktop\AI\DA_code-main\Data\noise\processed'; % Replace with your folder path

plotTitle = 'Error for High Number of Observed Modes';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*AOT_*_complete.mat'); % Change to whatever pattern you need.
theFiles_AOT = dir(filePattern);

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*Synchronization_*_complete.mat'); % Change to whatever pattern you need.
theFiles_Synch = dir(filePattern);

filePattern = fullfile(myFolder, '*Hybrid_*_complete.mat'); % Change to whatever pattern you need.
theFiles_Hybrid = dir(filePattern);

theFiles = [theFiles_AOT; theFiles_Synch; theFiles_Hybrid];

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

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as calling load() to load the data.
    data = load(fullFileName);

    if data.var.assimilate_type == "Synchronization"
        mu = inf;
    elseif data.var.assimilate_type == "Hybrid" 
        mu = -inf;
    else
        mu = data.var.mu;
    end

    if isKey(mu_files,mu)
        mu_files(mu) = [mu_files(mu), {fullFileName}];
    else
        mu_files(mu) = {fullFileName};
    end


    %    missing_mu = ~ismember(mu, mu_values);
    %     if missing_mu
    %         mu_values = [mu_values, mu];
    %     end
end

mu_values = cell2mat(keys(mu_files));

mu_values = sort(mu_values);



% Specify specific mu values for plotting
% target_mus = [1,1000,1000000,1000000000,3];
target_mus = mu_values;

colors = jet(length(target_mus)); % Generates unique colors

% Define the names for the 4 archetypes
% names = {'\mu = 1', '\mu = 1000', '\mu = 1000000', '\mu = 1000000000','Synchornization'};



% Set up coloring for all figures for each mu value

% Initialize the handles array for the 4 archetypes
h1 = gobjects(length(target_mus),1);
h2 = gobjects(length(target_mus),1);
h3 = gobjects(length(target_mus),1);

mu_values = target_mus;


% figure(plot_fig1);
% hold on;
for i = 1: length(mu_values)
    if fix(mu_values(i))|| mu_values(i) == 0
        name = sprintf('$\\mu = %d$', mu_values(i)); % Replace 'name' with your property name
    elseif mu_values(i) == inf
        name = "Synchronization";
    elseif mu_values(i) == -inf
        name = "Adaptive $\\mu";
    else
        name = sprintf('$\\mu = %.2f$', mu_values(i)); % Replace 'name' with your property name
    end
    figure(plot_fig1);
    h1(i) = plot(NaN, NaN, 'Color', colors(i,:), 'DisplayName', name, 'LineWidth',2);
    figure(plot_fig2);
    h2(i) = plot(NaN, NaN, 'Color', colors(i,:), 'DisplayName', name,'LineWidth',2);
    figure(plot_fig3);
    h3(i) = plot(NaN, NaN, 'Color', colors(i,:), 'DisplayName', name,'LineWidth',2);

end

figure(plot_fig1);
legend(h1);

figure(plot_fig2);
legend(h2);

figure(plot_fig3);
legend(h3);





% k = 0;
% i = 1;
% for i = 1 : length(mu_values)
% mu = 0;
% files_for_mu = mu_files(mu);
% colors = jet(length(files_for_mu));
for i = 1 : length(mu_values)
    files_for_mu = mu_files(mu_values(i));
    color = colors(i,:);
    for j = 1: length(files_for_mu)

        file_name = files_for_mu{j};

        %     k = mod(k,length(theFiles)) + 1;
        %     missing_mu = true;
        %     while(missing_mu)
        %         baseFileName = theFiles(k).name;
        %         fullFileName = fullfile(theFiles(k).folder, baseFileName);
        %     fprintf(1, 'Now reading %s\n', fullFileName);
        % Now do whatever you want with this file name,
        % such as calling load() to load the data.
        data = load(file_name);

        if (data.p.T == 25 && data.var.IC_type == "Projection")



            if isfield(data, 'var') && isprop(data.var, 'error') % checking if the 'var' and 'error' fields exist




                % Compute offset for this series
                offset = mod(j, baseOffset)+1;

                % Compute indices of points to be marked
                markedIndices = offset:markerSpacing:length(data.p.t);

                data.var.error_high = max(tol,data.var.error_high);
                data.var.error_low = max(tol,data.var.error_low);
                data.var.error = max(tol,data.var.error);


                figure(plot_fig1);

                p1 = semilogy(data.p.t,data.var.error_low(1:end-1),'LineWidth',2,'Color',color); % plot error value

                figure(plot_fig2);

                p2 = semilogy(data.p.t,data.var.error_high(1:end-1),'LineWidth',2,'Color',color); % plot error value

                figure(plot_fig3);

                p3 = semilogy(data.p.t,data.var.error(1:end-1),'LineWidth',2,'Color',color); % plot error value


                % Set 'MarkerIndices' property of the plot
                p1.MarkerIndices = markedIndices;
                p2.MarkerIndices = markedIndices;
                p3.MarkerIndices = markedIndices;

                % Set 'Marker' property of the plot
                p1.Marker = 'o';
                p2.Marker = 'o';
                p3.Marker = 'o';

                hold on;
                switch data.var.assimilate_type
                    case "AOT"
                        if fix(data.var.mu)|| data.var.mu == 0
                            legendLabels{i} = sprintf('$\\mu = %d$', data.var.mu); % Replace 'name' with your property name
                        else
                            legendLabels{i} = sprintf('$\\mu = %.2f$', data.var.mu); % Replace 'name' with your property name
                        end
                        % Set 'Marker' property of the plot
                        p1.Marker = 'o';
                    case "Synchronization"
                        legendLabels{i} = 'Synchronization';
                        % Set 'Marker' property of the plot
                        p1.Marker = '+';
                        p2.Marker = '+';
                        p3.Marker = '+';
                    case "Hybrid"
                        legendLabels{i} = 'Hybrid';
                        % Set 'Marker' property of the plot
                        p1.Marker = '+';
                        p2.Marker = '+';
                        p3.Marker = '+';
                end
            else
                fprintf('The file %s does not contain var.error field\n', baseFileName);
            end
            missing_mu = false;

        end
        k = mod(k,length(theFiles)) + 1;

    end
end

figure(plot_fig1);
set(gca, 'YScale', 'log');
% legend(legendLabels,'Interpreter', 'latex');
legend(h1,'Interpreter', 'latex');
title(plotTitle);
xlabel('Time');
ylabel('L^2 norm of u-v on observed modes')
% axis[];
hold off; % Releases hold on the plot

% Get current y-axis limits
y_limits = ylim;

% Set new y-axis limits, preserving the max and setting the min to 1e-15
% ylim([1e-15, y_limits(2)]);

figure(plot_fig2);
set(gca, 'YScale', 'log');
legend(h2,'Interpreter', 'latex');
title(plotTitle);
xlabel('Time');
ylabel('L^2 norm of u-v on unobserved modes')
% axis[];
hold off; % Releases hold on the plot

% Get current y-axis limits
y_limits = ylim;

% Set new y-axis limits, preserving the max and setting the min to 1e-15
% ylim([1e-15, y_limits(2)]);

figure(plot_fig3);
set(gca, 'YScale', 'log');
legend(h3,'Interpreter', 'latex');
title(plotTitle);
xlabel('Time');
ylabel('L^2 norm of u-v on all modes')
% axis[];
hold off; % Releases hold on the plot

% Get current y-axis limits
y_limits = ylim;

% Set new y-axis limits, preserving the max and setting the min to 1e-15
% ylim([1e-15, y_limits(2)]);
