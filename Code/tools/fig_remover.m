close all;

% Specify the folder where the files live.
myFolder = 'C:\Users\Collin\Desktop\AI\DA_code-main\Data\noise'; % Replace with your folder path
saveFolder = 'C:\Users\Collin\Desktop\AI\DA_code-main\Data\noise\processed';

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

if ~isfolder(saveFolder)
    mkdir(saveFolder);
end


% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*_Hybrid_*complete.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

% plot_fig = figure;

% hold on; % Holds the plot so all data can be plotted on same graph


% legendLabels = {};

for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(theFiles(k).folder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Now do whatever you want with this file name,
  % such as calling load() to load the data.
  data = load(fullFileName);

    % Remove the figure handle if it exists
  if isfield(data, 'var') % Replace 'h' with the name of your figure handle if it's different
%       if isfield(data.var, 'varfigure')
          data.var.varfigure = [];
%       end
  end
  newFilePath = fullfile(saveFolder, baseFileName);
  save(newFilePath, '-struct', 'data');



end
close all;