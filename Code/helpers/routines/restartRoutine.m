
function [parameters, settings, ref_soln, vars] = restartRoutine(runFolder, settings)
% runDir = fullfile('runs', [settings.save_name, '_incomplete']);
if nargin== 0
    
    runDir = getMostRecentDir();
settingsPath = fullfile(runDir, 'settings.mat');
    load(settingsPath,'settings');
elseif nargin == 2
    % settingsPath = fullfile(runDir, 'settings.mat');
    runDir = fullfile('runs', runFolder);

else

    runDir = fullfile('runs', runFolder);
    settingsPath = fullfile(runDir, 'settings.mat');
    load(settingsPath,'settings');

end
if ~exist(runDir,'dir')
    runDir  = strcat(runDir, "_incomplete");
end
paramsPath = fullfile(runDir,'parameters.mat');
refPath = fullfile(runDir,'ref_soln.mat');
daPath = fullfile(runDir,'DA');

load(paramsPath,'parameters');
load(refPath,'ref_soln');


vars = repelem(DA_obs(parameters),1,1);
daFiles = dir(fullfile(daPath, 'DA*.mat'));

numFiles = numel(daFiles);

for i = 1:numFiles
    filePath = fullfile(daPath, daFiles(i).name);
    tempData = load(filePath);
    vars(i) = tempData.var;
end



%Edit and extend all of the vars and whatnot
if parameters.time_current == parameters.time_final || parameters.time_final ~= settings.final_time
    [parameters, settings, vars] = restartExtend(parameters,settings,vars);
end

end

