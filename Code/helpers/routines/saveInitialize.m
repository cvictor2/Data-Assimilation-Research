
function saveInitialize(settings, params, plots, ref_soln, vars)
    runDir_final = fullfile('runs', settings.save_name);
    runDir = fullfile('runs', [settings.save_name + "_incomplete"]);
    if exist(runDir_final, 'dir') %If restarting there may not be an _incomplete tag where it should be saved, so it should be added.
        movefile(runDir_final, runDir)
    elseif ~exist( runDir, 'dir')
        mkdir(runDir);
    end

    if settings.save_var
        daDir = fullfile(runDir, 'DA');
        if ~exist(daDir, 'dir')
            mkdir(daDir);
        end
    end

    if settings.save_plots
        plotDir = fullfile(runDir,'figures');
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end


    end
    saveRoutine(settings, params, plots, ref_soln, vars);

end

