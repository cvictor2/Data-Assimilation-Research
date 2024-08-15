function saveRoutine(settings, parameters, plots, ref_soln, vars)
    runDir = fullfile('runs', [settings.save_name + "_incomplete"]);
    settingsPath = fullfile(runDir, 'settings.mat');
    paramsPath = fullfile(runDir,'parameters.mat');
    refPath = fullfile(runDir,'ref_soln.mat');
    daPath = fullfile(runDir,'DA');
    plotDir = fullfile(runDir,'figures');

    save(settingsPath,'settings');
    save(paramsPath,'parameters');
    
    if(settings.save_ref)
        save(refPath, 'ref_soln');
    end
    if(settings.save_var)
        for i = 1: length(vars)
            daFile = fullfile(daPath,["DA" + num2str(i) + ".mat"]);
            var = vars(i);
            ensemble_temp = var.ensemble;
            var.ensemble = [];
            save(daFile,'var');
            var.ensemble = ensemble_temp;
        end
    end
    if settings.save_plots
        plots.save(plotDir);

    end

end
