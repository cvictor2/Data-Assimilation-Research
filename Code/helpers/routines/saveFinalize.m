
function saveFinalize(settings, parameters, plots, ref_soln, vars)
saveRoutine(settings, parameters, plots, ref_soln, vars);

runDir_old = fullfile('runs', [settings.save_name + "_incomplete"]);
runDir_new = fullfile('runs', [settings.save_name]);
movefile(runDir_old, runDir_new);

end
