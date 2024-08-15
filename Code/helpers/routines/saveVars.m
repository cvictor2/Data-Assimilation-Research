function saveVars(vars, p, message)
for i = 1:p.size_vars
    date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS'); %#ok<*TNOW1,*DATST>
    var_type = 'Misc.'; %#ok<NASGU>
    var = vars(i);
    if ispc
        save(sprintf('Data/%s_vars_%s_%s_%s_%s.mat',date_string,var.observer_type,var.assimilate_type,p.note,message),'var', 'p');
    else
        save(sprintf('Data\\%s_vars_%s_%s_%s_%s.mat',date_string,var.observer_type,var.assimilate_type,p.note,message),'var', 'p');
    end
    pause(1);
end
end

