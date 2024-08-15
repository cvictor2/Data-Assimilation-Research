function saveData(psi_hat_data, p)
date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
% save(sprintf('Data/%s_data_%s.mat',date_string,"restart"),'psi_hat', 'p');
if ispc
    save(sprintf('Data/%s_data_%s_%s.mat',date_string,p.note,"restart"),'psi_hat_data', 'p', '-v7.3');
else
    save(sprintf('Data\\%s_data_%s_%s.mat',date_string,p.note,"restart"),'psi_hat_data', 'p', '-v7.3');
end
end

