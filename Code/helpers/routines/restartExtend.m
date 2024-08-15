
function [p, s, vars] = restartExtend(p, s, vars)
% This function extends the matrices that are based on time when a restart
% is called

% %If settings has a different final time extend p to match
% if s.time_final ~= p.time_final
if p.time_final == s.final_time
    s.final_time = s.final_time + s.restart_extend;
end
p = p.restart(s);
for i = 1: length(vars)
    vars(i) = vars(i).restart(p);
end
% %If settings and p agree, then we need to extend both by settings.restart_extend
% else


end

