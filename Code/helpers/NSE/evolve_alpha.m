
function u_next = evolve_alpha(p, ev, u_prev) %#ok<*DEFNU>

%Forward Euler
[~,k1] = compute_rhs_hat(p.time, u_prev, p);
u_next = ev.E2_alpha.*u_prev + ev.dt*ev.E2.*k1;
end

