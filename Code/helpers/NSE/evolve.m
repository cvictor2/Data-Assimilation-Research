function [energy,u_next] = evolve(p, ev, u_prev)

%Forward Euler
[energy,k1] = compute_rhs_hat(p.time_current, u_prev, p);
u_next = ev.E2.*u_prev + ev.dt*ev.E2.*k1;
end

