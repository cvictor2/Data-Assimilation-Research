function [energy,rhs_hat]= compute_rhs_hat(t,psi_hat,p) %#ok<*INUSL>
% Compute the fft of the right-hand side of the equation.

% Dealias by setting the "middle" Fourier coefficients to zero.
% This is where Matlab stores the highest frequencies.
psi_hat(p.dealias_modes) = 0;

% enstrophy_sq = u_x_hat*u_x_hat'/p.Nsq;

% psi_x_hat = bsxfun(@times,p.ikx,psi_hat);
% psi_y_hat = bsxfun(@times,p.iky,psi_hat);
psi_x_hat = p.IKX.*psi_hat;
psi_y_hat = p.IKY.*psi_hat;
energy = (norm(psi_x_hat)^2 + norm(psi_y_hat)^2)/p.Nx/p.Ny;
%energy = 0;

% % Compute the nonlinear terms in physical space using dealiased versions.
% psi_x      = ifftn(              psi_x_hat,'symmetric'); % psi_x = -v
% psi_y      = ifftn(              psi_y_hat,'symmetric'); % psi_y =  u
% psi_x_lap  = ifftn(p.k_lap.*psi_x_hat,'symmetric');
% psi_y_lap  = ifftn(p.k_lap.*psi_y_hat,'symmetric');
% % rhs_hat = p.helm_inv.*p.k_lap_inv.*fft2(-psi_y.*psi_x_lap + psi_x.*psi_y_lap + p.f_hat);
% rhs_hat = p.k_lap_inv.*(fft2(-psi_y.*psi_x_lap + psi_x.*psi_y_lap) + p.f_hat);

% Basdevant formula (4 fft/ifft operations, 6 products). See:
%   P. Emami and J. C. Bowman,"On the Global Attractor of 2D Incompressible Turbulence with Random Forcing",
%   Journal of Differential Equations. 264, 4036-4066 (2018).
u1 =  ifftn(psi_y_hat,'symmetric');
u2 = -ifftn(psi_x_hat,'symmetric');
% u1       =  ifftn(bsxfun(@times,p.iky,psi_hat),'symmetric'); % psi_y =  u
% u2       = -ifftn(bsxfun(@times,p.ikx,psi_hat),'symmetric'); % psi_x = -v
rhs_hat = p.k_lap_inv.*(p.k1sq_m_k2sq.*fftn(u1.*u2) + p.k1k2.*fftn(u2.^2 - u1.^2) + p.f_hat);
% rhs_hat = p.helm_inv.*p.k_lap_inv.*(p.k1sq_m_k2sq.*fftn(u1.*u2) + p.k1k2.*fftn(u2.^2 - u1.^2) + p.f_hat);

% if p.compute_norms == 1
%     % Compute expensive norms while we have the data.
%     p.compute_norms = 0;
%     u_Linf = sqrt(max(u1(:).^2 + (u2(:)).^2));
% else
%     u_Linf = 0;
% end


end % =========== End function rhs_hat ============

