function energy = generate_energy(p, psi_hat)

% psi_hat(p.dealias_modes) = 0;

% enstrophy_sq = u_x_hat*u_x_hat'/p.Nsq;

% psi_x_hat = bsxfun(@times,p.ikx,psi_hat);
% psi_y_hat = bsxfun(@times,p.iky,psi_hat);

psi_x_hat = p.IKX.*psi_hat;
psi_y_hat = p.IKY.*psi_hat;
energy = sqrt(norm(psi_x_hat,'fro')^2 + norm(psi_y_hat,'fro')^2)/p.Nx/p.Ny;


% sqrt(norm(u1,'fro')^2 + norm(u2,'fro')^2)/p.Nx/p.Ny;

end