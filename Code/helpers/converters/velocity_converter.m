function [psi] = velocity_converter(u1,u2,p)
% Converts velocity in physical space to stream function psi
%(given in Fourier space)
% curl(u1,u2) = -Lap psi
u1_hat = fft2(u1);
u2_hat = fft2(u2);
u2_x_hat = p.IKX.*u2_hat;
u1_y_hat = p.IKY.*u1_hat;
% u2_x_hat = bsxfun(@times,p.ikx,u2_hat);
% u1_y_hat = bsxfun(@times,p.iky,u1_hat);
psi =  -p.k_lap_inv.*(u2_x_hat - u1_y_hat);
end