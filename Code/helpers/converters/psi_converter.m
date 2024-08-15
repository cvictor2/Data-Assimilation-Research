function [u1,u2] = psi_converter(psi_hat, p)
% Converts stream function psi to velocity components u1 and u2
% u1 =  ifftn(bsxfun(@times,p.iky,psi_hat),'symmetric'); % psi_y =  u
% u2 = -ifftn(bsxfun(@times,p.ikx,psi_hat),'symmetric'); % psi_x = -v
u1 = p.IKY.*psi_hat;
u2 = -p.IKX.*psi_hat;
u1 = ifftn(u1,'symmetric');
u2 = ifftn(u2,'symmetric');
end