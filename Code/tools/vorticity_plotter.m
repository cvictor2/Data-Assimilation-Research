function vorticity_plotter(psi_hat, p)

omega = ifftn( -p.k_lap.*psi_hat,'symmetric');
x = p.x;
y = p.y;
%             t = p.t;
            
figname = "Vorticity plotted from data";
% figname = sprintf('Observer Type: %s, Smoothing Type: %s', var.observer_type,  var.assimilate_type);
varfigure = figure('Position',[1 1 1100 500],'Color',[1 1 1], 'Name', figname);
            %                 subplot(2,2,1);
            %                         varfigures(i) =  tightfig(varfigures(i));
            %
pcolor(x,y,omega);
axis('square');
axis tight;
colormap jet; % Bone, copper, jet
shading interp; % flat, interp, faceted
lighting phong;
colorbar;
%         caxis([-2,2]);
%         axmin = max(max(max(-p.k_lap.*psi_hat)),1);
%         caxis([-axmin axmin]);

%         caxis([-.25,.25]);
% title(sprintf('Reference Vorticity at t = %1.2f',0));
title("Vorticity");

end