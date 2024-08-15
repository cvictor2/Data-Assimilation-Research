classdef DA_plotter < handle
    properties
        var
        figure_window
        plots_array
        settings
        figname


    end
    methods
        function obj = DA_plotter(p, settings, psi_hat, var)

            if nargin == 0
                return
            end

            % function var = initialize_plots(var, p, psi_hat)
            omega = ifftn( -p.k_lap.*psi_hat,'symmetric');
            omega_v = ifftn(-p.k_lap.*var.v_hat,'symmetric');

            x = p.x;
            y = p.y;
            t = p.t;

            obj.var = var;

            obj.figname = sprintf('Observer Type: %s, Smoothing Type: %s', var.observer_type,  var.assimilate_type);
            obj.figure_window = figure('Position',[1 1 1100 500],'Color',[1 1 1], 'Name', obj.figname);
            %                 subplot(2,2,1);
            %                         varfigures(i) =  tightfig(varfigures(i));
            %
            subplot(2,2,1);
            obj.plots_array(1) = pcolor(x,y,omega);
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
            title(sprintf('Reference Vorticity at t = %1.2f',0));

            %                 title(sprintf('Vorticity at t = %1.2f',t(1)));
            %     get(p1,'cdata');
            set(obj.plots_array(1),'cdata',ifftn( -p.k_lap.*psi_hat,'symmetric'));
            %         drawnow;
            hold on;
            %             if(options(7))
            %                 obj.plots_array(5) = scatter(var.nodes_coordinates(:,2),var.nodes_coordinates(:,1),10, 'filled');
            %             end


            subplot(2,2,2);
            %         intdiff = DA_interp(psi_hat, var, p);

            %         intdiff = intdiff/(norm(intdiff(:))*sqrt(Lx*Ly)/Nx/Ny);
            obj.plots_array(2) = pcolor(x,y,omega);
            set(obj.plots_array(2),'cdata',ifftn( -p.k_lap.*var.u_obs,'symmetric'));
            axis('square');
            %         colormap jet;
            colormap jet; % Bone, copper, jet
            shading interp; % flat, interp, faceted
            lighting phong;
            colorbar;
            %         caxis([-10,10]);

            title(sprintf('Interpolated Vorticity at t = %1.2f',0));

            %         hold on;
            %         varplots(i,5) = scatter(var.nodes_coordinates(:,2),var.nodes_coordinates(:,1),5,'filled');
            %


            %         var.i_nodes_x;
            %         scatter3(X,Y,Z);
            %         stem3(X,Y,V,'.','color','k','MarkerSize',15)

            omega_umv = ifftn(-p.k_lap.*(obj.var.v_hat - psi_hat),'symmetric');
            subplot(2,2,3);
            obj.plots_array(3) = pcolor(x,y,abs(omega_umv));
            hold on;
            axis('square');
            % colormap copper; % Bone, copper, jet
            shading interp; % flat, interp, faceted
            lighting phong;
            %         caxis([-3,3]);
            title(sprintf('Absolute Vorticity Difference at t = %1.2f',p.time_current));
            clim([1e-5 1e+0]);

            % colormap(hot);
            colorbar;
            ax = gca;
            ax.ColorScale = 'log';
            %         caxis([-.25,.25]);

            % obj.plots_array(3) = semilogy(0:p.dt:p.time_final+p.dt,var.error,'--','LineWidth',1);
            % hold on;
            % semilogy(0:p.dt:p.time_final+p.dt,var.error_low,'--','LineWidth',1);
            % semilogy(0:p.dt:p.time_final+p.dt,var.error_high,'--','LineWidth',1);
            % hold off;
            %         % varplots(i,3) = semilogy(0:dt:T,var.ens_umv,'--','LineWidth',1);
            % title('Error for Simulated Solution');
            % xlabel('Time');
            % axis([0 p.time_final 1e-20 1e-0]);

            %         ylabel('Enstrophy of u-v');
            %         axis([0 T 1e-20 1e-0]);
            %         hold on;
            %         text = sprintf('vars(%d).error',e);
            %         set(errorPlots(e), 'YDataSource',text);
            %         varplots(i,3) = pcolor(x,y,var.v - ifftn( -p.k_lap.*psi_hat,'symmetric'));
            %         hold on;
            %         axis('square');
            %         % colormap copper; % Bone, copper, jet
            %         shading interp; % flat, interp, faceted
            %         lighting phong;
            %         %         caxis([-3,3]);
            % %         title(sprintf('Vorticity Differences at t = %1.2f',t(1)));
            %         get(varplots(i,3),'cdata');
            %         caxis([0,.25]);
            %         title(sprintf('Absolute Vorticity Difference at t = %1.2f',0));


            %         drawnow;
            subplot(2,2,4);
            obj.plots_array(4) = pcolor(x,y,omega_v);
            hold on;
            axis('square');
            % colormap copper; % Bone, copper, jet
            shading interp; % flat, interp, faceted
            lighting phong;
            %         caxis([-3,3]);
            title(sprintf('Vorticity at t = %1.2f',t(1)));
            get(obj.plots_array(4),'cdata');
            colorbar;
            %         caxis([-.25,.25]);
            title(sprintf('Simulated Vorticity at t = %1.2f',0));

            drawnow;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Update Routine %%

        function obj = update_figures(obj, p, settings, psi_hat)
            figure(obj.figure_window);
            %             figure(fh1);
            subplot(2,2,1);
            set(obj.plots_array(1),'cdata',ifftn( -p.k_lap.*psi_hat,'symmetric'));
            title(sprintf('Reference Vorticity at t = %1.2f',p.time_current));
            colorbar;
            %             axmin = max(max(max(psi)),0.25);
            %             caxis([-axmin axmin]);
            %             drawnow;
            %             if(options(7))
            %                 set(varplot(5),'ydata',var.nodes_coordinates(:,1));
            %                 set(varplots(i,5),'xdata',var.nodes_coordinates(:,2));
            %             end

            subplot(2,2,2);
            %             intdiff = var.mu*DA_interp(psi_hat, var,p);
            %             intdiff/(norm(intdiff(:))*sqrt(Lx*Ly)/Nx/Ny);
            %             axmax = max(max(abs(ifftn( -p.k_lap.*psi_hat,'symmetric'))));
            set(obj.plots_array(2),'cdata',ifftn( -p.k_lap.*obj.var.u_obs,'symmetric'));
            title(sprintf('Interpolated Vorticity at t = %1.2f',p.time_current));
            colorbar;
            %             caxis([-axmax,axmax]);


            %             axmin = max(max(max(psi)),0.25);
            %             caxis([-axmin axmin]);
            %             drawnow;


            subplot(2,2,3);
            set(obj.plots_array(3),'cdata', ( abs(ifftn(-p.k_lap.*(obj.var.v_hat - psi_hat),'symmetric')) ));
            % set(obj.plots_array(3),'cdata',( abs(ifftn(-p.k_lap.*fftn(obj.var.v_hat), 'symmetric' ) - ifftn( -p.k_lap.*fftn(psi_hat),'symmetric')) ));
            % % %              set(varplots(i,3),'cdata',( ifftn(-p.k_lap.*fftn(var.v), 'symmetric' )));
            title(sprintf('Absolute Vorticity Difference at t = %1.2f',p.time_current));
            colorbar;
            % %                 caxis([-.5,.5]);
            % %              hold on;
            % %             varplots(i,3) = semilogy(0:dt:T,var.error,'-','LineWidth',1);
            % %             varplots(i,3) = semilogy(0:dt:T,var.ens_umv,'-','LineWidth',1);
            % semilogy(0:p.dt:p.time_final+p.dt,obj.var.error,'-','LineWidth',1);
            % hold on;
            % semilogy(0:p.dt:p.time_final+p.dt,obj.var.error_low,'--','LineWidth',1);
            % semilogy(0:p.dt:p.time_final+p.dt,obj.var.error_high,'--','LineWidth',1);
            % hold off;
            %
            % %             hold off;
            % %                         legend('Enstrophy Error', 'Psi L^2 Error');
            % legend('Psi L^2 Error (total)', 'Psi L^2 Error (low modes)', 'Psi L^2 Error (high modes)');
            %
            % title('Error of Simulated Solution');
            % xlabel('Time');
            % ylabel('Error');
            % axis([0 p.time_final 1e-20 1e-0]);



            subplot(2,2,4);
            set(obj.plots_array(4),'cdata',( ifftn(-p.k_lap.*obj.var.v_hat, 'symmetric' )));
            title(sprintf('Simulated Vorticity at t = %1.2f',p.time_current));
            colorbar;
            %                             caxis([-30,30]);
            %              hold on;


            drawnow;
            %
            %
        end


    end
end