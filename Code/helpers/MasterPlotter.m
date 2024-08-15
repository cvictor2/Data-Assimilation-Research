classdef MasterPlotter
    properties
        settings
        parameters

        da_plotters

        ref_figure
        ref_plot

        en_figure
        en_plots

        spec_figure
        spec_plots

        error_figure_low
        error_figure_high
        error_figure_total


        fontsize = 20;

        linewidth = 2;

    end
    methods
        function obj = MasterPlotter(settings, parameters, psi_hat, vars, energy)
            obj.settings = settings;
            obj.parameters = parameters;

            if isempty(vars)
                obj.settings.plot_vars = false;

            end




            if settings.plot_ref
                obj.ref_figure = figure();

                figure(obj.ref_figure)
                k_lap = parameters.k_lap;
                % function var = initialize_plots(var, p, psi_hat)
                omega = ifftn( -k_lap.*psi_hat,'symmetric');


                % figname = sprintf('Observer Type: %s, Smoothing Type: %s', var.observer_type,  var.assimilate_type);
                % obj.figure_window = figure('Position',[1 1 1100 500],'Color',[1 1 1], 'Name', obj.figname);
                %                 subplot(2,2,1);
                %                         varfigures(i) =  tightfig(varfigures(i));
                %
                % subplot(2,2,1);
                obj.ref_plot = pcolor(parameters.x,parameters.y,omega);
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
                % title(sprintf('Reference Vorticity at t = %1.2f',parameters.time_current));

                %                 title(sprintf('Vorticity at t = %1.2f',t(1)));
                %     get(p1,'cdata');
                set(obj.ref_plot,'cdata',ifftn( -k_lap.*psi_hat,'symmetric'));
                %         drawnow;
                % hold on;



            end

            if settings.plot_energy
                obj.en_figure = figure();


                figure(obj.en_figure);
                hold on;

                % title()
                % title("Approximate Slope of Error", "FontSize", obj.fontsize, "FontUnits",'points');
                % title("Kinetic Energy", "FontSize", obj.fontsize, "FontUnits",'points');
                xlabel("Time","FontSize", obj.fontsize);
                ylabel("Magnitude", "FontSize", obj.fontsize);
                set(gca, 'LineWidth', obj.linewidth);

                t = parameters.t;

                obj.en_plots = [];


                obj.en_plots(1) = plot(t, parameters.energy, "LineWidth",obj.linewidth);
                grid on;
                energy = parameters.energy;
                set(obj.en_plots(1), 'YDataSource','energy');
                % refreshdata(obj.en_plots(1));
                energy_vars = [];
                for i = 1:length(vars)
                    energy_vars(i,:) = vars(i).energy;
                    obj.en_plots(i+1) = plot(t, vars(i).energy, "LineWidth", obj.linewidth);
                    set(obj.en_plots(i+1), 'YDataSource','energy_vars(i,:)');

                end
    
                % if ~settings.zero_initial_data
                %     axis([parameters.time_initial parameters.time_final energy(1)*.8 energy(1)*1.2]);
                % end







                hold off;
            end

            if settings.plot_spectrum
                obj.spec_figure = figure();

                figure(obj.spec_figure);
                hold on;

                % title(sprintf('Spectrum at %.2f', parameters.time_current), "FontSize", obj.fontsize, "FontUnits",'points');
                xlabel("Wave Modes","FontSize", obj.fontsize);
                ylabel("Magnitude", "FontSize", obj.fontsize);
                set(gca, 'LineWidth', obj.linewidth);


                spec = [];
                modes = 1:parameters.Nx/2;

                spec_temp = generate_spectrum(psi_hat);
                spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;
                spec(1,:) = spec_temp;
                obj.spec_plots(1) = loglog(modes, spec(1,:), 'LineWidth', obj.linewidth);
                grid on;
                % obj.spec_plots(1) = plotSpectrum(psi_hat);
                count = 0;
                for i = 1: length(vars)
                    if vars(i).assimilate_type ~= "EnKF"
                        % M = 1
                        spec_temp = generate_spectrum(vars(i).v_hat);
                        spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;
                        spec(i+1 + count,:) = spec_temp;
                        obj.spec_plots(i+1 + count) = loglog(modes, spec(i+1 + count,:), 'LineWidth', obj.linewidth);

                    else
                        spec_temp = generate_spectrum(vars(i).v_hat);
                        spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;
                        spec(i+1 + count,:) = spec_temp;
                        obj.spec_plots(i+1 + count) = loglog(modes, spec(i+1 + count,:), 'LineWidth', obj.linewidth);
                        M = vars(i).ensemble_count;
                        if settings.plot_spectrum_ensembles
                        for j = 1:M
                            count = count+1;
                            spec_temp = generate_spectrum(vars(i).ensemble(j).analysis);
                            spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;
                            spec(i+1+count,:) = spec_temp;
                            obj.spec_plots(i+1 + count) = loglog(modes, spec(i+1 + count,:), 'LineWidth', obj.linewidth);
                            % obj.spec_plots(i+1) = plotSpectrum(vars(i).v_hat);
                            hold on;
                        end
                        end
                       

                    end
                end


                dealias_line = (2/3)*(parameters.Nx/2);
                line([dealias_line dealias_line], [1e-20 1e+2], 'Color', 'r', 'LineWidth', 2);  % plots a vertical dashed red line

                set(gca,'XScale','log');
                set(gca,'YScale','log');

                axis([0 parameters.Nx/2 1e-20 1e+2]);

                hold off;
            end

            obj.da_plotters = [];
            if obj.settings.plot_vars
                da_plotters = repelem(DA_plotter,1,1);
                for i = 1:length(vars)
                    da_plotters(i) = DA_plotter(parameters, settings, psi_hat, vars(i));
                    % vars(i) = vars(i).initialize_plots(p,psi_hat);
                end
                obj.da_plotters = da_plotters;
            end



            if settings.plot_error
                obj.error_figure_low = figure();
                figure(obj.error_figure_low);
                hold off;
                for i = 1: length(vars)
                    semilogy(0:parameters.dt:parameters.time_final+parameters.dt,vars(i).error_low,'-','LineWidth',1);
                    hold on;
                end
                grid on;
                title('Error of Simulated Solutions on Observed Modes');
                xlabel('Time');
                ylabel('Error');
                axis([parameters.time_initial parameters.time_final 1e-20 1e-0]);

                obj.error_figure_high = figure();
                figure(obj.error_figure_high);
                hold off;
                for i = 1: length(vars)
                    semilogy(0:parameters.dt:parameters.time_final+parameters.dt,vars(i).error_high,'-','LineWidth',1);
                    hold on;
                end
                grid on;
                title('Error of Simulated Solutions on Unobserved Modes');
                xlabel('Time');
                ylabel('Error');
                axis([parameters.time_initial parameters.time_final 1e-20 1e-0]);


                obj.error_figure_total = figure();
                figure(obj.error_figure_total);
                hold off;
                for i = 1: length(vars)
                    semilogy(0:parameters.dt:parameters.time_final+parameters.dt,vars(i).error,'-','LineWidth',1);
                    hold on;
                end
                grid on;
                title('Error of Simulated Solutions');
                xlabel('Time');
                ylabel('Error');
                axis([parameters.time_initial parameters.time_final 1e-20 1e-0]);


            % semilogy(0:p.dt:p.time_final+p.dt,obj.var.error,'-','LineWidth',1);
            % hold on;
            % semilogy(0:p.dt:p.time_final+p.dt,obj.var.error_low,'--','LineWidth',1);
            % semilogy(0:p.dt:p.time_final+p.dt,obj.var.error_high,'--','LineWidth',1);


            end

            drawnow;

            % hold off;
        end

        function obj = update(obj, parameters, psi_hat, vars)


            if obj.settings.plot_ref
                figure(obj.ref_figure);
                k_lap = obj.parameters.k_lap;
                set(obj.ref_plot,'cdata',ifftn( -k_lap.*psi_hat,'symmetric'));
                title(sprintf('Reference Vorticity at t = %1.2f',parameters.time_current));

            end


            if(obj.settings.plot_spectrum)
                figure(obj.spec_figure);
                hold on;

                % title(sprintf('Spectrum at %.2f', parameters.time_current),"FontSize", obj.fontsize, "FontUnits",'points');




                spec = [];
                % modes = 1:p.Nx/2;

                spec_temp = generate_spectrum(psi_hat);

                spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;
                spec(1,:) = spec_temp;
                set(obj.spec_plots(1), 'ydata', spec(1,:));
                count = 0;
                % obj.spec_plots(1) = loglog(modes, spec(1,:));
                % obj.spec_plots(1) = plotSpectrum(psi_hat);
                for i = 1: length(vars)

                    if vars(i).assimilate_type ~= "EnKF"
                    spec_temp = generate_spectrum(vars(i).v_hat);
                    spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;

                    spec(i+1+count,:) = spec_temp;
                    set(obj.spec_plots(i+1+count), 'YData', spec(i+1+count,:));
                    % obj.spec_plots(i+1) = loglog(modes, spec(i+1,:));
                    % obj.spec_plots(i+1) = plotSpectrum(vars(i).v_hat);
                    else
                        spec_temp = generate_spectrum(vars(i).v_hat);
                        spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;

                        spec(i+1+count,:) = spec_temp;
                        set(obj.spec_plots(i+1+count), 'YData', spec(i+1+count,:));
                        M = vars(i).ensemble_count;
                        if obj.settings.plot_spectrum_ensembles
                        for j = 1:M
                            count = count+1;
                            spec_temp = generate_spectrum(vars(i).ensemble(j).analysis);
                            spec_temp(spec_temp < 1e-20 | isnan(spec_temp)) = 1e-20;

                            spec(i+1+count,:) = spec_temp;
                            set(obj.spec_plots(i+1+count), 'YData', spec(i+1+count,:));

                        end
                        end


                    end
                end
                hold off;


            end

            if obj.settings.plot_energy
                % This routine updates the plots of the kinetic energy.
                % DO NOT TRY TO SIMPLIFY THIS!
                % refreshdata cannot access properties of classes, so you need temporary
                % variables. Also the 'caller' ensures that the data is pulled from the
                % workspace of THIS function.
                figure(obj.en_figure);
                % title(sprintf('Kinetic Energy%.2f', parameters.time_current),"FontSize", obj.fontsize, "FontUnits",'points');

                energy = parameters.energy;

                refreshdata(obj.en_plots(1),'caller')
                energy_vars = [];
                for i = 1:length(vars)
                    energy_vars(i,:) = vars(i).energy;

                    refreshdata(obj.en_plots(i+1),'caller')
                end

                % refreshdata(obj.en_plots(1))

                % hold on;
                %
                % title(sprintf('Kinetic Energy at %.2f', parameters.time_current), "FontSize", obj.fontsize, "FontUnits",'points');
                % xlabel("Time","FontSize", obj.fontsize);
                % ylabel("Magnitude", "FontSize", obj.fontsize);
                % set(gca, 'LineWidth', obj.linewidth);
                %
                % t = parameters.t;
                %
                % obj.en_plots = [];
                %
                % obj.en_plots(1) = plot(t, parameters.energy, "LineWidth",obj.linewidth);
                % % set(obj.en_plots(1), 'YDataSource','parameters.energy');
                % % refreshdata(obj.en_plots(1));
                % for i = 1:length(vars)
                %     obj.en_plots(i+1) = plot(t, vars(i).energy, "LineWidth", obj.linewidth);
                %     % set(obj.en_plots(i+1), 'YDataSource','vars(i).energy');
                %
                % end





            end

            if obj.settings.plot_vars
                for i = 1:length(vars)
                    % Note: each of these plotter objects have a var
                    % associated with them. var is passed by reference, so
                    % the updating can occur without needing to pass values
                    % back and forth.
                    obj.da_plotters(i).update_figures(parameters,obj.settings,psi_hat);
                end


            end
            if obj.settings.plot_error
                            % obj.error_figure_low = figure();
                figure(obj.error_figure_low);
                hold off;
                for i = 1: length(vars)
                    semilogy(0:parameters.dt:parameters.time_final+parameters.dt,vars(i).error_low,'-','LineWidth',1);
                    hold on;
                end
                title('Error of Simulated Solutions on Observed Modes');
                xlabel('Time');
                ylabel('Error');
                axis([parameters.time_initial parameters.time_final 1e-20 1e-0]);

                % obj.error_figure_high = figure();
                figure(obj.error_figure_high);
                hold off;
                for i = 1: length(vars)
                    semilogy(0:parameters.dt:parameters.time_final+parameters.dt,vars(i).error_high,'-','LineWidth',1);
                    hold on;
                end
                title('Error of Simulated Solutions on Unobserved Modes');
                xlabel('Time');
                ylabel('Error');
                axis([parameters.time_initial parameters.time_final 1e-20 1e-0]);



                % obj.error_figure_total = figure();
                figure(obj.error_figure_total);
                hold off;
                for i = 1: length(vars)
                    semilogy(0:parameters.dt:parameters.time_final+parameters.dt,vars(i).error,'-','LineWidth',1);
                    hold on;
                end
                title('Error of Simulated Solutions');
                xlabel('Time');
                ylabel('Error');
                axis([parameters.time_initial parameters.time_final 1e-20 1e-0]);

            end


            drawnow;

        end

        function save(obj, saveDir)


            if obj.settings.plot_energy
                saveName = 'energy';
                savePathFig = fullfile(saveDir, [saveName, '.fig']);
                savePathJPG = fullfile(saveDir, [saveName, '.jpg']);
                saveas(obj.en_figure, savePathFig);
                saveas(obj.en_figure, savePathJPG);
            end

            if obj.settings.plot_ref
                saveName = 'ref_soln';
                savePathFig = fullfile(saveDir, [saveName, '.fig']);
                savePathJPG = fullfile(saveDir, [saveName, '.jpg']);
                saveas(obj.ref_figure, savePathFig);
                saveas(obj.ref_figure, savePathJPG);
            end

            if obj.settings.plot_spectrum
                saveName = 'spectrum';
                savePathFig = fullfile(saveDir, [saveName, '.fig']);
                savePathJPG = fullfile(saveDir, [saveName, '.jpg']);
                saveas(obj.spec_figure, savePathFig);
                saveas(obj.spec_figure, savePathJPG);
            end

            if obj.settings.plot_error
                saveName = 'error_low';
                savePathFig = fullfile(saveDir, [saveName, '.fig']);
                savePathJPG = fullfile(saveDir, [saveName, '.jpg']);
                saveas(obj.error_figure_low, savePathFig);
                saveas(obj.error_figure_low, savePathJPG);

                saveName = 'error_high';
                savePathFig = fullfile(saveDir, [saveName, '.fig']);
                savePathJPG = fullfile(saveDir, [saveName, '.jpg']);
                saveas(obj.error_figure_high, savePathFig);
                saveas(obj.error_figure_high, savePathJPG);


                saveName = 'error_total';
                savePathFig = fullfile(saveDir, [saveName, '.fig']);
                savePathJPG = fullfile(saveDir, [saveName, '.jpg']);
                saveas(obj.error_figure_total, savePathFig);
                saveas(obj.error_figure_total, savePathJPG);
            end

            


            if obj.settings.plot_vars
                for i = 1:length(obj.da_plotters)
                    saveName = 'DA';
                    savePathFig = fullfile(saveDir,[saveName, num2str(i), '.fig']);
                    savePathJPG = fullfile(saveDir,[saveName, num2str(i), '.jpg']);
                    saveas(obj.da_plotters(i).figure_window, savePathFig);
                    saveas(obj.da_plotters(i).figure_window, savePathJPG);
                end
            end








        end
    end
end