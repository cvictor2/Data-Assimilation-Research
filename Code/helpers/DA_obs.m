classdef DA_obs < handle
    properties
        params
        observer_type
        assimilate_type
        assimilate_subtype
        IC_type

        mu

        v
        v_hat
        v_2
        v_3
        u_obs
        v_obs
        obs_noise

        trunc
        trunc_array
        trunc_index
        trunc_index_comp


        switch_times_mu
        mu_values
        mu_index

        obs_locations
        obs_indices
        obs_length

        obs_interpolant

        error
        error_low
        error_high

        error_observed

        ensemble_count
        ensemble

        energy

        cpu_time


    end
    methods
        %         function var = DA_obs(p, type_obs, type_assimilate, type_IC, mu, trunc_num, trunc_array, mu_values, mu_switch_times)
        function var = DA_obs(p, varargin)

            var.params = p;

            ip = inputParser;
            addParameter(ip, 'assimilate_type', "AOT");
            addParameter(ip, 'assimilate_subtype', "NA");
            addParameter(ip, 'observer_type', "Projection");
            addParameter(ip, 'IC_type', "Projection");
            addParameter(ip, 'trunc', p.trunc);
            addParameter(ip, 'mu', 100);
            addParameter(ip, 'switch_times_mu', [inf]);
            addParameter(ip, 'mu_values', [100]);
            addParameter(ip, 'ensemble_count',50);


            parse(ip, varargin{:});

            fields = fieldnames(ip.Results);
            for i = 1:length(fields)
                var.(fields{i}) = ip.Results.(fields{i});
            end






            %             var.params = p;
            %             var.observer_type = type_obs;
            %             var.assimilate_type = type_assimilate;
            %             var.IC_type = type_IC;
            %             var.mu = mu;
            %             var.trunc = trunc_num;

            var.v = zeros([p.Nx,p.Ny]);
            var.v_hat = zeros([p.Nx,p.Ny]);
            var.v_2 = zeros([p.Nx,p.Ny]);
            var.v_3 = zeros([p.Nx,p.Ny]);


            var.u_obs = zeros([p.Nx,p.Ny]);
            var.v_obs = zeros([p.Nx,p.Ny]);


            %By default the trunc number is given as in parameters. If this
            %value is overwritten, then regenerate the correct truncation
            %array.
            %% To do: Fix this, the current iteration will include the Nyquist frequency, which should be excluded.
            if p.trunc ~= var.trunc
                var.trunc_array = ones(p.Nx,p.Ny);
                kx    =      ([0:p.Nx/2, -p.Nx/2+1:-1]*(2*pi/p.Lx));
                ky    =    [0:p.Ny/2, -p.Ny/2+1:-1]*(2*pi/p.Ly).';
                for j = 1:p.Ny
                    for i = 1:p.Nx
                        if(kx(i)^2 + ky(j)^2 >= var.trunc^2)
                            var.trunc_array(i,j) = 0;
                        end
                    end
                end
                var.trunc_array(1,1) = 0;
            else
                var.trunc_array = p.trunc_array;
            end

            var.trunc_index = find(var.trunc_array == 1);
            var.trunc_index_comp = find(var.trunc_array == 0);


            if var.assimilate_type == "EnKF"
                
                var.ensemble = [];    
                var.ensemble = repelem(Ensemble(var, p),1,1);
                M = var.ensemble_count;
                for i = 1:M 
                    var.ensemble(i) = Ensemble(var, p);

                end


            else
                var.ensemble_count = 0;
                var.ensemble = [];
            end


            var.error = nan(1,p.num_timesteps + 1);
            var.error_low = nan(1,p.num_timesteps + 1);
            var.error_high = nan(1,p.num_timesteps + 1);
            var.error_observed = nan(1, p.num_timesteps +1);

            var.energy = p.t.*nan;

            var.cpu_time = zeros(1,p.num_timesteps);
%             if nargin >= 7
                %
                %                 var.obs_length = h;
                %                 index_length_x = floor(p.Lx/h);
                %                 index_length_y = floor(p.Ly/h);
                %                 if(index_length_x > p.Nx)
                %                     index_length_x = p.Nx;
                %                 end
                %                 if(index_length_y > p.Ny)
                %                     index_length_y = p.Ny;
                %                 end
                %                 %This gives a specified length scale between observations
                %                 %                 indices_x = 1:index_length_x:p.Nx;
                %                 %                 indices_y = 1:index_length_y:p.Ny;
                %                 % This gives a specified number of observed points
                %                 indices_x = floor(linspace(1, p.Nx, index_length_x));
                %                 indices_y = floor(linspace(1, p.Ny, index_length_y));
                %
                %                 [A,B] = meshgrid(indices_x,indices_y);
                %                 c=cat(2,A',B');
                %                 var.obs_indices =reshape(c,[],2);
                %                 var.obs_locations = var.obs_indices;
                %                 var.obs_locations(:,1) = p.x(var.obs_indices(:,1));
                %                 var.obs_locations(:,2) = p.y(var.obs_indices(:,2));
                %
                %
                %                 %             u_data2 = [];
                %                 r = var.obs_indices(:,1);
                %                 c = var.obs_indices(:,2);
                %                 %             r = var.i_nodes_coordinates(:,1);
                %                 %             c = var.i_nodes_coordinates(:,2);
                %                 x_pts_per = p.x(r)';
                %                 y_pts_per = p.y(c)';
                %
                %                 obs_data_1 = zeros(size(x_pts_per));
                %                 %             obs_data_2 = u2(sub2ind(size(u2),r,c));
                %                 %         u_obs1 = zeros(p.Nx, p.Ny);
                %                 %         u_obs2 = zeros(p.Nx, p.Ny);
                %
                %                 %             for(i = -1:1)
                %                 %                 for(j=-1:1)
                %                 %                     x_pts_per = [x_pts_per; x_pts + p.Lx*i];
                %                 %                     y_pts_per = [y_pts_per; y_pts + p.Ly*j];
                %                 %
                %                 %                     u_data1 = [u_data1; obs_data_1];
                %                 %                     %                     u_data2 = [u_data2; obs_data_2];
                %                 %                 end
                %                 %             end
                %                 %             [X,Y] = meshgrid(p.x(var.i_nodes_x), p.y(var.i_nodes_y));
                %                 %             interp2(X,Y, u1(var.i_nodes_x, var.i_nodes_y)' - v1(var.i_nodes_x, var.i_nodes_y)', p.X, p.Y, 'linear', 0);
                %                 var.obs_interpolant = scatteredInterpolant(x_pts_per,y_pts_per, obs_data_1);
                %                 %             F2 = scatteredInterpolant(x_pts_per,y_pts_per, u_data2);
                %                 %         F1.Method = 'nearest';
                %                 %         F2.Method = 'nearest';
                %
                %







                %                 var.i_nodes_array = var.v;
                %                 for j = 1: var.int_nodes_x
                %                     for k = 1: var.int_nodes_y
                %                         var.i_nodes_array(var.i_nodes_x(j),var.i_nodes_y(k)) = 1;
                %                     end
                %                 end
                %                 var.i_nodes_array = sparse(var.i_nodes_array);
                %                 [a,b,~] = find(var.i_nodes_array);
                %                 var.i_nodes_coordinates = [a,b];
                %                 var.nodes_coordinates = [p.x(var.i_nodes_coordinates(:,1))', p.y(var.i_nodes_coordinates(:,2))'];
                %

%                 var.mu_switch_times = [];
%                 var.mu_values = [];
if var.assimilate_type == "Hybrid"
    mu_index = 1;
    %                 if nargin == 9
    %                     var.mu_values = mu_values;
    %                     var.mu_switch_times = mu_switch_times;
    var.mu = var.mu_values(1);
    var.mu_index = 1;
end
if var.assimilate_subtype == "Adaptive"
    var.switch_times_mu = [0];
    var.mu_values = [var.mu];
end

% end
%             else
%                 var.obs_length = 0;
%             end
        end
        function var = initialize(var, psi_hat)
            switch var.IC_type
                case "Zero"
                    var.v_hat = zeros(var.params.Nx,var.params.Ny);
                case "Projection"
                    switch var.observer_type
                        case "Projection"
                            var.v_hat(var.trunc_index) = psi_hat(var.trunc_index);

                        case "COT"
                            psi = ifftn(psi_hat,'symmetric');
                            %                             [u1,u2] = psi_converter(psi_hat, var.params);
                            %                             u_obs = var.physical_extractor(u1,u2);
                            obs = var.obs_extractor(psi);
                            obs_hat = fftn(obs);
                            obs_hat(var.trunc_index_comp) = 0;
                            obs_hat(1,1) = 0;

                            var.v_hat = obs_hat;

                    end
            end
            var.u_obs = var.v_hat;
            total_diff = var.v_hat-psi_hat;
            low_diff = norm(total_diff(var.trunc_index))/sqrt(var.params.Nx*var.params.Ny);
            high_diff = norm(total_diff(var.trunc_index_comp))/sqrt(var.params.Nx*var.params.Ny);
            low_diff = max(low_diff, var.params.eps);

            var.error_low(1) = low_diff;
            var.error_high(1) = high_diff;
            var.error(1) = low_diff + high_diff;

            var.error_observed(1) = low_diff;


        end
        function var = observe(var, psi_hat_old, psi_hat)
            %AOT need old data to predict with, otherwise it will be off by
            %a time-step.
            %Synchronization needs new data to enforce equality in the next
            %time-step.
            switch var.observer_type
                case "Projection"
                    switch var.assimilate_type
                        case "AOT"
                            if var.mu < 1
                            %Forward Euler
                                                        var.u_obs(var.trunc_index) = psi_hat_old(var.trunc_index);
                                                        var.v_obs(var.trunc_index) = var.v_hat(var.trunc_index);
                            else
                            %Backward Euler
                            var.u_obs(var.trunc_index) = psi_hat(var.trunc_index);
                            end


                            var.u_obs(var.trunc_index_comp) = 0;
                            var.v_obs(var.trunc_index_comp) = 0;

                        case "Synchronization"
                            var.u_obs(var.trunc_index) = psi_hat(var.trunc_index);
                            var.u_obs(var.trunc_index_comp) = 0;
                            var.v_obs(:) = 0;
                            %                             var.u_obs(var.trunc_index) = var.u_obs;
                            %                             var.u_obs(var.trunc_index_comp) = 0;
                        case "Hybrid"
                            %Forward Euler
                            %                             var.u_obs(var.trunc_index) = psi_hat_old(var.trunc_index);
                            %                             var.v_obs(var.trunc_index) = var.v_hat(var.trunc_index);

                            %Backward Euler
                            var.u_obs(var.trunc_index) = psi_hat(var.trunc_index);


                            var.u_obs(var.trunc_index_comp) = 0;
                            var.v_obs(var.trunc_index_comp) = 0;
                    end
                case "COT"
                    switch var.assimilate_type
                        case "AOT"
                            [u1,u2] = psi_converter(psi_hat_old, var.params);
                            var.u_obs = var.physical_extractor(u1,u2);
                            %                             var.u_obs(var.trunc_index_comp) = 0;


                            [v1,v2] = psi_converter(var.v_hat, var.params);
                            var.v_obs = var.physical_extractor(v1,v2);
                            %                             var.v_obs(var.trunc_index_comp) = 0;

                        case "Synchronization"
                            %This does it at the velocity level
                            %                             [u1,u2] = psi_converter(psi_hat, var.params);
                            %                             var.u_obs = var.physical_extractor(u1,u2);
                            %                             var.u_obs(var.trunc_index_comp) = 0;

                            %This does it at the stream function level
                            psi = ifftn(psi_hat,'symmetric');
                            %                             var.obs_interpolant.Value = [];
                            obs = var.obs_extractor(psi);
                            obs_hat = fftn(obs);
                            obs_hat(var.trunc_index_comp) = 0;
                            obs_hat(1,1) = 0;
                            var.u_obs = obs_hat;


                            var.v_obs(:) = 0;

                    end
            end

        end
        function var = assimilate(var, ev, v_pred, t)
            %Prediction: Note that the predicted evolution is passed into
            %this function
            %Correction (and assimilation)
            switch var.assimilate_type
                case "AOT"
                    if var.mu < 1
                                        var.v_hat = v_pred + ev.dt*var.mu*(var.u_obs - var.v_obs); %This is forward euler.
                    else
                    %This is backward Euler
                    var.v_hat(var.trunc_index) = (v_pred(var.trunc_index) + ev.dt*var.mu*var.u_obs(var.trunc_index))/(1 + ev.dt*var.mu);
                    end
                    var.v_hat(var.trunc_index_comp) = v_pred(var.trunc_index_comp);
                    var.v_obs(var.trunc_index) = var.v_hat(var.trunc_index);
                    
                case "Synchronization"
                    var.v_hat(var.trunc_index_comp) = v_pred(var.trunc_index_comp);
                    var.v_hat(var.trunc_index) = var.u_obs(var.trunc_index);
                case "Hybrid"
                    if (t >= var.switch_times_mu(var.mu_index)) && (var.switch_times_mu(var.mu_index) + var.params.dt < t)
                        var.mu_index = var.mu_index + 1;
                        var.mu = var.mu_values(var.mu_index);
                    end
                    if var.mu == inf
                        var.v_hat(var.trunc_index_comp) = v_pred(var.trunc_index_comp);
                        var.v_hat(var.trunc_index) = var.u_obs(var.trunc_index);
                    else
                        var.v_hat(var.trunc_index) = (v_pred(var.trunc_index) + ev.dt*var.mu*var.u_obs(var.trunc_index))/(1 + ev.dt*var.mu);
                        var.v_hat(var.trunc_index_comp) = v_pred(var.trunc_index_comp);
                        var.v_obs(var.trunc_index) = var.v_hat(var.trunc_index);
                    end
                    %                     var.v_hat = var.u_obs;
                case "EnKF"


            
            end
        end
        % % function var = update_figures(var, p, psi_hat)
        %     figure(var.varfigure);
        %     %             figure(fh1);
        %     subplot(2,2,1);
        %     set(var.varplot(1),'cdata',ifftn( -p.k_lap.*psi_hat,'symmetric'));
        %     title(sprintf('Reference Vorticity at t = %1.2f',p.time_current));
        %     colorbar;
        %     %             axmin = max(max(max(psi)),0.25);
        %     %             caxis([-axmin axmin]);
        %     %             drawnow;
        %     %             if(options(7))
        %     %                 set(varplot(5),'ydata',var.nodes_coordinates(:,1));
        %     %                 set(varplots(i,5),'xdata',var.nodes_coordinates(:,2));
        %     %             end
        % 
        %     subplot(2,2,2);
        %     %             intdiff = var.mu*DA_interp(psi_hat, var,p);
        %     %             intdiff/(norm(intdiff(:))*sqrt(Lx*Ly)/Nx/Ny);
        %     %             axmax = max(max(abs(ifftn( -p.k_lap.*psi_hat,'symmetric'))));
        %     set(var.varplot(2),'cdata',ifftn( -p.k_lap.*var.u_obs,'symmetric'));
        %     title(sprintf('Interpolated Vorticity at t = %1.2f',p.time_current));
        %     colorbar;
        %     %             caxis([-axmax,axmax]);
        % 
        % 
        %     %             axmin = max(max(max(psi)),0.25);
        %     %             caxis([-axmin axmin]);
        %     %             drawnow;
        % 
        % 
        %     subplot(2,2,3);
        %     %              set(varplots(i,3),'cdata',( abs(ifftn(-p.k_lap.*fftn(var.v), 'symmetric' ) - ifftn( -p.k_lap.*psi_hat,'symmetric')) ));
        %     % %              set(varplots(i,3),'cdata',( ifftn(-p.k_lap.*fftn(var.v), 'symmetric' )));
        %     %              title(sprintf('Absolute Vorticity Difference at t = %1.2f',t(ti)));
        %     %              colorbar;
        %     %                 caxis([-.5,.5]);
        %     %              hold on;
        %     %             varplots(i,3) = semilogy(0:dt:T,var.error,'-','LineWidth',1);
        %     %             varplots(i,3) = semilogy(0:dt:T,var.ens_umv,'-','LineWidth',1);
        %     semilogy(0:p.dt:p.time_final+p.dt,var.error,'-','LineWidth',1);
        %     hold on;
        %     semilogy(0:p.dt:p.time_final+p.dt,var.error_low,'--','LineWidth',1);
        %     semilogy(0:p.dt:p.time_final+p.dt,var.error_high,'--','LineWidth',1);
        %     hold off;
        % 
        %     %             hold off;
        %     %                         legend('Enstrophy Error', 'Psi L^2 Error');
        %     legend('Psi L^2 Error (total)', 'Psi L^2 Error (low modes)', 'Psi L^2 Error (high modes)');
        % 
        %     title('Error of Simulated Solution');
        %     xlabel('Time');
        %     ylabel('Error');
        %     %             axis([0 T 1e-20 1e-0]);
        % 
        % 
        % 
        %     subplot(2,2,4);
        %     set(var.varplot(4),'cdata',( ifftn(-p.k_lap.*var.v_hat, 'symmetric' )));
        %     title(sprintf('Simulated Vorticity at t = %1.2f',p.time_current));
        %     colorbar;
        %     %                             caxis([-30,30]);
        %     %              hold on;
        % 
        % 
        %     drawnow;
        %     %
        %     %
        % end
        function obs = physical_extractor(var, u1, u2)
            p = var.params;

            x_pts_per = [];
            y_pts_per = [];
            u_data1 = [];
            u_data2 = [];
            r = var.obs_indices(:,1);
            c = var.obs_indices(:,2);
            %             r = var.i_nodes_coordinates(:,1);
            %             c = var.i_nodes_coordinates(:,2);
            x_pts = p.x(r)';
            y_pts = p.y(c)';

            obs_data_1 = u1(sub2ind(size(u1),r,c));
            obs_data_2 = u2(sub2ind(size(u2),r,c));
            %         u_obs1 = zeros(p.Nx, p.Ny);
            %         u_obs2 = zeros(p.Nx, p.Ny);

            for(i = -1:1)
                for(j=-1:1)
                    x_pts_per = [x_pts_per; x_pts + p.Lx*i];
                    y_pts_per = [y_pts_per; y_pts + p.Ly*j];

                    u_data1 = [u_data1; obs_data_1];
                    u_data2 = [u_data2; obs_data_2];
                end
            end




            %         [X,Y] = meshgrid(p.x(var.i_nodes_x), p.y(var.i_nodes_y));
            %         int_diff = interp2(X,Y, psi(var.i_nodes_x, var.i_nodes_y)' - var.v(var.i_nodes_x, var.i_nodes_y)', p.X, p.Y, 'linear', 0);
            %         psi = fft2(psi_hat);
            %         int_diff = zeros(size(psi));
            %         int_diff(:,var.i_nodes_y) = psi(:,var.i_nodes_y) - var.v(:, var.i_nodes_y);
            %         int_diff(var.i_nodes_x,:) = psi(var.i_nodes_x,:) - var.v( var.i_nodes_x,:);




            F1 = scatteredInterpolant(x_pts_per,y_pts_per, u_data1);
            F2 = scatteredInterpolant(x_pts_per,y_pts_per, u_data2);
            %         F1.Method = 'nearest';
            %         F2.Method = 'nearest';
            u_obs1 = F1(p.X,p.Y)';
            u_obs2 = F2(p.X,p.Y)';
            obs = velocity_converter(u_obs1,u_obs2, p);


        end
        function obs = obs_extractor(var, psi)
            p = var.params;

            %             x_pts_per = [];
            %             y_pts_per = [];
            %             u_data1 = [];
            %             u_data2 = [];
            r = var.obs_indices(:,1);
            c = var.obs_indices(:,2);
            %             r = var.i_nodes_coordinates(:,1);
            %             c = var.i_nodes_coordinates(:,2);
            x_pts_per = p.x(r)';
            y_pts_per = p.y(c)';

            obs_data_1 = psi(sub2ind(size(psi),r,c));
            %             obs_data_2 = u2(sub2ind(size(u2),r,c));
            %         u_obs1 = zeros(p.Nx, p.Ny);
            %         u_obs2 = zeros(p.Nx, p.Ny);

            %             for(i = -1:1)
            %                 for(j=-1:1)
            %                     x_pts_per = [x_pts_per; x_pts + p.Lx*i];
            %                     y_pts_per = [y_pts_per; y_pts + p.Ly*j];
            %
            %                     u_data1 = [u_data1; obs_data_1];
            %                     %                     u_data2 = [u_data2; obs_data_2];
            %                 end
            %             end
            %             [X,Y] = meshgrid(p.x(var.i_nodes_x), p.y(var.i_nodes_y));
            %             interp2(X,Y, u1(var.i_nodes_x, var.i_nodes_y)' - v1(var.i_nodes_x, var.i_nodes_y)', p.X, p.Y, 'linear', 0);
            %             F1 = scatteredInterpolant(x_pts_per,y_pts_per, obs_data_1);
            var.obs_interpolant.Values = obs_data_1;


            %             F2 = scatteredInterpolant(x_pts_per,y_pts_per, u_data2);
            %         F1.Method = 'nearest';
            %         F2.Method = 'nearest';
            obs = var.obs_interpolant(p.X,p.Y)';
            %             u_obs2 = F2(p.X,p.Y)';
            %             obs = velocity_converter(u_obs1,u_obs2, p);


        end
        function var = restart(var, p)
            var.params = p;
            % error_temp = var.error;
            error_new = ones(1,p.num_timesteps + 1)*nan;
            error_new(1:length(var.error)) = var.error;
            var.error = error_new;
            
            error_new_low = ones(1,p.num_timesteps + 1)*nan;
            error_new_low(1:length(var.error_low)) = var.error_low;
            var.error_low = error_new_low;

            error_new_high = ones(1,p.num_timesteps + 1)*nan;
            error_new_high(1:length(var.error_high)) = var.error_high;
            var.error_high = error_new_high;

        end
    end
end