classdef Parameters
    properties
        note
        evolution_scheme = "Implicit Euler"

        trunc = 50;
        trunc_array
        Nx
        Ny
        Nx_sq
        Ny_sq
        resolution = 2^10;
        time_initial
        time_current
        time_final = 10;
        t
        dt
        ti = 1;
        num_timesteps
        G = [];
        Lx
        Ly

        psi_hat_initial;

        kx
        ikx
        ky
        iky
        IKX
        IKY

        k_lap
        k_lap_inv
        k1k2
        k1sq_m_k2sq

        dealias_mask
        dealias_modes
        f_hat
        f1
        f2
        % helm_inv
        eps = eps

        coordinates
        L

        parseval

        dx
        dy
        x
        y
        X
        Y
        energy
        
        error %debugging term

        nu    %Viscosity (should be 0 for Euler)
        gamma %Damping coefficient

    end

    properties (Dependent)
        alpha
    end

    methods
        function params = Parameters(settings,varargin)
            p = inputParser;

            % Adding parameters
            addParameter(p, 'resolution', params.resolution);
            addParameter(p, 'G', params.G);
            addParameter(p, 'time_initial', params.time_initial);
            addParameter(p, 'time_current', params.time_current);
            addParameter(p, 'dt', params.time_current);
            addParameter(p, 'viscosity', -inf);  % corrected the typo
            addParameter(p, 'gamma', 0);
            addParameter(p, 'trunc', 50);

            params.time_final = settings.final_time;

            parse(p, varargin{:});

            fields = fieldnames(p.Results);

            for i = 1:length(fields)
                field_name = fields{i};

                if strcmp(field_name, 'viscosity')
                    if p.Results.viscosity ~= -inf  % corrected the typo
                        params.nu = p.Results.viscosity;
                    end
                else
                    params.(field_name) = p.Results.(field_name);
                end
            end

            % p = inputParser;
            % addParameter(p, 'resolution', params.resolution);
            % addParameter(p, 'G', params.G);
            % addParameter(p, 'time_initial', params.time_initial);
            % addParameter(p, 'time_current', params.time_current);
            % addParameter(p, 'dt', params.time_current);
            % addParameter(p, 'viscosity', -inf);
            %
            % % addParameter(p, 'time_final', params.time_final);
            %
            % params.time_final = settings.final_time;
            %
            % parse(p, varargin{:});
            %
            % fields = fieldnames(p.Results);
            % for i = 1:length(fields)
            %     if (fields{i}) == 'viscosity'
            %         if p.Results.viscocity ~= -inf
            %             params.nu = p.Results.viscosity;
            %         end
            %     else
            %     % if (fields{1}) ~= 'visocity'
            %         params.(fields{i}) = p.Results.(fields{i});
            %     % end
            %     end
            % end

            % Initialize other properties based on resolution
            params = initializeProperties(params,settings);
            if settings.zero_initial_data
                params.psi_hat_initial = params.psi_hat_initial.*0;
            end
        end

        % function nu = get.nu(obj)
        %     switch obj.resolution
        %         case 2^10
        %             nu = 0.0001;
        %         case 2^9
        %             nu = 0.0005;
        %     end
        % end

        function alpha = get.alpha(obj)
            switch obj.resolution
                case 2^10
                    alpha = 32/sqrt(obj.Nx^2+obj.Ny^2);
                case 2^9
                    alpha = 0.005;
            end
        end

        function obj = initializeProperties(obj,settings)

            switch obj.resolution
                case 2^10
                    obj.Nx = 2^10;
                    obj.Ny = 2^10;
                    % obj.dt = 0.001;
                    obj.psi_hat_initial = obj.loadPsiHatInitial();
                    if isempty(obj.G)
                        obj.G = 1e6;
                    end
                    if isempty(obj.dt)
                        obj.dt = 0.001;
                    end

                case 2^9
                    obj.Nx = 2^9;
                    obj.Ny = 2^9;
                    % obj.dt = 0.01;
                    obj.psi_hat_initial = obj.loadPsiHatInitial();
                    if isempty(obj.G)
                        obj.G = 1e5;
                    end
                      if isempty(obj.dt)
                        obj.dt = 0.01;
                      end
                case 2^6
                    obj.Nx = 2^6;
                    obj.Ny = 2^6; 
                    obj.psi_hat_initial = obj.loadPsiHatInitial();
                    if     isempty(obj.G)
                        obj.G = 1e5;
                    end
                    if isempty(obj.dt)
                        obj.dt = 0.005;
                    end

                otherwise
                    obj.Nx = obj.resolution;
                    obj.Ny = obj.resolution;
                    if isempty(obj.dt)
                        obj.dt = 0.01;
                    end
                    if isempty(obj.psi_hat_initial)
                        obj.psi_hat_initial = zeros(obj.Nx,obj.Ny);
                    end
                    if isempty(obj.G)
                        obj.G = 1;
                    end

            end
            if isempty(obj.nu)
                switch obj.resolution
                    case 2^10
                        obj.nu = 0.0001;
                    case 2^9
                        obj.nu = 0.0005;
                end
            end

            if isempty(obj.time_initial)
                obj.time_initial = 0;
            end
            if isempty(obj.time_current)
                obj.time_current = obj.time_initial;
            end
            if isempty(obj.time_final)
                obj.time_final = 10;
            end
            if isempty(obj.t)
                obj.t = obj.time_initial:obj.dt:obj.time_final;
                obj.num_timesteps = length(obj.t);
            end

            obj.energy = obj.t.*nan;
            
            obj.error = ones(obj.num_timesteps +1,1).*NaN;
   

            % Make a mask for dealiasing
            % dealias_modes_x = ceil(Nx/3):(Nx - floor(Nx/3) + 1);
            % dealias_modes_y = ceil(Ny/3):(Ny - floor(Ny/3) + 1);
            cutoff = 1/3;
            dealias_mask = zeros(obj.Nx,obj.Ny);
            dealias_mask(ceil(cutoff*obj.Nx):(obj.Nx - floor(cutoff*obj.Nx) + 1),:) = 1;
            dealias_mask(:,ceil(cutoff*obj.Ny):(obj.Ny - floor(cutoff*obj.Ny) + 1)) = 1;
            obj.dealias_modes = [1;find(dealias_mask == 1)]; % Also enforce mean-zero




            x_min = -pi;
            x_max = pi;
            y_min = -pi;
            y_max = pi;

            Lx = x_max - x_min;
            Ly = y_max - y_min;
            L = sqrt(Lx*Ly);

            obj.Lx = Lx;
            obj.Ly = Ly;
            obj.L = L;

            obj.parseval = sqrt(Lx*Ly)/obj.Nx/obj.Ny; % Multiply by norm(u_hat(:)) to get L^2 norm

            % Physical space
            dx = (x_max - x_min)/obj.Nx; % Assume uniform grid.
            x  = x_min + (0:(obj.Nx-1))*dx;
            dy = (y_max - y_min)/obj.Ny; % Assume uniform grid.
            y  = y_min + (0:(obj.Ny-1))*dy;

            obj.dx = dx;
            obj.dy = dy;
            obj.x = x;
            obj.y = y;
            [obj.X,obj.Y] = meshgrid(x,y);


            x_pts = reshape(obj.X, obj.Nx^2,1);
            y_pts = reshape(obj.Y, obj.Ny^2,1);
            obj.coordinates = [x_pts, y_pts];
            % x_pts_per = [];
            % y_pts_per = [];
            % for(i = -1:1)
            %     for(j=-1:1)
            %         x_pts_per = [x_pts_per; x_pts + Lx*i];
            %         y_pts_per = [y_pts_per; y_pts + p.Ly*j];
            %     end
            % end

            % p.coordinates_periodic = [x_pts_per,y_pts_per];




            % Wave numbers k (we also multiply k by i for convenience).
            kx    =      ([0:obj.Nx/2, -obj.Nx/2+1:-1]*(2*pi/obj.Lx));
            obj.ikx = (1i*[0:((obj.Nx/2)-1) 0 ((-obj.Nx/2)+1):-1]*(2*pi/obj.Lx)).';
            obj.kx = ([0:((obj.Nx/2)-1) 0 ((-obj.Nx/2)+1):-1]*(2*pi/obj.Lx)).';
            kx_sq =    ([0:obj.Nx/2, -obj.Nx/2+1:-1].^2*(2*pi/obj.Lx)^2).';

            % Wave numbers k (we also multiply k by i for convenience).
            ky    =    [0:obj.Ny/2, -obj.Ny/2+1:-1]*(2*pi/obj.Ly).';
            obj.iky = 1i*[0:((obj.Ny/2)-1) 0 ((-obj.Ny/2)+1):-1]*(2*pi/obj.Ly);
            obj.ky = [0:((obj.Ny/2)-1) 0 ((-obj.Ny/2)+1):-1]*(2*pi/obj.Ly);
            ky_sq =    [0:obj.Ny/2, -obj.Ny/2+1:-1].^2*(2*pi/obj.Ly)^2;

            % Fourier Truncation of high wave modes
            % obj.trunc_num = 50; %Arbitrary choice of truncation numer
            obj.trunc_array = ones(obj.Nx,obj.Ny);
            for j = 1:obj.Ny
                for i = 1:obj.Nx
                    if(kx(i)^2 + ky(j)^2 >= obj.trunc^2)
                        obj.trunc_array(i,j) = 0;
                        % [i,j]
                    end
                end
            end
            obj.trunc_array(1,1) = 0;
            [obj.IKX, obj.IKY] = meshgrid(obj.ikx, obj.iky);
            obj.IKX = obj.IKX.';
            obj.IKY = obj.IKY.';
            % [KX, KY] = meshgrid(kx, ky);
            %
            % trunc_mask = (KX.^2 + KY.^2 >= obj.trunc^2);
            % obj.trunc_array(trunc_mask) = 0;

            % spy(obj.trunc_array)

            % Wave numbers of Laplacian operator.
            obj.k_lap     = zeros(obj.Nx,obj.Ny);
            obj.k_lap_inv = zeros(obj.Nx,obj.Ny);
            obj.k1sq_m_k2sq = zeros(obj.Nx,obj.Ny);
            % obj.helm_inv = zeros(obj.Nx,obj.Ny);

            % [KX_SQ, KY_SQ] = meshgrid(kx_sq, ky_sq);
            %
            % norm_k_sq = KX_SQ + KY_SQ;
            % obj.k_lap = -norm_k_sq;
            % obj.k_lap_inv = -1 ./ norm_k_sq;
            %
            % obj.k1sq_m_k2sq = KX_SQ - KY_SQ;
            % k1k2_matrix = KX_SQ.* KY_SQ;
            % obj.k1k2 = k1k2_matrix;
            %
            % helm_matrix = 1./(1+obj.alpha^2*norm_k_sq);
            % obj.helm_inv = helm_matrix;

            for j = 1:obj.Ny
                for i = 1:obj.Nx
                    norm_k_sq = (kx_sq(i) + ky_sq(j));
                    obj.k_lap(i,j) = -norm_k_sq;
                    obj.k_lap_inv(i,j) = -1/norm_k_sq;

                    obj.k1sq_m_k2sq(i,j) = kx_sq(i)-ky_sq(j);
                    obj.k1k2(i,j) = kx(i)*ky(j);

                    % obj.helm_inv(i,j) = 1./(1+obj.alpha^2*norm_k_sq);
                end
            end
            obj.k_lap_inv(abs(obj.k_lap_inv) == Inf) = 0;
            obj.k_lap_inv(isnan(obj.k_lap_inv)) = 0;
            obj.k_lap_inv(1,1) = 0;

            % if (obj.alpha<sqrt(obj.eps))
                % obj.helm_inv = 1;
            % end

            obj.Nx_sq = obj.Nx^2;
            obj.Ny_sq = obj.Ny^2;



            if ~settings.zero_forcing
                if obj.resolution > 64
            f_hat = OlsonTitiForcing(obj.Nx,obj.Ny,obj.nu,obj.G);
            % obj.f_hat = OlsonTitiForcing(obj.Nx,obj.Ny,obj.gamma,obj.G);
            % obj.f_hat = OlsonTitiForcing(obj.Nx,obj.Ny,1,obj.G);
            % obj.f_hat = RandomForcing(obj.Nx,obj.Ny,obj.G, 10, 64);
            % [f_hat, f1, f2] = OlsonTitiForcing(obj.Nx,obj.Ny,1,obj.G);
            % f_hat = zeros(obj.Nx,obj.Ny);
            % f1 = f_hat;
            % f2 = f_hat;
            % f_hat = zeros(obj.Nx,obj.Ny);

            % f1 = 1/2.*(sin(2*obj.X));
            % f2 = 1/2.*(sin(2*obj.Y));
            [f1,f2] = psi_converter(f_hat, obj);
            % f_hat = velocity_converter(f1,f2,obj);
            % 1/2sin(2x)
            % 1/2sin(2y)
                else
                    f_hat = cos(5.*(obj.X + obj.Y));
                    f_norm = norm(f_hat, 'fro');
                    f_hat = (f_hat./f_norm).*(obj.G*(obj.nu^2));
                    f_hat = fftn(f_hat);

                    [f1, f2] = psi_converter(f_hat, obj);

                end
            else
                f_hat = zeros(obj.Nx,obj.Ny);
                f1 = zeros(obj.Nx,obj.Ny);
                f2 = zeros(obj.Nx,obj.Ny);
            end
            obj.f_hat = f_hat;
            obj.f1 = f1;
            obj.f2 = f2;

        end

        function psi_hat = loadPsiHatInitial(obj)
            if ispc
                prefix = '.\restart files\';
            else
                prefix = './restart files/';
            end
            switch obj.resolution
                case 2^10
                    file = 'restart1024.mat';
                case 2^9
                    file = 'restart512_12100.mat';
                case 2^6
                    file = 'restart64.mat';
            end
            data = load([prefix file], 'psi_hat_initial');
            psi_hat = data.psi_hat_initial;
        end
        function params = restart(params, settings)
            params.time_final = settings.final_time;
            params.t = params.time_initial:params.dt:params.time_final;
            Lt = length(params.t);
            en_temp = params.energy;
            params.energy = params.t .* nan;
            params.energy(1:length(en_temp)) = en_temp;
            params.num_timesteps = length(params.t);
        end
    end
end
