% class holding the parameters of the waveguide
% The width method describes the parameterization of one edge
%   of the waveguide; the waveguide is assumed to be symmetric across the z
%   axis.
classdef Waveguide
    properties
        eps           %epsilon of the waveguide material (environment has n=1)
        k           %free-space wavenumber
        upper_edge_func  %function handle describing upper edge of waveguide
        lower_edge_func  %function handle describing lower edge of waveguide
        upper_edge_dz    %function handle describing upper edge slope
        lower_edge_dz    %function handle describing lower edge slope
    end
    
    methods
        % commonly used specs: (sqrt(12), 2*pi/1.55, 0.25, 0.04)
    function obj = Waveguide(eps, k, lower_edge_func, upper_edge_func)
            obj.eps = eps;
            obj.k = k;
            
            obj.upper_edge_func = upper_edge_func;
            obj.lower_edge_func = lower_edge_func;
            
            % for convenience, we do symbolic differentiation from the
            % waveguide profile functions, then convert symbolic exprs
            % back to function handles because they evaluate faster
            syms f(x)
            f(x) = obj.upper_edge_func(x);
            fprime = diff(f, x);
            obj.upper_edge_dz = matlabFunction(fprime);
            
            syms g(x)
            g(x) = obj.lower_edge_func(x);
            gprime = diff(g, x);
            obj.lower_edge_dz = matlabFunction(gprime);
    end
           
    function d = half_width(obj, z)
           d = (obj.upper_edge_func(z) - obj.lower_edge_func(z)) / 2;
    end
        
        
        function betas = getbeta_position(obj, z, max_modes)
            d = obj.half_width(z);
            if (nargin == 3)
                betas = obj.getbeta_width(d, max_modes);
            else
                betas = obj.getbeta_width(d);
            end
        end
        
        function betas = getbeta_width(obj, d, max_modes)
            V = d * obj.k * sqrt(obj.eps-1);

            if (nargin == 3)
                n_modes = max_modes;
            else
                n_modes = ceil(2*V/pi);    
            end

            betas = [];
            
            for i = 0 : n_modes-1
                if (i >= ceil(2*V/pi)) %happens if we demand a mode that doesn't exist
                    tmp_beta = 0;
                elseif (mod(i, 2) == 0) %symm modes
                    i_sym = i/2;
                    tmp_U = fminbnd(@(U) sym_modes(U, V), i_sym*pi, i_sym*pi+pi/2);
                    tmp_beta = sqrt(obj.eps*obj.k^2-(tmp_U/d)^2);
                elseif (mod(i, 2) == 1) %asymm modes
                    i_asym = ceil(i/2)-1;
                    tmp_U = fminbnd(@(U) asym_modes(U, V), i_asym*pi + pi/2, i_asym*pi + pi);
                    tmp_beta = sqrt(obj.eps*obj.k^2-(tmp_U/d)^2);
                end
                betas = [betas, tmp_beta];
            end
        end
  
        function fcn = eigenmode_function(obj, z, order)
            local_betas = obj.getbeta_position(z);
            if (order > length(local_betas))
                fcn = @(x) 0*x;
                return
            end
            
            d = obj.half_width(z);
            u = sqrt(obj.eps*obj.k^2-local_betas(order)^2);
            w = sqrt(local_betas(order)^2-obj.k^2);
            
            if (mod(order, 2) == 1) %symmetrical modes
                fcn = @(x) symmetrical_eigenmode(x, d, obj.k, local_betas(order), u, w);
            elseif (mod(order, 2) == 0) %asymmetrical modes
                fcn = @(x) asymmetrical_eigenmode(x, d, obj.k, local_betas(order), u, w);
            end
        end
        
        function fcn = eigenmode_h(obj, z, order)
            local_betas = obj.getbeta(z);
            electric_field = obj.eigenmode_function(z, order);
            fcn = @(x) (1/obj.k) * local_betas(order) * electric_field(x);
        end
        
        function plot_eigenmodes(obj, mode_orders, z, num_of_points)
            if (nargin < 4)
                num_of_points = 1000;
            end
            
            d = obj.half_width(z);
            X = linspace(-3*d, 3*d, num_of_points);
            
            local_betas = obj.getbeta_position(z);
                        
            if strcmp(mode_orders, 'all')
                mode_orders = 1 : length(local_betas);
            end

            fig = figure;
            hold on
            grid on
            for i = mode_orders
                mode_function = obj.eigenmode_function(z, i);
                plot(X, arrayfun(mode_function, X));
            end

            plot([-obj.half_width(z) -obj.half_width(z)], [-0.75 0.75], 'k')
            plot([obj.half_width(z) obj.half_width(z)], [-0.75 0.75], 'k')

            hold off
            figure(fig)
        end
        
        function visualize_waveguide(obj, z_range)
            z_array = z_range(1) : 0.05 : z_range(2);
            
            hold on
            grid on
            plot(z_array, arrayfun(@(z) obj.upper_edge_func(z), z_array))
            plot(z_array, arrayfun(@(z) obj.lower_edge_func(z), z_array))
            hold off        
        end
        
        
    end
end


%Begin eigenfunctions for solving propagation constants
function err=sym_modes(U, V)
%dispersion equations
W=sqrt(V^2-U^2);

%Mode equations
err=abs(U*tan(U)-W);
end

function err=asym_modes(U, V)
%dispersion equations
W=sqrt(V^2-U^2);

%Mode equations
err=abs(-U*cot(U)-W);
end

function y = symmetrical_eigenmode(x, d, k, beta, u, w)
normalization = sqrt(k/(2*beta)) / sqrt((cos(d*u)^2)/w + d + sin(2*d*u)/(2*u));
    C = cos(u*d)*exp(w*d);
    if x < -d
        y = normalization * C*exp(w*x);
    elseif x > d
        y = normalization * C*exp(-w*x);
    else
        y = normalization * cos(u*x);
    end
end

function y = asymmetrical_eigenmode(x, d, k, beta, u, w)
normalization = sqrt(k/(2*beta)) / sqrt( (sin(d*u)^2)/w + d - sin(2*d*u)/(2*u));
 C = sin(u*d)*exp(w*d);
 if x < -d
     y = -normalization * C*exp(w*x);
 elseif x > d
     y = normalization * C*exp(-w*x);
 else
     y = normalization * sin(u*x);
 end
end