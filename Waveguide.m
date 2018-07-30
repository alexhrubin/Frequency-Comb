% class holding the parameters of the waveguide
% The width method describes the parameterization of one edge
%   of the waveguide; the waveguide is assumed to be symmetric across the z
%   axis.
classdef Waveguide
    properties
        n           %index of the waveguide material (environment has n=1)
        k           %free-space wavenumber
        d_zero      %this is the half-width at z=0
        a           %controls the rate at which the waveguide expands
        betas       %holds the propagation constants of the modes
        
    end
    
    methods
        % commonly used specs: (sqrt(12), 2*pi/1.55, 0.2, 0.0015)
        function obj = Waveguide(n, k, d_zero, a)
            obj.n = n;
            obj.k = k;
            obj.d_zero = d_zero;
            obj.a = a;
        end
        
        function d = width(obj, z)
            %d = obj.d_zero * exp(obj.a * z);  %exponentially expanding
            d = obj.d_zero + obj.a * z;  %linearly expanding
            %d = obj.d_zero + obj.d_zero*exp(obj.a*z) +  0.1*sin(2*pi*z/5);
        end
        
        function p = dedz(obj, z)
           % p = obj.a * obj.d_zero * exp(obj.a * z);
           p = obj.a;
           % p = obj.a*obj.d_zero*exp(obj.a*z) + 0.1*2*pi*cos(2*pi*z/5) / 5;
        end
        
        function beta = getbeta(obj, z)
            d = obj.width(z);
            V = d * obj.k * sqrt(obj.n^2-1);
            n_sym= floor(V/pi);
            x_sym=[];
            for i_sym = 0:n_sym
                tmp_U = fminbnd(@(U) sym_modes(U, V), i_sym*pi, i_sym*pi+pi/2);
                tmp_beta = sqrt((obj.n*obj.k)^2-(tmp_U/d)^2);
                x_sym = [x_sym, tmp_beta];
            end

            n_asym= round(V/pi);
            x_asym=[];
            if n_asym>0
                for i_asym = 1:n_asym
                    tmp_U = fminbnd(@(U) asym_modes(U, V), pi/2+(i_asym-1)*pi, pi+(i_asym-1)*pi);
                    tmp_beta = sqrt((obj.n*obj.k)^2-(tmp_U/d)^2);
                    x_asym = [x_asym, tmp_beta];
                end
            end
            beta = sort([x_sym, x_asym], 'descend');
        end
  
        function fcn = eigenmode_function(obj, z, order)
            local_betas = obj.getbeta(z);
            if (order > length(local_betas))
                fcn = @(x) 0*x;
                return
            end
            
            d = obj.width(z);
            u = sqrt((obj.n*obj.k)^2-local_betas(order)^2);
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
            
            d = obj.width(z);
            X = linspace(-5*d, 5*d, num_of_points);
            
            local_betas = obj.getbeta(z);
                        
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

            plot([-obj.width(z) -obj.width(z)], [-1 1], 'k')
            plot([obj.width(z) obj.width(z)], [-1 1], 'k')

            hold off
            figure(fig)
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