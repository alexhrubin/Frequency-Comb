% class holding the parameters of the waveguide
% The width method describes the parameterization of one edge
%   of the waveguide; the other edge is assumed to be along the z-axis.
classdef Waveguide
    properties
        n
        d_zero %this is the half-width at z=0
        k
        a %controls the rate at which the waveguide expands
    end
    methods
        function obj = Waveguide(n, d_zero, k, a) %half-width at z
            obj.n = n;
            obj.d_zero = d_zero;
            obj.k = k;
            obj.a = a;
        end
        function d = width(obj, z)
            d = obj.d_zero * exp(obj.a*z);
        end
        function p = dedz(obj, z)
            p = obj.a * width(obj, z);
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

        function profiles = mode_profiles(obj, z)
            num_of_points = 5000;
            d = obj.width(z); %half-width of waveguide at z
            X = linspace(-5*d, 5*d, num_of_points);

            betas = getbeta(obj, z);
            profiles = [];

            for i = 1:length(betas)
                u = sqrt((obj.n*obj.k)^2-betas(i)^2);
                w = sqrt(betas(i)^2-obj.k^2);
                mode = zeros(1, num_of_points);

                if mod(i,2) == 1 %symmetrical modes
                    for j = 1:length(X)
                        x = X(j);
                        C = cos(u*d)*exp(w*d);
                        if x < -d
                            mode(j) = C*exp(w*x);
                        elseif x > d
                            mode(j) = C*exp(-w*x);
                        else
                            mode(j) = cos(u*x);
                        end
                    end

                elseif mod(i,2) == 0 %asymmetrical modes
                    for j = i:length(X)
                        x = X(j);
                        C = sin(u*d)*exp(w*d);
                        if x < -d
                            mode(j) = -C*exp(w*x);
                        elseif x > d
                            mode(j) = C*exp(-w*x);
                        else
                            mode(j) = sin(u*x);
                        end

                    end
                end
                profiles = cat(1, profiles, mode);
            end
            profiles = cat(1, profiles, X);
        end

        
        
        function plot_modes(obj, mode_orders, z)
            profiles = obj.mode_profiles(z);

            if mode_orders == "all"
                mode_orders = (1:size(profiles,1)-1)
            end

            fig = figure;
            hold on
            grid on
            for i = mode_orders(1):mode_orders(end)
                plot(profiles(end,:), profiles(i,:))
            end

            plot([-obj.width(z) -obj.width(z)], [-1 1], 'k')
            plot([obj.width(z) obj.width(z)], [-1 1], 'k')

            hold off
            figure(fig)
        end
        
    end
end