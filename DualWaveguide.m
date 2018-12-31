classdef DualWaveguide
    properties
        wg1_lower_edge %lower waveguide lower edge
        wg1_upper_edge %lower waveguide upper edge
        wg2_lower_edge %upper waveguide lower edge
        wg2_upper_edge %upper waveguide upper edge
        
        wg1  %Waveguide object describing the lower waveguide
        wg2  %Waveguide object describing the upper waveguide
        
        sep     %separation of waveguide centers at z=0
        eps
        k
    end
    
    methods
        function obj = DualWaveguide(eps, k, sep, wg1_lower_edge, wg1_upper_edge, wg2_lower_edge, wg2_upper_edge)
            obj.wg1_lower_edge = @(z) wg1_lower_edge(z);
            obj.wg1_upper_edge = @(z) wg1_upper_edge(z);
            obj.wg2_lower_edge = @(z) wg2_lower_edge(z) + sep;
            obj.wg2_upper_edge = @(z) wg2_upper_edge(z) + sep;
            
            obj.wg1 = Waveguide(eps, k, wg1_lower_edge, wg1_upper_edge);
            obj.wg2 = Waveguide(eps, k, wg2_lower_edge, wg2_upper_edge);
            
            obj.sep = sep;
            obj.eps = eps;
            obj.k = k;
        end
        
        function d = edge_sep(obj, z)
            d = obj.sep + obj.wg2.lower_edge_func(z) - obj.wg1.upper_edge_func(z);
        end
        
        %Provides an initial guess as to the dual-wg betas
        function betas = naked_wg_betas(obj, z)
            betas = sort([obj.wg1.getbeta_position(z), obj.wg2.getbeta_position(z)]);
            betas = [betas*0.98, betas*1.01];
%             for i = 1 : length(betas)
%                 if (mod(i,2) == 1)
%                     betas(i) = 0.98 * betas(i);
%                 else
%                     betas(i) = 1.01 * betas(i);
%                 end
%             end  
        end
        
        function betas = getbeta_position(obj, z, starting_guess)
            a1 = 2 * obj.wg1.half_width(z);
            a2 = 2 * obj.wg2.half_width(z);
            d = obj.edge_sep(z);
           
            eigenfunc = @(b) obj.boundary_conds(z, b, 'det');
          %  eigenfunc = @(b) eigenproblem(b, obj.k, obj.eps, a1, a2, d);
            
          %  plot(obj.k:0.001:14.4, arrayfun(eigenfunc, obj.k:0.001:14.4));
            plot(0:0.001:20, arrayfun(eigenfunc, 0:0.001:20));
            grid on
            
            options = optimset('Display', 'off');
            
            betas = [];
            for i = starting_guess
                next_beta = fsolve(eigenfunc, i, options);
                % next_beta = fsolve(@(b) eigenproblem(b, obj.k, obj.eps, a1, a2, d), i, options)
                if (imag(next_beta) ~= 0)
                    continue
                end
                if (next_beta - obj.k < 1e-10)
                    next_beta = obj.k;
                end
                betas = [betas, next_beta];
                
            end
        end
        
        function f = eigenmode_function(obj, z, beta)
            a1 = 2 * obj.wg1.half_width(z);
            a2 = 2 * obj.wg2.half_width(z);
            d = obj.edge_sep(z);
            consts = obj.boundary_conds(z, beta, 'null');
            
            f = @(y) modefcn(y, consts, a1, a2, d, 2*pi/1.55, 12, beta);
        end
        
        function v = boundary_conds(obj, z, beta, mode)
            gam = sqrt(obj.eps*obj.k^2 - beta^2);
            eta = sqrt(beta^2 - obj.k^2);
            
            a1 = 2 * obj.wg1.half_width(z);
            a2 = 2 * obj.wg2.half_width(z);
            d = obj.edge_sep(z);

           
             M = [1 0 -1 0 0 0 0 0;
                 eta -gam 0 0 0 0 0 0;
                  0, sin(gam*a1), cos(gam*a1), -1, -1, 0, 0, 0;
                  0, gam*cos(gam*a1), -gam*sin(gam*a1), -eta, eta, 0, 0, 0;
                  0, 0, 0, exp(eta*d), exp(-eta*d), 0, -1, 0;
                  0, 0, 0, eta*exp(eta*d), -eta*exp(-eta*d), -gam, 0, 0;
                  0, 0, 0, 0, 0, sin(gam*a2), cos(gam*a2), -1;
                  0, 0, 0, 0, 0, gam*cos(gam*a2), -gam*sin(gam*a2), eta];
              
            if strcmp(mode, 'det')
                v = det(M);
            elseif strcmp(mode, 'null')
                v = null(M);
            else
                disp("ERROR")
            end
        end
        

        
        function plot_eigenmodes(obj, mode_orders, z)
            nb = obj.naked_wg_betas(z);
            local_betas = obj.getbeta_position(z, nb);
            
            if strcmp(mode_orders, 'all')
                mode_orders = 1 : length(local_betas);
            end
            
            Y = -0.5 : 0.0001 : 2;
            figure
            hold on
            grid on
            
            for i = mode_orders
                modefcn = obj.eigenmode_function(z, local_betas(i));
                plot(Y, arrayfun(modefcn, Y))
            end
            
            a = obj.wg1_lower_edge(z);
            b = obj.wg1_upper_edge(z);
            c = obj.wg2_lower_edge(z);
            d = obj.wg2_upper_edge(z);
            
            plot([a a], [-0.5 0.5], 'k')
            plot([b b], [-0.5 0.5], 'k')
            plot([c c], [-0.5 0.5], 'k')
            plot([d d], [-0.5 0.5], 'k')
            
            hold off
            
        end
        
        function visualize_waveguide(obj, z_range)           
            zq = z_range(1) : 0.05 : z_range(2);
            
            figure
            hold on
            grid on

            plot(zq, arrayfun(@(z) obj.wg1_lower_edge(z), zq), 'Color', 'blue')
            plot(zq, arrayfun(@(z) obj.wg1_upper_edge(z), zq), 'Color', 'blue')
            
            plot(zq, arrayfun(@(z) obj.wg2_lower_edge(z), zq), 'Color', 'red')
            plot(zq, arrayfun(@(z) obj.wg2_upper_edge(z), zq), 'Color', 'red')
            
            hold off
        end
    end
end



function val = modefcn(y, consts, a1, a2, d, k, eps, beta)
    gam = sqrt(eps*k^2 - beta^2);
    eta = sqrt(beta^2 - k^2);
    norm = 1;

    if (y <= 0)
        val = norm * consts(1) * exp(eta * y);
    elseif (0 < y && y <= a1)
        val = norm * (consts(2) * sin(gam*y) + consts(3) * cos(gam*y));
    elseif (a1 < y && y <= a1 + d)
        val = norm * (consts(4) * exp(eta*(y-a1)) + consts(5) * exp(-eta*(y-a1)));
    elseif (a1 + d < y && y <= a1 + d + a2)
        val = norm * (consts(6) * sin(gam*(y-a1-d)) + consts(7) * cos(gam*(y-a1-d)));
    elseif (a1 + d + a2 < y)
        val = norm * consts(8) * exp(-eta*(y-a1-a2-d));
    end
end

function err = eigenproblem(beta, k, eps, a1, a2, sep)
    gam = sqrt(eps*k^2 - beta^2);
    eta = sqrt(beta^2 - k^2);

    left_side = exp(-2*eta*sep) * (eta^2+gam^2)^2 * sin(gam*a1)*sin(gam*a2);

    right_side = (2*eta*gam*cos(gam*a1) + (eta^2-gam^2)*sin(gam*a1)) * (2*eta*gam*cos(gam*a2) + (eta^2-gam^2)*sin(gam*a2));

    err = left_side - right_side;
end