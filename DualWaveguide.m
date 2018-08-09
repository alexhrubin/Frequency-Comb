classdef DualWaveguide
    properties
        sep     %separation of waveguide centers
        wg1     %lower waveguide (center is z-axis
        wg2     %upper waveguide (center is z = d
    end
    
    methods
        function obj = DualWaveguide(sep, wg1, wg2)
            obj.sep = sep;
            obj.wg1 = wg1;
            obj.wg2 = wg2;
        end
        
        function getbeta_position(obj, z, max_modes)
            disp('finish this function!!!')
        end
        
        function f = eigenmode_function(obj, z, beta)
            f = @(x) modefcn(x, obj.a1, a2, sep, k, eps, beta);
        end
        
        function v = boundary_conds(z, beta)
            gam = sqrt(obj.eps*obj.k^2 - beta^2);
            eta = sqrt(beta^2 - obj.k^2);
            
            a1  %what is the best way to relate these to z thru wg1,2????
            a2
            d
            
            M = [[1, 0, -1, 0, 0, 0, 0, 0];
                 [eta, -gam, 0, 0, 0, 0, 0, 0];
                 [0, sin(gam*a1), cos(gam*a1), -1, -1, 0, 0, 0];
                 [0, gam*cos(gam*a1), -gam*sin(gam*a1), -eta, eta, 0, 0, 0];
                 [0, 0, 0, exp(eta*d), exp(-eta*d), 0, -1, 0];
                 [0, 0, 0, eta*exp(eta*d), -eta*exp(-eta*d), -gam, 0, 0];
                 [0, 0, 0, 0, 0, sin(gam*a2), cos(gam*a2), -1];
                 [0, 0, 0, 0, 0, gam*cos(gam*a2), -gam*sin(gam*a2), eta]];
            
            v = null(M);
        end  
        
        function plot_eigenmodes(obj, mode_orders, z)
            disp('not completed')
        end
        
        function visualize_waveguide(obj, z_range)           
            zq = z_range(1) : 0.05 : z_range(2);
            
            figure
            hold on
            grid on

            plot(zq, obj.wg1.upper_edge_func(zq), 'Color', 'blue')
            plot(zq, obj.wg1.lower_edge_func(zq), 'Color', 'blue')
            
            plot(zq, obj.wg2.upper_edge_func(zq) + obj.sep, 'Color', 'red')
            plot(zq, obj.wg2.lower_edge_func(zq) + obj.sep, 'Color', 'red')
            
            hold off
        end
    end
end

function y = modefcn(x, a1, a2, sep, k, eps, beta)
    gam = sqrt(eps*k^2 - beta^2);
    eta = sqrt(beta^2 - k^2);
    
    if (x < 0)
        continue
    elseif (0 < x) && (x < a1)
        continue
    elseif (a1 < x) && (x < a1 + sep)
        continue
    elseif (a1 + sep < x) && (x < a1 + sep + a2)
        continue
    elseif (a1 + sep + a2 < x)
        continue
    end
end
        
        
        
        
        
        
        
        
        
        
        