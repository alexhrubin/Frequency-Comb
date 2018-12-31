clear all; close all; clear classes; clc;
 
lam = 1;  % wavelength
%% Set flags.
inspect_only = false;
withuniformgrid=1;
source=2; % 1 for point source, 2 for modal source
% Geometric Parameters
lx=30;
ly=2;
air=1; %distance between wvg and the edge;
pmlthick = [0.2 0.2 0];  % PML thickness
dl = [0.01 0.01 0.01];  % grid size
dl_wvg = 0.005;
wg_rad = 0.25;

%rec = Rectangle(Axis.x, -15, [-lx/2 lx/2; (-ly/3-wg_rad/2) (-ly/3+wg_rad/2 + 0.0574)], 0.005);
domain = [-lx/2 lx/2; -ly/2-air ly/2+air; 0 dl(3)];
%wvg = Box([-lx/2 lx/2; (-ly/3-wg_rad/2)  (-ly/3+wg_rad/2) ; 0 dl(3)]);
wvg = Box([-lx/2 lx/2; -wg_rad  wg_rad ; 0 dl(3)]);
 
% Materials
vacuum = Material('vacuum', 'none', 1.0);
SiO2 = Material('SiO2', 'k', 12.0);
  
% Source
if (source==1)
    polarization = Axis.z;
    src = PointSrc(polarization, [-lx/2+0.4 -ly/3 dl(3)/2]);
elseif (source==2)
    normal_axis = Axis.x;
    intercept=-9.6;
    opts.clue = 'order';
    opts.order = 3;
    src = ModalSrc(normal_axis, intercept, opts);
end
% solution for uniform wvg.
[E, H, obj_array, src_array, extra,solveinfo] = maxwell_run(...
    'OSC', 1, lam, ...
    'DOM', vacuum, domain, dl, BC.p, pmlthick, ...
    'OBJ', SiO2, wvg,...
    'SRCJ', src, ...
    inspect_only);
 
%% Visualize the solution.
eps=solveinfo{4};
[E, H] = solve_eq_direct(solveinfo{1}, solveinfo{2}, solveinfo{3}, eps, solveinfo{5}, solveinfo{6}, solveinfo{7}, solveinfo{8}, solveinfo{9});
subplot(2,2,1);
pcolor(transpose(eps{1}));shading flat;colorbar
subplot(2,2,2);
pcolor(real(transpose(E{3})));shading flat;colorbar
colormap hot
 
% update eps for adiabatic taper;
eps=solveinfo{4};
epsnew;
[E2, H2] = solve_eq_direct(solveinfo{1}, solveinfo{2}, solveinfo{3}, eps, solveinfo{5}, solveinfo{6}, solveinfo{7}, solveinfo{8}, solveinfo{9});
 
subplot(2,2,3);
pcolor(transpose(eps{1}));shading flat;colorbar
subplot(2,2,4);
pcolor(real(transpose(E2{3})));shading flat;colorbar