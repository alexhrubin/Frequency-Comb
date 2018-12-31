y0=0.5/0.01;
y1=2.9/0.01;
[Nx,Ny]=size(eps{1});

% N_center_bottom = round(7*(Ny-200)/10+100);
% N_center_top = round(7*(Ny-200)/10 + 100);
N_center = round(Ny/2);

%r_y = round(Ny/2.23)-0.22/2;
%rectangle = Rectangle(Axis.x, -15, [0 Nx; r_y r_y+45], 0.005);

eps{1} = ones(Nx, Ny);
eps{2} = ones(Nx, Ny);
eps{3} = ones(Nx, Ny);

for ii=1:Nx
%     jj_bottom = round(ii/Nx * (y1-y0) + y0);
%     jj_top = round(ii/Nx * (y0-y1) + y1);

    jj = round(ii/Nx * (y1-y0) + y0);
    eps{1}(ii, N_center-round(jj/2) : N_center+round(jj/2))=12;
    eps{2}(ii, N_center-round(jj/2) : N_center+round(jj/2))=12;
    eps{3}(ii, N_center-round(jj/2) : N_center+round(jj/2))=12;

    
    
%     eps{1}(ii,N_center_top-round(jj_top/2) : N_center_top+round(jj_top/2))=12;
%     eps{2}(ii,N_center_top-round(jj_top/2) : N_center_top+round(jj_top/2))=12;
%     eps{3}(ii,N_center_top-round(jj_top/2) : N_center_top+round(jj_top/2))=12;
    
%     eps{1}(ii,N_center_bottom-round(jj_bottom/2) : N_center_bottom+round(jj_bottom/2))=12;
%     eps{2}(ii,N_center_bottom-round(jj_bottom/2) : N_center_bottom+round(jj_bottom/2))=12;
%     eps{3}(ii,N_center_bottom-round(jj_bottom/2) : N_center_bottom+round(jj_bottom/2))=12;
    
end

% for ii = 1:Nx
%      jj=exp(a*Nx)+y0;
%      eps{1}(ii,Nyo-round(exp(a*jj/d)):Nyo+round(exp(a*jj/d))) = 12;
%      eps{2}(ii,Nyo-round(exp(a*jj/d)):Nyo+round(exp(a*jj/d))) = 12;
%      eps{3}(ii,Nyo-round(exp(a*jj/d)):Nyo+round(exp(a*jj/d))) = 12;
% end