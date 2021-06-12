function [pathU,pathL,Nx,Ny] = calc_FBZ_path_2d(Ubranch,Lbranch)

Nx = numel(Ubranch(1,:));
Ny = numel(Ubranch(:,1));

% % % gamma_2_x_U = Ubranch(Ny/2+1, Nx/2+1:1:end);
% % % x_2_s_U = Ubranch(Ny/2+2:1:end, end)';
% % % s_2_y_U = Ubranch(end, end-1:-1:Nx/2+1);
% % % y_2_gamma_U = Ubranch(end-1:-1:Ny/2+1, Nx/2+1)';

gamma_2_x_U = Ubranch(Ny/2+1, Nx/2+1:-1:1);
x_2_s_U = Ubranch(Ny/2+1-1:-1:1, 1);
s_2_y_U = Ubranch(1, 2:1:Nx/2+1);
y_2_gamma_U = Ubranch(2:1:Ny/2+1, Nx/2+1);

gamma_2_x_U = reshape(gamma_2_x_U, 1,[]);
x_2_s_U = reshape(x_2_s_U, 1,[]);
s_2_y_U = reshape(s_2_y_U, 1,[]);
y_2_gamma_U = reshape(y_2_gamma_U, 1,[]);

% % % gamma_2_x_L = Lbranch(Ny/2+1, Nx/2+1:1:end);
% % % x_2_s_L = Lbranch(Ny/2+2:1:end, end)';
% % % s_2_y_L = Lbranch(end, end-1:-1:Nx/2+1);
% % % y_2_gamma_L = Lbranch(end-1:-1:Ny/2+1, Nx/2+1)';

gamma_2_x_L = Lbranch(Ny/2+1, Nx/2+1:-1:1);
x_2_s_L = Lbranch(Ny/2+1-1:-1:1, 1);
s_2_y_L = Lbranch(1, 2:1:Nx/2+1);
y_2_gamma_L = Lbranch(2:1:Ny/2+1, Nx/2+1);

gamma_2_x_L = reshape(gamma_2_x_L, 1,[]);
x_2_s_L = reshape(x_2_s_L, 1,[]);
s_2_y_L = reshape(s_2_y_L, 1,[]);
y_2_gamma_L = reshape(y_2_gamma_L, 1,[]);

pathU = [gamma_2_x_U, x_2_s_U, s_2_y_U, y_2_gamma_U];
pathL = [gamma_2_x_L, x_2_s_L, s_2_y_L, y_2_gamma_L];

end

