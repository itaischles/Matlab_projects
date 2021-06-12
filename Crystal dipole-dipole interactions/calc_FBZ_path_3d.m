function [pathU,pathL,Nx,Ny,Nz] = calc_FBZ_path_3d(Ubranch,Lbranch)

tic

Nx = numel(Ubranch(1,:,1));
Ny = numel(Ubranch(:,1,1));
Nz = numel(Ubranch(1,1,:));

% % % gamma_2_x_U = Ubranch(Ny/2+1, Nx/2+1:1:end, Nz/2+1);
% % % x_2_s_U = Ubranch(Ny/2+1+1:1:end, end, Nz/2+1);
% % % s_2_y_U = Ubranch(end, end-1:-1:Nx/2+1, Nz/2+1);
% % % y_2_gamma_U = Ubranch(end-1:-1:Ny/2+1, Nx/2+1, Nz/2+1);
% % % gamma_2_z_U = Ubranch(Ny/2+1, Nx/2+1, Nz/2+1+1:1:end);
% % % z_2_u_U = Ubranch(Ny/2+1, Nx/2+1+1:1:end, end);
% % % u_2_r_U = Ubranch(Ny/2+1+1:1:end, end, end);
% % % r_2_t_U = Ubranch(end, end-1:-1:Nx/2+1, end);
% % % t_2_z_U = Ubranch(end-1:-1:Ny/2+1, Nx/2+1, end);

gamma_2_x_U = Ubranch(Ny/2+1, Nx/2+1:-1:1, Nz/2+1);
x_2_s_U = Ubranch(Ny/2+1-1:-1:1, 1, Nz/2+1);
s_2_y_U = Ubranch(1, 2:1:Nx/2+1, Nz/2+1);
y_2_gamma_U = Ubranch(2:1:Ny/2+1, Nx/2+1, Nz/2+1);
gamma_2_z_U = Ubranch(Ny/2+1, Nx/2+1, Nz/2+1-1:-1:1);
z_2_u_U = Ubranch(Ny/2+1, Nx/2+1-1:-1:1, 1);
u_2_r_U = Ubranch(Ny/2+1-1:-1:1, 1, 1);
r_2_t_U = Ubranch(1, 2:1:Nx/2+1, 1);
t_2_z_U = Ubranch(2:1:Ny/2+1, Nx/2+1, 1);

% % % gamma_2_x_U = Ubranch(1, 1:1:end/2+1, 1);
% % % x_2_s_U = Ubranch(2:1:end/2+1, end/2+1, 1);
% % % s_2_y_U = Ubranch(end/2+1, end/2+1-1:-1:1, 1);
% % % y_2_gamma_U = Ubranch(end/2+1-1:-1:1, 1, 1);
% % % gamma_2_z_U = Ubranch(1, 1, 2:1:end/2+1);
% % % z_2_u_U = Ubranch(1, 2:1:end/2+1, end/2+1);
% % % u_2_r_U = Ubranch(2:1:end/2+1, end/2+1, end/2+1);
% % % r_2_t_U = Ubranch(end/2+1, end/2+1-1:-1:1, end/2+1);
% % % t_2_z_U = Ubranch(end/2+1-1:-1:1, 1, end/2+1);

gamma_2_x_U = reshape(gamma_2_x_U, 1,[]);
x_2_s_U = reshape(x_2_s_U, 1,[]);
s_2_y_U = reshape(s_2_y_U, 1,[]);
y_2_gamma_U = reshape(y_2_gamma_U, 1,[]);
gamma_2_z_U = reshape(gamma_2_z_U, 1,[]);
z_2_u_U = reshape(z_2_u_U, 1,[]);
u_2_r_U = reshape(u_2_r_U, 1,[]);
r_2_t_U = reshape(r_2_t_U, 1,[]);
t_2_z_U = reshape(t_2_z_U, 1,[]);

% % % gamma_2_x_L = Lbranch(Ny/2+1, Nx/2+1:1:end, Nz/2+1);
% % % x_2_s_L = Lbranch(Ny/2+1+1:1:end, end, Nz/2+1);
% % % s_2_y_L = Lbranch(end, end-1:-1:Nx/2+1, Nz/2+1);
% % % y_2_gamma_L = Lbranch(end-1:-1:Ny/2+1, Nx/2+1, Nz/2+1);
% % % gamma_2_z_L = Lbranch(Ny/2+1, Nx/2+1, Nz/2+1+1:1:end);
% % % z_2_u_L = Lbranch(Ny/2+1, Nx/2+1+1:1:end, end);
% % % u_2_r_L = Lbranch(Ny/2+1+1:1:end, end, end);
% % % r_2_t_L = Lbranch(end, end-1:-1:Nx/2+1, end);
% % % t_2_z_L = Lbranch(end-1:-1:Ny/2+1, Nx/2+1, end);

gamma_2_x_L = Lbranch(Ny/2+1, Nx/2+1:-1:1, Nz/2+1);
x_2_s_L = Lbranch(Ny/2+1-1:-1:1, 1, Nz/2+1);
s_2_y_L = Lbranch(1, 2:1:Nx/2+1, Nz/2+1);
y_2_gamma_L = Lbranch(2:1:Ny/2+1, Nx/2+1, Nz/2+1);
gamma_2_z_L = Lbranch(Ny/2+1, Nx/2+1, Nz/2+1-1:-1:1);
z_2_u_L = Lbranch(Ny/2+1, Nx/2+1-1:-1:1, 1);
u_2_r_L = Lbranch(Ny/2+1-1:-1:1, 1, 1);
r_2_t_L = Lbranch(1, 2:1:Nx/2+1, 1);
t_2_z_L = Lbranch(2:1:Ny/2+1, Nx/2+1, 1);

% % % gamma_2_x_L = Lbranch(1, 1:1:end/2+1, 1);
% % % x_2_s_L = Lbranch(2:1:end/2+1, end/2+1, 1);
% % % s_2_y_L = Lbranch(end/2+1, end/2+1-1:-1:1, 1);
% % % y_2_gamma_L = Lbranch(end/2+1-1:-1:1, 1, 1);
% % % gamma_2_z_L = Lbranch(1, 1, 2:1:end/2+1);
% % % z_2_u_L = Lbranch(1, 2:1:end/2+1, end/2+1);
% % % u_2_r_L = Lbranch(2:1:end/2+1, end/2+1, end/2+1);
% % % r_2_t_L = Lbranch(end/2+1, end/2+1-1:-1:1, end/2+1);
% % % t_2_z_L = Lbranch(end/2+1-1:-1:1, 1, end/2+1);

gamma_2_x_L = reshape(gamma_2_x_L, 1,[]);
x_2_s_L = reshape(x_2_s_L, 1,[]);
s_2_y_L = reshape(s_2_y_L, 1,[]);
y_2_gamma_L = reshape(y_2_gamma_L, 1,[]);
gamma_2_z_L = reshape(gamma_2_z_L, 1,[]);
z_2_u_L = reshape(z_2_u_L, 1,[]);
u_2_r_L = reshape(u_2_r_L, 1,[]);
r_2_t_L = reshape(r_2_t_L, 1,[]);
t_2_z_L = reshape(t_2_z_L, 1,[]);

pathU = [gamma_2_x_U, x_2_s_U, s_2_y_U, y_2_gamma_U, gamma_2_z_U, z_2_u_U, u_2_r_U, r_2_t_U, t_2_z_U];
pathL = [gamma_2_x_L, x_2_s_L, s_2_y_L, y_2_gamma_L, gamma_2_z_L, z_2_u_L, u_2_r_L, r_2_t_L, t_2_z_L];

elapsedTime = toc;
disp(strcat("FBZ path calculation completed in ", num2str(elapsedTime), " seconds"));

end

