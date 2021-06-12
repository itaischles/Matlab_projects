% For the theory behind this calculation see document stored in folder
% 'Theory'. To plot intermediate stages of calculation uncomment blocks
% marked with 'DEBUG'.

clearvars
clc

crystal_size = [20 20 20]; % [# unit cells in x, # unit cells in y, # unit cells in z] *** ALL MUST BE EVEN!!! ***
Ax = 10; % primitive lattice spacing in x (Angstroms)
Ay = 10; % primitive lattice spacing in y (Angstroms)
Az = 10; % primitive lattice spacing in z (Angstroms)
dA = [0,0,0]; % vector pointing to location of dipole 'A' within unit cell (Angstroms)
dB = [Ax/2,0,0]; % vector pointing to location of dipole 'B' within unit cell (Angstroms)
muA = 1; % strength of dipole 'A' (Debye)
muB = 1; % strength of dipole 'B' (Debye)
epsA = 1242/530; % energy of A dipole (eV)
epsB = 1242/530; % energy of B dipole (eV)
angleA = [80,45] * pi/180; % dipole 'A' orientation (phi (x,y) plane, theta (z axis))
angleB = [80,45] * pi/180; % dipole 'B' orientation (phi (x,y) plane, theta (z axis))
dielconst = 3; % dielectric constant
eigenstate_to_plot = [0,0,0]; % plot only this eigenstate (using indices, e.g. [2,0,0] is kx=2*delta(kx), ky,kz=0)

% calculate subset for plotting purposes
rngX = crystal_size(1)/2-3:1:crystal_size(1)/2+5;
rngY = crystal_size(2)/2-3:1:crystal_size(2)/2+5;
rngZ = crystal_size(3)/2-3:1:crystal_size(3)/2+5;
plot_subset = {rngX,rngY,rngZ};

input_params = {crystal_size, Ax, Ay, Az, dA, dB, muA, muB, angleA, angleB};

% create dipole arrays
[XA,YA,ZA,XB,YB,ZB,muXA,muYA,muZA,muXB,muYB,muZB] = create_crystal_3d(input_params);

% calculate interaction matrices
[JAA,JBB,JAB] = calc_interaction_3d(XA,YA,ZA,XB,YB,ZB,dA,dB,muXA,muYA,muZA,muXB,muYB,muZB,dielconst);

% calculate band structure
[KX,KY,KZ,Ubranch,Lbranch,JAA_k,JBB_k,JAB_k] = calc_band_structure_3d(XA,YA,ZA,dA,Ax,Ay,Az,JAA,JBB,JAB,epsA,epsB);

% calculate 1D path in First Brillouin Zone
[pathU,pathL,Nx,Ny,Nz] = calc_FBZ_path_3d(Ubranch,Lbranch);

% calculate eigenstates (eigen- wave functions)
[muU,muL,macro_dipole_U,macro_dipole_L] = calc_eigenstates_3d(KX,KY,KZ,XA,YA,ZA,dA,dB,JAA_k,JBB_k,JAB_k,1e-2,Ubranch,Lbranch,epsA,epsB,muXA,muYA,muZA,muXB,muYB,muZB,eigenstate_to_plot,plot_subset);

%% Crystal setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_params_plotting = {[6,6,6], Ax, Ay, Az, dA, dB, muA, muB, angleA, angleB};
[XAplot,YAplot,ZAplot,XBplot,YBplot,ZBplot,muXAplot,muYAplot,muZAplot,muXBplot,muYBplot,muZBplot] = create_crystal_3d(input_params_plotting);
figure(11); set(gcf, 'Position', [344 -223 455 420]); cla; hold on; axis equal; view(45,20); camproj('perspective')
scatter3(XAplot(:),YAplot(:),ZAplot(:),10,'r','filled');
scatter3(XBplot(:),YBplot(:),ZBplot(:),10,'b','filled');
quiver3(XAplot(:),YAplot(:),ZAplot(:),muXAplot(:),muYAplot(:),muZAplot(:),0.5,'Color','r');
quiver3(XBplot(:),YBplot(:),ZBplot(:),muXBplot(:),muYBplot(:),muZBplot(:),0.5,'Color','b');
title('Dipole array'); xlabel('X'); ylabel('Y'); zlabel('Z'); xlim(1.1*Ax*[-2,2]); ylim(1.1*Ay*[-2,2]); zlim(1.1*Az*[-2,2]); set(gca,'Box','on')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Interaction matrices

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure(22); set(gcf, 'Position', [72 162 560 420]); cla; hold on
% % % imagesc(JAA(:,:,crystal_size(3)/2+1))
% % % % % % slc = slice(XA-dA(1),YA-dA(2),ZA-dA(3),JAA,0,0,0);
% % % % % % title('J_{AA}'); axis equal; set(slc,'EdgeColor','none'); clear slc
% % % 
% % % figure(33); set(gcf, 'Position', [633 162 560 420]); cla; hold on
% % % imagesc(JBB(:,:,crystal_size(3)/2+1))
% % % % % % slc = slice(XA-dA(1),YA-dA(2),ZA-dA(3),JBB,0,0,0);
% % % % % % title('J_{BB}'); axis equal; set(slc,'EdgeColor','none'); clear slc
% % % 
% % % figure(44); set(gcf, 'Position', [1195 162 560 420]); cla; hold on
% % % imagesc(JAB(:,:,crystal_size(3)/2+1))
% % % % % % slc = slice(XA-dA(1),YA-dA(2),ZA-dA(3),JAB,0,0,0);
% % % % % % title('J_{AB}'); axis equal; set(slc,'EdgeColor','none'); clear slc
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fourier transform of interaction matrices

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure(55); set(gcf, 'Position', [72 162 560 420]); cla; hold on
% % % % % % slc = slice(KX*Ax/pi,KY*Ay/pi,KZ*Az/pi,JAA_k,0,0,0);
% % % % % % title('J_{AA}(k)'); axis equal; xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi'); zlabel('k_zA_z/\pi'); set(slc,'EdgeColor','none'); clear slc
% % % imagesc(KX*Ax/pi,KY*Ay/pi,JAA_k(:,:,crystal_size(3)/2+1)); xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi')
% % % title('J_{AA}(k)'); axis equal;
% % % figure(66); set(gcf, 'Position', [633 162 560 420]); cla; hold on
% % % % % % slc = slice(KX*Ax/pi,KY*Ay/pi,KZ*Az/pi,JBB_k,0,0,0);
% % % % % % title('J_{BB}(k)'); axis equal; xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi'); zlabel('k_zA_z/\pi'); set(slc,'EdgeColor','none'); clear slc
% % % imagesc(KX*Ax/pi,KY*Ay/pi,JBB_k(:,:,crystal_size(3)/2+1)); xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi')
% % % title('J_{BB}(k)'); axis equal;
% % % figure(77); set(gcf, 'Position', [1195 162 560 420]); cla; hold on
% % % % % % slc = slice(KX*Ax/pi,KY*Ay/pi,KZ*Az/pi,abs(JAB_k),0,0,0);
% % % % % % title('|J_{AB}(k)|'); axis equal; xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi'); zlabel('k_zA_z/\pi'); set(slc,'EdgeColor','none'); clear slc
% % % imagesc(KX*Ax/pi,KY*Ay/pi,abs(JAB_k(:,:,crystal_size(3)/2+1))); xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi')
% % % title('|J_{AB}(k)|'); axis equal;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First Briloine zone surfaces

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure(88); set(gcf, 'Position', [317 73 560 420]); cla; hold on; axis equal; view(40,35);
% % % pU = slice(KX*Ax/pi,KY*Ay/pi,KZ*Az/pi,Ubranch,0,0,0);
% % % title('Dispersion relation, U branch'); xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi'); zlabel('k_zA_z/\pi'); set(pU,'EdgeColor','none')
% % % 
% % % figure(99); set(gcf, 'Position', [875 70 560 420]); cla; hold on; axis equal; view(40,35)
% % % pL = slice(KX*Ax/pi,KY*Ay/pi,KZ*Az/pi,Lbranch,0,0,0);
% % % title('Dispersion relation, L branch'); xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi'); zlabel('k_zA_z/\pi'); set(pL,'EdgeColor','none')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Brilloine zone 1D path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(111); set(gcf, 'Position', [1442 69 454 363]); hold on
ax = gca;
pathU_plot = plot(pathU);
plot(pathL,'Color',pathU_plot.Color);
clear pathU_plot
plot(epsA*ones(size(pathU)),'Color',[0.5,0.5,0.5],'LineStyle','--')
plot(epsB*ones(size(pathU)),'Color',[0.5,0.5,0.5],'LineStyle','--')
xticks([1,...
        Nx/2+1,...
        Nx/2+Ny/2+1,...
        Nx/2+Ny/2+Nx/2+1,...
        Nx/2+Ny/2+Nx/2+Ny/2+1,...
        Nx/2+Ny/2+Nx/2+Ny/2+Nz/2+1,...
        Nx/2+Ny/2+Nx/2+Ny/2+Nz/2+Nx/2+1,...
        Nx/2+Ny/2+Nx/2+Ny/2+Nz/2+Nx/2+Ny/2+1,...
        Nx/2+Ny/2+Nx/2+Ny/2+Nz/2+Nx/2+Ny/2+Nx/2+1,...
        Nx/2+Ny/2+Nx/2+Ny/2+Nz/2+Nx/2+Ny/2+Nx/2+Ny/2+1]);
xticklabels({'\Gamma','X','S','Y','\Gamma','Z','U','R','T','Z'})
title('Band diagrams'); ylabel('Energy (a.u.)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eigenstate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kxind = eigenstate_to_plot(1) + numel(KX)/2+1;
kyind = eigenstate_to_plot(2) + numel(KY)/2+1;
kzind = eigenstate_to_plot(3) + numel(KZ)/2+1;

XA_subset = XA(rngX,rngY,rngZ);
YA_subset = YA(rngX,rngY,rngZ);
ZA_subset = ZA(rngX,rngY,rngZ);
XB_subset = XB(rngX,rngY,rngZ);
YB_subset = YB(rngX,rngY,rngZ);
ZB_subset = ZB(rngX,rngY,rngZ);

figure(222); set(gcf, 'Position', [395 1 560 420]); cla; hold on; axis equal; camproj('orthographic')
scatter3(XA_subset(:),YA_subset(:),ZA_subset(:),10,'r','filled');
scatter3(XB_subset(:),YB_subset(:),ZB_subset(:),10,'b','filled');
quiver3(XA_subset(:),YA_subset(:),ZA_subset(:), muU{1}(:), muU{2}(:), muU{3}(:), 1, 'Color', 'r');
quiver3(XB_subset(:),YB_subset(:),ZB_subset(:), muU{4}(:), muU{5}(:), muU{6}(:), 1, 'Color', 'b');
quiver3(0, 0, 0, macro_dipole_U(1), macro_dipole_U(2), macro_dipole_U(3), 10, 'Color', 'm', 'LineWidth', 3);
title(strcat('k=[',num2str(KX(kxind)),',',num2str(KY(kyind)),',',num2str(KZ(kzind)),']','   E=',num2str(Ubranch(kyind,kxind,kzind)))); xlabel('X'); ylabel('Y'); zlabel('Z');

figure(333); set(gcf, 'Position', [956 1 560 420]); cla; hold on; axis equal; camproj('orthographic')
scatter3(XA_subset(:),YA_subset(:),ZA_subset(:),10,'r','filled');
scatter3(XB_subset(:),YB_subset(:),ZB_subset(:),10,'b','filled');
quiver3(XA_subset(:),YA_subset(:),ZA_subset(:), muL{1}(:), muL{2}(:), muL{3}(:), 1, 'Color', 'r');
quiver3(XB_subset(:),YB_subset(:),ZB_subset(:), muL{4}(:), muL{5}(:), muL{6}(:), 1, 'Color', 'b');
quiver3(0, 0, 0, macro_dipole_L(1), macro_dipole_L(2), macro_dipole_L(3), 10, 'Color', 'm', 'LineWidth', 3);
title(strcat('k=[',num2str(KX(kxind)),',',num2str(KY(kyind)),',',num2str(KZ(kzind)),']','   E=',num2str(Lbranch(kyind,kxind,kzind)))); xlabel('X'); ylabel('Y'); zlabel('Z');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Macroscopic dipole moment

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % abs_macro_dipole_U = sqrt(macro_dipole_U(1).^2+macro_dipole_U(2).^2+macro_dipole_U(3).^2);
% % % abs_macro_dipole_L = sqrt(macro_dipole_L(1).^2+macro_dipole_L(2).^2+macro_dipole_L(3).^2);
% % % disp('*******************************************************************************')
% % % disp(strcat('Norm. macroscopic dipole moment of selected U-eigenstate U branch: ',num2str(abs_macro_dipole_U)))
% % % disp(strcat('Norm. macroscopic dipole moment of selected L-eigenstate U branch: ',num2str(abs_macro_dipole_L)))
% % % disp('*******************************************************************************')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%