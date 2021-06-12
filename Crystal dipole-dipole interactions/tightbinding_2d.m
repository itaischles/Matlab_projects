% For the theory behind this calculation see document stored in folder
% 'Theory'. To plot intermediate stages of calculation uncomment blocks
% marked with 'DEBUG'.

clearvars
clc

crystal_size = [100 100]; % [# unit cells in x, # unit cells in y] *** BOTH MUST BE EVEN!!! ***
Ax = 20; % primitive lattice spacing in x (Angstroms)
Ay = 20; % primitive lattice spacing in y (Angstroms)
dA = [0,0]; % vector pointing to location of dipole 'A' within unit cell (Angstroms)
dB = [10,10]; % vector pointing to location of dipole 'B' within unit cell (Angstroms)
muA = 8; % strength of dipole 'A' (Debye)
muB = 8; % strength of dipole 'B' (Debye)
epsA = 1242/550; % energy of A dipole (eV)
epsB = 1242/550; % energy of B dipole (eV)
thetaA = 45 * pi/180; % dipole 'A' orientation
thetaB = 100 * pi/180; % dipole 'B' orientation
eigenstate_to_plot = [0,0]; % plot only this eigenstate (values may range from -pi to pi)

input_params = {crystal_size, Ax, Ay, dA, dB, muA, muB, thetaA, thetaB};

% create dipole arrays
[XA,YA,XB,YB,muXA,muYA,muXB,muYB] = create_crystal_2d(input_params);

% calculate interaction matrices
[JAA,JBB,JAB] = calc_interaction_2d(XA,YA,XB,YB,dA,dB,muXA,muYA,muXB,muYB);

% calculate band structure
[KX,KY,Ubranch,Lbranch,JAA_k,JBB_k,JAB_k] = calc_band_structure_2d(Ax,Ay,JAA,JBB,JAB,epsA,epsB);

% calculate 1D path in First Brillouin Zone
[pathU,pathL,Nx,Ny] = calc_FBZ_path_2d(Ubranch,Lbranch);

% calculate eigenstates (eigen- wave functions)
[muU,muL,macro_dipole_U,macro_dipole_L] = calc_eigenstates_2d(KX,KY,Ax,Ay,XA,YA,dA,dB,JAA_k,JBB_k,JAB_k,1e-2,Ubranch,Lbranch,epsA,epsB,muXA,muYA,muXB,muYB,eigenstate_to_plot);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11); set(gcf, 'Position', [344 -223 455 420]); cla; hold on; axis equal
scatter(XA(:),YA(:),20,'r','filled');
scatter(XB(:),YB(:),20,'b','filled');
quiver(XA(:),YA(:), muXA(:), muYA(:), 0.3, 'Color', 'r');
quiver(XB(:),YB(:), muXB(:), muYB(:), 0.3, 'Color', 'b');
title('Dipole array'); xlabel('X (Angstrom)'); ylabel('Y (Angstrom)'); xlim(1.1*Ax*[-4,4]); ylim(1.1*Ay*[-4,4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure(22); set(gcf, 'Position', [72 162 560 420]); cla; hold on
% % % imagesc(JAA);
% % % title('J_{AA}'); axis equal
% % % figure(33); set(gcf, 'Position', [633 162 560 420]); cla; hold on
% % % imagesc(JBB);
% % % title('J_{BB}'); axis equal
% % % figure(44); set(gcf, 'Position', [1195 162 560 420]); cla; hold on
% % % imagesc(JAB);
% % % title('J_{AB}'); axis equal
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % figure(55); set(gcf, 'Position', [72 162 560 420]); cla; hold on
% % % imagesc(KX*Ax/pi,KY*Ay/pi,JAA_k);
% % % title('J_{AA}(k)'); axis equal; xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi')
% % % figure(66); set(gcf, 'Position', [633 162 560 420]); cla; hold on
% % % imagesc(KX*Ax/pi,KY*Ay/pi,JBB_k);
% % % title('J_{BB}(k)'); axis equal; xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi')
% % % figure(77); set(gcf, 'Position', [1195 162 560 420]); cla; hold on
% % % imagesc(KX*Ax/pi,KY*Ay/pi,abs(JAB_k));
% % % title('|J_{AB}(k)|'); axis equal; xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(88); set(gcf, 'Position', [354 50 560 420]); cla; hold on
surfc(KX*Ax/pi,KY*Ay/pi,Lbranch,'EdgeColor','none','FaceLighting','gouraud','SpecularStrength',0);
surfc(KX*Ax/pi,KY*Ay/pi,Ubranch,'EdgeColor','none','FaceLighting','gouraud','SpecularStrength',0);
title('Dispersion relation'); xlabel('k_xA_x/\pi'); ylabel('k_yA_y/\pi'); zlabel('Energy (eV)'); grid('on')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(99); set(gcf, 'Position', [1007 143 321 251]); hold on
ax = gca;
ci = ax.ColorOrderIndex;
colors = ax.ColorOrder;
plot(pathU,'Color',colors(ci,:));
plot(pathL,'Color',colors(ci,:));
plot(epsA*ones(size(pathU)),'Color',[0.5,0.5,0.5],'LineStyle','--')
plot(epsB*ones(size(pathU)),'Color',[0.5,0.5,0.5],'LineStyle','--')
xticks([1,...
        Nx/2+1,...
        Nx/2+Ny/2+1,...
        Nx/2+Ny/2+Nx/2+1,...
        Nx/2+Ny/2+Nx/2+Ny/2+1]);
xticklabels({'\Gamma','X','S','Y','\Gamma'})
title('Band diagrams'); ylabel('Energy (eV)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,kxind] = min(abs(KX*Ax-eigenstate_to_plot(1)));
[~,kyind] = min(abs(KY*Ay-eigenstate_to_plot(2)));

figure(111); set(gcf, 'Position', [395 1 560 420]); cla; hold on; axis equal
% % % scatter(XA(:),YA(:),10,'r','filled');
% % % scatter(XB(:),YB(:),10,'b','filled');
quiver(XA(:),YA(:), muU{1}(:), muU{2}(:), 0.3, 'Color', 'r');
quiver(XB(:),YB(:), muU{3}(:), muU{4}(:), 0.3, 'Color', 'b');
quiver(0, 0, macro_dipole_U(1), macro_dipole_U(2), 10, 'Color', 'm', 'LineWidth', 3);
title(strcat('KA=[',num2str(KX(kxind)*Ax),',',num2str(KY(kyind)*Ay),']','   E=',num2str(Ubranch(kyind,kxind)))); xlabel('X (Angstrom)'); ylabel('Y (Angstrom)');

figure(222); set(gcf, 'Position', [956 1 560 420]); cla; hold on; axis equal
% % % scatter(XA(:),YA(:),10,'r','filled');
% % % scatter(XB(:),YB(:),10,'b','filled');
quiver(XA(:),YA(:), muL{1}(:), muL{2}(:), 0.3, 'Color', 'r');
quiver(XB(:),YB(:), muL{3}(:), muL{4}(:), 0.3, 'Color', 'b');
quiver(0, 0, macro_dipole_L(1), macro_dipole_L(2), 10, 'Color', 'm', 'LineWidth', 3);
title(strcat('KA=[',num2str(KX(kxind)*Ax),',',num2str(KY(kyind)*Ay),']','   E=',num2str(Lbranch(kyind,kxind)))); xlabel('X (Angstrom)'); ylabel('Y (Angstrom)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
