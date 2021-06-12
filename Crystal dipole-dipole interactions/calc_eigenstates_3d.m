function [muU,muL,macro_dipole_U,macro_dipole_L] = calc_eigenstates_3d(KX,KY,KZ,XA,YA,ZA,dA,dB,JAA_k,JBB_k,JAB_k,JAB_cutoff,Ubranch,Lbranch,epsA,epsB,muXA,muYA,muZA,muXB,muYB,muZB,eigenstate_to_plot,subset)

tic

delta_d = dB-dA;

% % % muXAU = cell(numel(KY),numel(KX),numel(KZ));
% % % muYAU = cell(numel(KY),numel(KX),numel(KZ));
% % % muZAU = cell(numel(KY),numel(KX),numel(KZ));
% % % muXBU = cell(numel(KY),numel(KX),numel(KZ));
% % % muYBU = cell(numel(KY),numel(KX),numel(KZ));
% % % muZBU = cell(numel(KY),numel(KX),numel(KZ));
% % % 
% % % muXAL = cell(numel(KY),numel(KX),numel(KZ));
% % % muYAL = cell(numel(KY),numel(KX),numel(KZ));
% % % muZAL = cell(numel(KY),numel(KX),numel(KZ));
% % % muXBL = cell(numel(KY),numel(KX),numel(KZ));
% % % muYBL = cell(numel(KY),numel(KX),numel(KZ));
% % % muZBL = cell(numel(KY),numel(KX),numel(KZ));

% % % macro_dipole_X_U = zeros(numel(KY),numel(KX),numel(KZ));
% % % macro_dipole_Y_U = zeros(numel(KY),numel(KX),numel(KZ));
% % % macro_dipole_Z_U = zeros(numel(KY),numel(KX),numel(KZ));
% % % macro_dipole_X_L = zeros(numel(KY),numel(KX),numel(KZ));
% % % macro_dipole_Y_L = zeros(numel(KY),numel(KX),numel(KZ));
% % % macro_dipole_Z_L = zeros(numel(KY),numel(KX),numel(KZ));

if ~isempty(eigenstate_to_plot) && numel(eigenstate_to_plot)==3
    
    i = eigenstate_to_plot(1) + numel(KX)/2+1;
    j = eigenstate_to_plot(2) + numel(KY)/2+1;
    k = eigenstate_to_plot(3) + numel(KZ)/2+1;
    
    if abs(JAB_k(j,i,k))>JAB_cutoff
        cA_k_U = 1;
        cB_k_U = (Ubranch(j,i,k)-epsA-JAA_k(j,i,k))/(JAB_k(j,i,k)*exp(1i*(KX(i)*delta_d(1)+KY(j)*delta_d(2)+KZ(k)*delta_d(3))));
        cA_k_L = 1;
        cB_k_L = (Lbranch(j,i,k)-epsA-JAA_k(j,i,k))/(JAB_k(j,i,k)*exp(1i*(KX(i)*delta_d(1)+KY(j)*delta_d(2)+KZ(k)*delta_d(3))));
    else
        if JAA_k(j,i,k)+epsA>JBB_k(j,i,k)+epsB
            cA_k_U = 1;
            cB_k_U = 0;
            cA_k_L = 0;
            cB_k_L = 1;
        else
            cA_k_U = 0;
            cB_k_U = 1;
            cA_k_L = 1;
            cB_k_L = 0;
        end
    end
    
    phase = exp(1i*(KX(i)*(XA-dA(1))+KY(j)*(YA-dA(2))+KZ(k)*(ZA-dA(3))));
    
    eigenstates_U_A = phase.*cA_k_U*exp(1i*(KX(i)*dA(1)+KY(j)*dA(2)+KZ(k)*dA(3)));
    eigenstates_L_A = phase.*cA_k_L*exp(1i*(KX(i)*dA(1)+KY(j)*dA(2)+KZ(k)*dA(3)));
    eigenstates_U_B = phase.*cB_k_U*exp(1i*(KX(i)*dA(1)+KY(j)*dB(2)+KZ(k)*dB(3)));
    eigenstates_L_B = phase.*cB_k_L*exp(1i*(KX(i)*dA(1)+KY(j)*dB(2)+KZ(k)*dB(3)));
    
    muXAU = real(eigenstates_U_A).*muXA;
    muYAU = real(eigenstates_U_A).*muYA;
    muZAU = real(eigenstates_U_A).*muZA;
    muXBU = real(eigenstates_U_B).*muXB;
    muYBU = real(eigenstates_U_B).*muYB;
    muZBU = real(eigenstates_U_B).*muZB;
    
    muXAL = real(eigenstates_L_A).*muXA;
    muYAL = real(eigenstates_L_A).*muYA;
    muZAL = real(eigenstates_L_A).*muZA;
    muXBL = real(eigenstates_L_B).*muXB;
    muYBL = real(eigenstates_L_B).*muYB;
    muZBL = real(eigenstates_L_B).*muZB;
    
    macro_dipole_X_U = sum(sum(sum(muXAU+muXBU)))/numel(XA);
    macro_dipole_Y_U = sum(sum(sum(muYAU+muYBU)))/numel(XA);
    macro_dipole_Z_U = sum(sum(sum(muZAU+muZBU)))/numel(XA);
    macro_dipole_X_L = sum(sum(sum(muXAL+muXBL)))/numel(XA);
    macro_dipole_Y_L = sum(sum(sum(muYAL+muYBL)))/numel(XA);
    macro_dipole_Z_L = sum(sum(sum(muZAL+muZBL)))/numel(XA);
    
end
% % %     
% % %     for i=1:numel(KX)
% % %         for j=1:numel(KY)
% % %             for k=1:numel(KZ)
% % %                 
% % %                 if abs(JAB_k(j,i,k))>JAB_cutoff
% % %                     cA_k_U = 1;
% % %                     cB_k_U = (Ubranch(j,i,k)-epsA-JAA_k(j,i,k))/(JAB_k(j,i,k)*exp(1i*(KX(i)*delta_d(1)+KY(j)*delta_d(2)+KZ(k)*delta_d(3))));
% % %                     cA_k_L = 1;
% % %                     cB_k_L = (Lbranch(j,i,k)-epsA-JAA_k(j,i,k))/(JAB_k(j,i,k)*exp(1i*(KX(i)*delta_d(1)+KY(j)*delta_d(2)+KZ(k)*delta_d(3))));
% % %                 else
% % %                     if JAA_k(j,i,k)+epsA>JBB_k(j,i,k)+epsB
% % %                         cA_k_U = 1;
% % %                         cB_k_U = 0;
% % %                         cA_k_L = 0;
% % %                         cB_k_L = 1;
% % %                     else
% % %                         cA_k_U = 0;
% % %                         cB_k_U = 1;
% % %                         cA_k_L = 1;
% % %                         cB_k_L = 0;
% % %                     end
% % %                 end
% % %                 
% % %                 phase = exp(1i*(KX(i)*(XA-dA(1))+KY(j)*(YA-dA(2))+KZ(k)*(ZA-dA(3))));
% % %                 
% % %                 eigenstates_U_A = phase.*cA_k_U*exp(1i*(KX(i)*dA(1)+KY(j)*dA(2)+KZ(k)*dA(3)));
% % %                 eigenstates_L_A = phase.*cA_k_L*exp(1i*(KX(i)*dA(1)+KY(j)*dA(2)+KZ(k)*dA(3)));
% % %                 eigenstates_U_B = phase.*cB_k_U*exp(1i*(KX(i)*dA(1)+KY(j)*dB(2)+KZ(k)*dB(3)));
% % %                 eigenstates_L_B = phase.*cB_k_L*exp(1i*(KX(i)*dA(1)+KY(j)*dB(2)+KZ(k)*dB(3)));
% % %                 
% % %                 muXAU{j,i,k} = real(eigenstates_U_A).*muXA;
% % %                 muYAU{j,i,k} = real(eigenstates_U_A).*muYA;
% % %                 muZAU{j,i,k} = real(eigenstates_U_A).*muZA;
% % %                 muXBU{j,i,k} = real(eigenstates_U_B).*muXB;
% % %                 muYBU{j,i,k} = real(eigenstates_U_B).*muYB;
% % %                 muZBU{j,i,k} = real(eigenstates_U_B).*muZB;
% % %                 
% % %                 muXAL{j,i,k} = real(eigenstates_L_A).*muXA;
% % %                 muYAL{j,i,k} = real(eigenstates_L_A).*muYA;
% % %                 muZAL{j,i,k} = real(eigenstates_L_A).*muZA;
% % %                 muXBL{j,i,k} = real(eigenstates_L_B).*muXB;
% % %                 muYBL{j,i,k} = real(eigenstates_L_B).*muYB;
% % %                 muZBL{j,i,k} = real(eigenstates_L_B).*muZB;
% % %                 
% % %                 macro_dipole_X_U(j,i,k) = sum(sum(sum(muXAU{j,i,k}(:)+muXBU{j,i,k}(:))))/numel(XA);
% % %                 macro_dipole_Y_U(j,i,k) = sum(sum(sum(muYAU{j,i,k}(:)+muYBU{j,i,k}(:))))/numel(XA);
% % %                 macro_dipole_Z_U(j,i,k) = sum(sum(sum(muZAU{j,i,k}(:)+muZBU{j,i,k}(:))))/numel(XA);
% % %                 macro_dipole_X_L(j,i,k) = sum(sum(sum(muXAL{j,i,k}(:)+muXBL{j,i,k}(:))))/numel(XA);
% % %                 macro_dipole_Y_L(j,i,k) = sum(sum(sum(muYAL{j,i,k}(:)+muYBL{j,i,k}(:))))/numel(XA);
% % %                 macro_dipole_Z_L(j,i,k) = sum(sum(sum(muZAL{j,i,k}(:)+muZBL{j,i,k}(:))))/numel(XA);
% % %                 
% % %             end
% % %         end
% % %     end
% % % end

muU = {muXAU(subset{1},subset{2},subset{3}),muYAU(subset{1},subset{2},subset{3}),muZAU(subset{1},subset{2},subset{3}),muXBU(subset{1},subset{2},subset{3}),muYBU(subset{1},subset{2},subset{3}),muZBU(subset{1},subset{2},subset{3})};
muL = {muXAL(subset{1},subset{2},subset{3}),muYAL(subset{1},subset{2},subset{3}),muZAL(subset{1},subset{2},subset{3}),muXBL(subset{1},subset{2},subset{3}),muYBL(subset{1},subset{2},subset{3}),muZBL(subset{1},subset{2},subset{3})};
macro_dipole_U = [macro_dipole_X_U, macro_dipole_Y_U, macro_dipole_Z_U];
macro_dipole_L = [macro_dipole_X_L, macro_dipole_Y_L, macro_dipole_Z_L];

elapsedTime = toc;
disp(strcat("Eigenstate calculation completed in ", num2str(elapsedTime), " seconds"));

end

