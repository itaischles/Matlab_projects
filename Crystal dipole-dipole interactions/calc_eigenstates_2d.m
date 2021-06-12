function [muU,muL,macro_dipole_U,macro_dipole_L] = calc_eigenstates_2d(KX,KY,Ax,Ay,XA,YA,dA,dB,JAA_k,JBB_k,JAB_k,JAB_cutoff,Ubranch,Lbranch,epsA,epsB,muXA,muYA,muXB,muYB,eigenstate_to_plot)


delta_d = dB-dA;

if ~isempty(eigenstate_to_plot) && numel(eigenstate_to_plot)==2
    
    [~,i] = min(abs(KX*Ax-eigenstate_to_plot(1)));
    [~,j] = min(abs(KY*Ay-eigenstate_to_plot(2)));
    
    if abs(JAB_k(j,i))>JAB_cutoff
        cA_k_U = 1;
        cB_k_U = (Ubranch(j,i)-epsA-JAA_k(j,i))/(JAB_k(j,i)*exp(1i*(KX(i)*delta_d(1)+KY(j)*delta_d(2))));
        cA_k_L = 1;
        cB_k_L = (Lbranch(j,i)-epsA-JAA_k(j,i))/(JAB_k(j,i)*exp(1i*(KX(i)*delta_d(1)+KY(j)*delta_d(2))));
    else
        if JAA_k(j,i)+epsA>JBB_k(j,i)+epsB
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
    
    phase = exp(1i*(KX(i)*(XA-dA(1))+KY(j)*(YA-dA(2))));
    
    eigenstates_U_A = phase.*cA_k_U*exp(1i*(KX(i)*dA(1)+KY(j)*dA(2)));
    eigenstates_L_A = phase.*cA_k_L*exp(1i*(KX(i)*dA(1)+KY(j)*dA(2)));
    eigenstates_U_B = phase.*cB_k_U*exp(1i*(KX(i)*dA(1)+KY(j)*dB(2)));
    eigenstates_L_B = phase.*cB_k_L*exp(1i*(KX(i)*dA(1)+KY(j)*dB(2)));
    
    muXAU = real(eigenstates_U_A).*muXA;
    muYAU = real(eigenstates_U_A).*muYA;
    muXBU = real(eigenstates_U_B).*muXB;
    muYBU = real(eigenstates_U_B).*muYB;
    
    muXAL = real(eigenstates_L_A).*muXA;
    muYAL = real(eigenstates_L_A).*muYA;
    muXBL = real(eigenstates_L_B).*muXB;
    muYBL = real(eigenstates_L_B).*muYB;
    
    macro_dipole_X_U = sum(sum(sum(muXAU+muXBU)))/numel(XA);
    macro_dipole_Y_U = sum(sum(sum(muYAU+muYBU)))/numel(XA);
    macro_dipole_X_L = sum(sum(sum(muXAL+muXBL)))/numel(XA);
    macro_dipole_Y_L = sum(sum(sum(muYAL+muYBL)))/numel(XA);
    
end

muU = {muXAU,muYAU,muXBU,muYBU};
muL = {muXAL,muYAL,muXBL,muYBL};
macro_dipole_U = [macro_dipole_X_U, macro_dipole_Y_U];
macro_dipole_L = [macro_dipole_X_L, macro_dipole_Y_L];


end

