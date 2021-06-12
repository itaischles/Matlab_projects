function [JAA,JBB,JAB,R] = calc_interaction_2d(XA_angstrom,YA_angstrom,XB_angstrom,YB_angstrom,dA_angstrom,dB_angstrom,muXA_debye,muYA_debye,muXB_debye,muYB_debye,dielconst)

eps0 = 8.8541878e-12;

muXA = muXA_debye*3.3356e-30;
muYA = muYA_debye*3.3356e-30;
muXB = muXB_debye*3.3356e-30;
muYB = muYB_debye*3.3356e-30;

XA = XA_angstrom*1e-10;
YA = YA_angstrom*1e-10;
XB = XB_angstrom*1e-10;
YB = YB_angstrom*1e-10;
dA = dA_angstrom*1e-10;
dB = dB_angstrom*1e-10;

R = sqrt((XA-dA(1)).^2+(YA-dA(2)).^2);
RAB = sqrt((XB-dA(1)).^2+(YB-dA(2)).^2);

% Dipole-dipole interaction
JAA = ((muXA.^2+muYA.^2).*R.^2 - 3*(muXA.*(XA-dA(1))+muYA.*(YA-dA(2))).^2)./R.^5 * 1/(4*pi*eps0*dielconst);
JAA(isnan(JAA))=0;
JBB = ((muXB.^2+muYB.^2).*R.^2 - 3*(muXB.*(XB-dB(1))+muYB.*(YB-dB(2))).^2)./R.^5 * 1/(4*pi*eps0*dielconst);
JBB(isnan(JBB))=0;
JAB = ((muXA.*muXB+muYA.*muYB).*RAB.^2-3*(muXA.*(XB-dA(1))+muYA.*(YB-dA(2))).*(muXB.*(XB-dA(1))+muYB.*(YB-dA(2))))./RAB.^5 * 1/(4*pi*eps0*dielconst);

% % % % nearest neighbor interaction
% % % NNdist_RAB = min(min(min(RAB)));
% % % NNdist_R = min(min(min(R(R~=0))));
% % % JAA = zeros(size(RAB));
% % % JBB = zeros(size(RAB));
% % % JAB = zeros(size(RAB));
% % % JAA(R==NNdist_R) = 1;
% % % JBB(R==NNdist_R) = 1;
% % % JAB(RAB==NNdist_RAB) = 1;

% % % % 1/R^n interaction
% % % JAA = 1./R.^6;
% % % JBB = 1./R.^6;
% % % JAB = 1./RAB.^6;
% % % JAA(isinf(JAA))=0;
% % % JBB(isinf(JBB))=0;
% % % JAB(isinf(JAB))=0;

% % % % exponential interaction
% % % R0 = 1;
% % % JAA = ones(size(R));
% % % JBB = ones(size(R));
% % % JAB = ones(size(RAB));
% % % JAA = JAA.*(1.0*exp(-R/R0));
% % % JBB = JBB.*(1.0*exp(-R/R0));
% % % JAB = JAB.*(1.0*exp(-RAB/R0));

JAA = JAA/1.6e-19;
JBB = JBB/1.6e-19;
JAB = JAB/1.6e-19;

end

