function [KX,KY,Ubranch,Lbranch,JAA_k,JBB_k,JAB_k] = calc_band_structure_2d(Ax,Ay,JAA,JBB,JAB,epsA,epsB)

JAA_clipped = fftshift(JAA(1:end-1,1:end-1));
JBB_clipped = fftshift(JBB(1:end-1,1:end-1));
JAB_clipped = fftshift(JAB(1:end-1,1:end-1));

% Fourier transform interaction matrices
JAA_k = real(fftn(JAA_clipped));
JBB_k = real(fftn(JBB_clipped));
JAB_k = fftn(JAB_clipped);

% calculate K vectors (-pi..pi)
KX = linspace(-pi/Ax,pi/Ax,numel(JAA_k(1,:))+1);
KY = linspace(-pi/Ay,pi/Ay,numel(JAA_k(:,1))+1);
KX = KX(1:end-1);
KY = KY(1:end-1);

% shift interaction matrices so that the k=0 point is in the center
JAA_k = fftshift(JAA_k);
JBB_k = fftshift(JBB_k);
JAB_k = fftshift(JAB_k);

% % % % naive FFT implementation
% % % X = XA-dA(1);
% % % Y = YA-dA(2);
% % % KX = linspace(-pi/Ax,pi/Ax,numel(X(1,:,1)));
% % % KY = linspace(-pi/Ay,pi/Ay,numel(X(:,1,1)));
% % % JAA_k = zeros(size(JAA));
% % % JBB_k = zeros(size(JBB));
% % % JAB_k = zeros(size(JAB));
% % % for i=1:numel(KX)
% % %     for j=1:numel(KY)
% % %         JAA_k(j,i) = real(sum(sum(sum(exp(1i*(KX(i)*X+KY(j)*Y)).*JAA))));
% % %         JBB_k(j,i) = real(sum(sum(sum(exp(1i*(KX(i)*X+KY(j)*Y)).*JBB))));
% % %         JAB_k(j,i) = sum(sum(sum(exp(1i*(KX(i)*X+KY(j)*Y)).*JAB)));
% % %     end
% % % end

Ubranch = 0.5*(epsA+epsB+JAA_k+JBB_k)+0.5*sqrt((epsA-epsB+JAA_k-JBB_k).^2+4*abs(JAB_k).^2);
Lbranch = 0.5*(epsA+epsB+JAA_k+JBB_k)-0.5*sqrt((epsA-epsB+JAA_k-JBB_k).^2+4*abs(JAB_k).^2);

end

