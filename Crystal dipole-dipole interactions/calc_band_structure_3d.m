function [KX,KY,KZ,Ubranch,Lbranch,JAA_k,JBB_k,JAB_k] = calc_band_structure_3d(XA,YA,ZA,dA,Ax,Ay,Az,JAA,JBB,JAB,epsA,epsB)

tic

JAA_clipped = fftshift(JAA(1:end-1,1:end-1,1:end-1));
JBB_clipped = fftshift(JBB(1:end-1,1:end-1,1:end-1));
JAB_clipped = fftshift(JAB(1:end-1,1:end-1,1:end-1));

% Fourier transform interaction matrices
JAA_k = real(fftn(JAA_clipped));
JBB_k = real(fftn(JBB_clipped));
JAB_k = fftn(JAB_clipped);

% calculate K vectors (-pi..pi)
KX = linspace(-pi/Ax,pi/Ax,numel(JAA_k(1,:,1))+1);
KY = linspace(-pi/Ay,pi/Ay,numel(JAA_k(:,1,1))+1);
KZ = linspace(-pi/Az,pi/Az,numel(JAA_k(1,1,:))+1);
KX = KX(1:end-1);
KY = KY(1:end-1);
KZ = KZ(1:end-1);

% shift interaction matrices so that the k=0 point is in the center
JAA_k = fftshift(JAA_k);
JBB_k = fftshift(JBB_k);
JAB_k = fftshift(JAB_k);

Ubranch = 0.5*(epsA+epsB+JAA_k+JBB_k)+0.5*sqrt((epsA-epsB+JAA_k-JBB_k).^2+4*abs(JAB_k).^2);
Lbranch = 0.5*(epsA+epsB+JAA_k+JBB_k)-0.5*sqrt((epsA-epsB+JAA_k-JBB_k).^2+4*abs(JAB_k).^2);

elapsedTime = toc;
disp(strcat("Band structure calculation completed in ", num2str(elapsedTime), " seconds"));

end

