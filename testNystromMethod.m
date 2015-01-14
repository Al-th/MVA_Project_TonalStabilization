%% Initialisation and ground truth
%Create matrice Image
%Create matrice W (pixel distance)
%Compute true USV for later comparison
clear all;

N = 10;

I = randn(N,N);

W = zeros(N*N,N*N);

sigma = 2;

[p,q] = meshgrid(1:N,1:N);
correspondances = [p(:),q(:)];
Icorr = I(sub2ind(size(I),correspondances(:,1),correspondances(:,2)));

for i = 1:N*N
    for j = 1:N*N
        W(i,j) = exp(-abs(Icorr(i)-Icorr(j))^2/sigma^2);
    end
end
imagesc(W);
pause(0.1);
tic
[U,S,V] = svd(W);
toc
%%
% Store all the possible correspondances
% Draw nbAff random correspondances
% Store the correspondances in randPixels
% Store the indexes
nbAff = 3;
tic
[U2,S2,V2] = getNystromApproximation(I,5,100,sigma);
toc
U2 = real(U2); 
S2 = real(S2);
V2 = real(V2);
% 
% hold on
% plot(V(:,2));   
% plot(-V2(2,:),'r--');


subplot(3,1,1)
imagesc(U*S*V');
subplot(3,1,2);
imagesc(U2*S2*V2);
subplot(3,1,3);
imagesc(real(U2*S2*V2) - U*S*V');