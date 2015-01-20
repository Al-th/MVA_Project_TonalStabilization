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
subplot(2,7,1)
imagesc(abs(U*S*V'));
freezeColors();
subplot(2,7,8);
imagesc(abs(U*S*V'-U*S*V'));

for i = 1:5
    clear S2 V2;
    nbAff = 100;
    nbEigen = ((i-1)*2)+1;
    tic
    [S2,V2] = getNystromApproximation(I,nbEigen,nbAff,sigma);
    toc

    S2 = real(S2);
    V2 = real(V2);

    subplot(2,7,i+1);
    imagesc(abs(V2*S2(1:nbEigen,1:nbEigen)*V2'));
    freezeColors();
    subplot(2,7,i+8);
    imagesc(abs(V2*S2(1:nbEigen,1:nbEigen)*V2' - U*S*V'));
    
    difference = single(V2*S2(1:nbEigen,1:nbEigen)*V2') - single(U*S*V');
    squaredError = difference .^ 2;
    meanSquaredError = sum(squaredError(:)) / numel(V2*S2(1:nbEigen,1:nbEigen)*V2');
    rmsError(i) = sqrt(meanSquaredError);
end