function [S,V] = getNystromApproximation(I,nbEigenvectors,nbSamples,sigma)
%1. Get samples
%2. Get correspondances between sampled and unsampled points
%3. Compute the approximation


%% 1
%disp('Computing affinity pairs matrix');
nbAff = nbSamples;
[M,N] = size(I);

if(nbSamples > M*N)
    error('Too much samples asked');
    return;
end

if(nbEigenvectors > nbSamples)
   error('Too much eigenvectors w.r.t number of samples');
   return;
end
A = zeros(nbAff, nbAff);

[p,q] = meshgrid(1:M,1:N);

correspondances = [p(:) q(:)];

sampleIndexes = [];
while((size(sampleIndexes,1)~=nbAff) || (size(sampleIndexes,1) ~= size(union(sampleIndexes,sampleIndexes),1)))
    missingIndexes = nbAff - size(sampleIndexes,1);
    sampleIndexes = [sampleIndexes; randi(M*N,missingIndexes,1)];
    sampleIndexes = union(sampleIndexes,sampleIndexes);
end
samplePixels = correspondances(sampleIndexes,:);
Isampled = I(sub2ind(size(I),samplePixels(:,1),samplePixels(:,2)));

unsampledIndexes = setdiff([1:M*N],sampleIndexes);
unsampledPixels = correspondances(unsampledIndexes,:);
Iunsampled = I(sub2ind(size(I),unsampledPixels(:,1),unsampledPixels(:,2)));

for i = 1:nbAff
   for j = 1:nbAff
       A(i,j) = exp(-abs(Isampled(i)-Isampled(j))^2/sigma^2);
   end
end
%disp('Done.');
%% 2
%disp('Computing affinity matrix with unsampled points');
B = zeros(nbAff,(M*N)-nbAff);

for i = 1:nbAff
    for j = 1:(M*N)-nbAff
        B(i,j) = exp(-abs(Isampled(i)-Iunsampled(j))^2/sigma^2);
    end
end
%disp('Done.');
%% 3
%disp('Computing eigenvectors');

% disp('Normalizing A and B for Laplacian...');
% B_T = B';
% d1 = sum(A, 2) + sum(B, 2);
% d2 = sum(B_T, 2) + B_T*(pinv(A)*sum(B, 2));
% dhat = sqrt(1./[d1; d2]);
% A = A .* (dhat(1:nbSamples)*dhat(1:nbSamples)');
% m = M*N - nbSamples;
% B1 = dhat(1:nbSamples)*dhat(nbSamples+(1:m))';
% B = B .* B1;


%disp('Orthogalizing and eigendecomposition...');
Asi = sqrtm(pinv(A));
B_T = B';
BBT = B*B_T;
W = single(zeros(size(A, 1)+size(B_T, 1), size(A, 2)));
W(1:size(A, 1), :) = A;
W(size(A, 1)+1:size(W, 1), :) = B_T;

% Calculate R = A + A^-1/2*B*B'*A^-1/2
R = A + Asi*BBT*Asi;
R = (R + R')/2; % Make sure R is symmetric, sometimes R can be non-symmetric because of numerical inaccuracy
[U L] = eig(R);
[val ind] = sort(diag(L), 'descend');
U = U(:, ind); % in decreasing order
L = L(ind, ind); % in decreasing order
W = W*Asi;
V = W*U(:, 1:nbEigenvectors)*pinv(sqrt(L(1:nbEigenvectors, 1:nbEigenvectors)));
S = L;



end