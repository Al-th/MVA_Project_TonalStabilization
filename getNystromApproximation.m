function [U,S,V] = getNystromApproximation(I,nbEigenvectors,nbSamples,sigma)
%1. Get samples
%2. Get correspondances between sampled and unsampled points
%3. Compute the approximation


%% 1
disp('Computing affinity pairs matrix');
nbAff = nbSamples;
[M,N] = size(I);

if(nbSamples > M*N)
    disp('Too much samples asked');
    return;
end

if(nbEigenvectors > nbSamples)
   disp('Too much eigenvectors w.r.t number of samples');
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
disp('Done.');
%% 2
disp('Computing affinity matrix with unsampled points');
B = zeros(nbAff,(M*N)-nbAff);

for i = 1:nbAff
    for j = 1:(M*N)-nbAff
        B(i,j) = exp(-abs(Isampled(i)-Iunsampled(j))^2/sigma^2);
    end
end
disp('Done.');
%% 3
disp('Computing eigenvectors');

[~,p] = chol(A);
if(p==0)
    disp('A is sym definite');
    
    [~,m] = size(B);
    [n,~] = size(A);
    d1 = sum([A;B'],1);
    d2 = sum(B,1) + sum(B',1)*pinv(A)*B;
    dhat = sqrt(1 ./ [d1 d2]);
%     A = A .* (dhat(1:n)*dhat(1:n)');
%     B = B .* (dhat(1:n)*dhat(n+(1:m))');
%     
    
    Asi = sqrtm(pinv(A));
    Q = A + Asi*(B*B')*Asi;
    [Us,Ss,Vs] = svd(Q);
    Vtemp = ([A;B']*(Asi*(Us*pinv(sqrtm(Ss)))));
    S = zeros(nbEigenvectors,nbEigenvectors);
    for i = 1:nbEigenvectors
       V(:,i) = Vtemp(:,i);
       S(i,i) = Ss(i);
    end
    U = V;

    disp('Done.');
else
    disp('A is not definite positive, using alternative method');
    
    Asi = sqrtm(pinv(A));
    Q = A + Asi*(B*B')*Asi;
    [Us,Ss,Vs] = svd(Q);
    Ub = [Us' pinv(Ss)*Us'*B];
    Z = Ub'*sqrtm(Ss);
    [Uz,Sz,Vz] = svd(Z'*Z);
    Szsi = sqrtm(pinv(Sz));
    Vtemp = Z*Uz*Szsi;
    S = zeros(nbEigenvectors,nbEigenvectors);
    for i = 1:nbEigenvectors
       V(i,:) = Vtemp(:,i)';
       S(i,i) = Sz(i);
    end
    U = V';
end
end