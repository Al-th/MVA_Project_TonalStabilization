function Aii = computeAdjustmentFrame(Ai,fi,fii,downsampleFactor)

downsampledFrame = permute(downsample(permute(downsample(fi,downsampleFactor),[2,1,3]),downsampleFactor),[2,1,3]);
downsampledFrame2 = permute(downsample(permute(downsample(fii,downsampleFactor),[2,1,3]),downsampleFactor),[2,1,3]);

[height,width,~] = size(downsampledFrame);
%%
%%Step 1 : bilateral filtering (spatial sigma = 10% of smaller dim, range
%%sigma = 10% of the values rnge
%disp('Bilateral filtering of images...');
spatialSigma = 0.1*min(width,height);
rangeSigma = 0.1;
filteredFrame = bfilter2(downsampledFrame,5,[spatialSigma,rangeSigma]);
filteredFrame2 = bfilter2(downsampledFrame2,5,[spatialSigma,rangeSigma]);
%disp('Done.');
%%
%%Step 2 : Set of correspondances
%%Ri/i+1 = {x s.t. |(Li(x)-mu(Li))?(Li+1(x)-mu(Li+1)| < 0.05}
%%L_{i} denote the luminance channel of the frame f_{i}
%%mu(L_{i}) represents the mean of the luminance channel L_{i}

%disp('Converting RBG to Lab ...');
labFrame = RGB2Lab(filteredFrame);
labFrame2 = RGB2Lab(filteredFrame2);
%disp('Done.');

%disp('Computing robust set of correspondances ...');
R = abs((labFrame(:,:,1) - mean(mean(labFrame(:,:,1)))) - (labFrame2(:,:,1) - mean(mean(labFrame2(:,:,1))))) < 0.5;

%disp('Done.');

%%
%%Step 3 : Initialize adjustment map
%%A_{i+1} = A_{i} + f_{i}(x) - f_{i+1}(x) if x in R, 0 otherwise


%disp('Computing initialisation of adjustment map');
A_init = zeros(height,width,3);
for i=1:height
    for j=1:width
        if(R(i,j))
            A_init(i,j,1) = Ai(i,j,1) + (labFrame(i,j,1) - labFrame2(i,j,1));
            A_init(i,j,2) = Ai(i,j,2) + (labFrame(i,j,2) - labFrame2(i,j,2));
            A_init(i,j,3) = Ai(i,j,3) + (labFrame(i,j,3) - labFrame2(i,j,3));
        end
    end
end
%disp('Done.');


%%
%disp('Fast interpolation using Nystrom method')
sigma1 = 5;
sigma2 = 5;
sigma3 = 5;

nbEigenValues = 10;

%disp('Computing SVD approximation for dimension 1...');
[S1,V1] = getNystromApproximation(labFrame(:,:,1),nbEigenValues,100,sigma1);
%disp('Computing SVD approximation for dimension 2...');
[S2,V2] = getNystromApproximation(labFrame(:,:,2),nbEigenValues,100,sigma2);
%disp('Computing SVD approximation for dimension 3...');
[S3,V3] = getNystromApproximation(labFrame(:,:,3),nbEigenValues,100,sigma3);
%disp('Done.');

V1 = V1';
V2 = V2';
V3 = V3';

%%

clear ChiACol AinitCol AfullCol Proj;

ChiACol(:,1) = reshape(R,width*height,1);
AinitCol(:,1) = reshape(A_init(:,:,1),width*height,1);

ChiACol(:,2) = reshape(R,width*height,1);
AinitCol(:,2) = reshape(A_init(:,:,2),width*height,1);

ChiACol(:,3) = reshape(R,width*height,1);
AinitCol(:,3) = reshape(A_init(:,:,3),width*height,1);

AfullCol(:,1) = (V1'*(S1(1:nbEigenValues,1:nbEigenValues)*V1*AinitCol(:,1))) ./ (V1'*(S1(1:nbEigenValues,1:nbEigenValues)*V1*ChiACol(:,1)));
AfullCol(:,2) = (V1'*(S1(1:nbEigenValues,1:nbEigenValues)*V1*AinitCol(:,2))) ./ (V1'*(S1(1:nbEigenValues,1:nbEigenValues)*V1*ChiACol(:,2)));
AfullCol(:,3) = (V1'*(S1(1:nbEigenValues,1:nbEigenValues)*V1*AinitCol(:,3))) ./ (V1'*(S1(1:nbEigenValues,1:nbEigenValues)*V1*ChiACol(:,3)));

Acompl(:,:,1) = reshape(AfullCol(:,1),height,width);
Acompl(:,:,2) = reshape(AfullCol(:,2),height,width);
Acompl(:,:,3) = reshape(AfullCol(:,3),height,width);

Aii = zeros(height,width,3);
for d = 1:3
    for i = 1:height
        for j = 1:width
            if(R(i,j))
                Aii(i,j,d) = A_init(i,j,d); 
            else
                Aii(i,j,d) = Acompl(i,j,d);
            end
        end
    end
end

%disp('Done computing adjustment map.');

end