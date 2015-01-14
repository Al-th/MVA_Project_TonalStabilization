clear;
clc;
movieTitle = 'data/entrance.AVI';
disp(['Starting algorithm, loading movie ' movieTitle]);
mov = VideoReader('data/entrance.AVI');
vidFrames = read(mov);
nbFrames = get(mov,'NumberOfFrames');
disp('Done.');
%%

%nbFrames = 2;

clear A;
clear frameCorrected;

disp('Loading specific frames');
A = zeros(120,160,3,nbFrames);


frame = double(vidFrames(:,:,:,1));
for i = 1:nbFrames
    fprintf('Computing adjustmentFrame %d%',100*(i/nbFrames));
    frame2 = double(vidFrames(:,:,:,i+1));

    frame = frame./max(frame(:));
    frame2 = frame2./max(frame2(:));
    
    A(:,:,:,i+1) = computeAdjustmentFrame(A(:,:,:,i),frame,frame2,4);
    
    frame = frame2;
end

%%Correcting movie with adjustmentMap

upsampledA= permute(upsample(permute(upsample(A,4),[2,1,3,4]),4),[2,1,3,4]);

for i = 1:nbFrames
    fprintf('Computing corrected frames %d%', 100*(i/nbFrames));
    frameCorrected(:,:,:,i) = Lab2RGB(RGB2Lab(double(vidFrames(:,:,:,i))/255.0) + upsampledA(:,:,:,i));
    
end
%%

