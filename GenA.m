close all
clear
clc

%% Construct the A matrix.
load('ParamUsed.mat');
load('VoxelSize.mat');
load('FixedFreqUsed.mat');
clear FreqInd

NxVoxel = length(x_v);
NyVoxel = length(y_v);
NzVoxel = length(z_v);

% Calculate all the distances needed.
VoxelCoord = combvec(x_v, y_v, z_v);
p_xVoxel = VoxelCoord(1, :);
p_yVoxel = VoxelCoord(2, :);
p_zVoxel = VoxelCoord(3, :);

p_tagx = zeros(TagNum, NxVoxel*NyVoxel*NzVoxel);
p_tagy = zeros(TagNum, NxVoxel*NyVoxel*NzVoxel);
p_tagz = zeros(TagNum, NxVoxel*NyVoxel*NzVoxel);
p_tagx(:, 1) = tagPosition(:, 1);
p_tagy(:, 1) = tagPosition(:, 2);
p_tagz(:, 1) = tagPosition(:, 3);
p_tagx = repmat(p_tagx(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_tagy = repmat(p_tagy(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_tagz = repmat(p_tagz(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);

p_recvx = zeros(RecvNum, NxVoxel*NyVoxel*NzVoxel);
p_recvy = zeros(RecvNum, NxVoxel*NyVoxel*NzVoxel);
p_recvz = zeros(RecvNum, NxVoxel*NyVoxel*NzVoxel);
p_recvx(:, 1) = rxPosition(:, 1);
p_recvy(:, 1) = rxPosition(:, 2);
p_recvz(:, 1) = rxPosition(:, 3);
p_recvx = repmat(p_recvx(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_recvy = repmat(p_recvy(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);
p_recvz = repmat(p_recvz(:, 1), [1, NxVoxel*NyVoxel*NzVoxel]);

DistTagVoxel = sqrt((p_tagx - repmat(p_xVoxel, [TagNum, 1])).^2 ...
    +(p_tagy - repmat(p_yVoxel, [TagNum, 1])).^2+(p_tagz - repmat(p_zVoxel, [TagNum, 1])).^2);

DistRecvVoxel = sqrt((p_recvx - repmat(p_xVoxel, [RecvNum, 1])).^2 ...
    +(p_recvy - repmat(p_yVoxel, [RecvNum, 1])).^2+(p_recvz - repmat(p_zVoxel, [RecvNum, 1])).^2);

DistTemp = repmat(DistTagVoxel, RecvNum, 1) + repelem(DistRecvVoxel, TagNum, 1);
DistUplink = repmat(DistTemp, length(freq), 1);

% Calculate the frequency constant.
FreqConst = (repelem(-1j*2*pi*freq/c, TagNum*RecvNum)).';
Freq_ratio = (repelem(freq/min(freq), TagNum*RecvNum).^2).';
% A = Freq_ratio.*exp(FreqConst.*DistUplink);
A = exp(FreqConst.*DistUplink);

save('A_mat.mat', 'A');

