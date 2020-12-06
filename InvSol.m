close all
clear
clc

%% Inverse solution.
load('ParamUsed.mat');
load('VoxelSize.mat');
load('ReconsData.mat');
load('A_mat.mat');

% % Matched filtering.
% load('A_mat.mat');
% tic;
% imgComplex = A'*b;
% tSol = toc;
% NormCal = norm(A*imgComplex-b);


% % Moore-Penrose pseudo-inverse.
% load('svdComp.mat');
% imgComplex = zeros(size(A, 2), 1);
% sVec = diag(S);
% K = length(sVec);
% for i = 1:K
%     uVec = U(:, i);
%     vVec = V(:, i);
%     imgComplex = imgComplex + (1/sVec(i)).*(uVec'*b).*vVec;
% end

% % Alternative implementation.
% load('A_pinv_mat.mat');
% imgComplex = A_pinv*b;
% NormCal = norm(A*imgComplex-b);


% % Truncated Moore-Penrose pseudo-inverse. (tol = 200 good, K = 162)
% load('svdComp.mat');
% imgComplex = zeros(size(A, 2), 1);
% sVec = diag(S);
% K = 162;
% for i = 1:K
%     uVec = U(:, i);
%     vVec = V(:, i);
%     imgComplex = imgComplex + (1/sVec(i)).*(uVec'*b).*vVec;
% end

% % Alternative implementation.
% load('A_pinv_trunc_mat.mat');
% imgComplex = A_pinv*b;
% NormCal = norm(A*imgComplex-b);


% % Moore-Penrose pseudo-inverse with zeroth-order Tikhonov
% % regularization.
% load('svdComp.mat');
% imgComplex = zeros(size(A, 2), 1);
% sVec = diag(S);
% K = length(sVec);
% lambda = 1*10^(-6.5);
% for i = 1:K
%     uVec = U(:, i);
%     vVec = V(:, i);
%     imgComplex = imgComplex + (sVec(i)/(sVec(i)^2 + lambda^2))*(uVec'*b)*vVec;
% end

% % Alternative implementation.
% load('A_ZeroTik_mat.mat');
% b_new = A.'*b;
% imgComplex = A_ZeroTik*b_new;
% NormCal = norm(A*imgComplex-b);


% Truncated Moore-Penrose pseudo-inverse with zeroth-order Tikhonov
% regularization. (alpha around 10^(2.5) is the critical point)
load('svdComp.mat');
imgComplex = zeros(size(A, 2), 1);
sVec = diag(S);
K = 300;
lambda = 1*10^(4);
for i = 1:K
    uVec = U(:, i);
    vVec = V(:, i);
    imgComplex = imgComplex + (sVec(i)/(sVec(i)^2 + lambda^2))*(uVec'*b)*vVec;
end

% % Alternative implementation.
% load('A_trunc_ZeroTik_mat.mat');
% b_new = A.'*b;
% imgComplex = A_trunc_ZeroTik*b_new;
% NormCal = norm(A*imgComplex-b);

%% Post-Processing.
% Apply threshold for reconstructed reflectivity intensity.
imgComplexAbs = abs(imgComplex).^PowerOrder;
% imgComplexAbs = abs(imgComplex);
ReconsDistMax = max(imgComplexAbs);
ReconsDistMin = min(imgComplexAbs);
ReconsDistNorm = (imgComplexAbs - ReconsDistMin)/(ReconsDistMax - ReconsDistMin);

% Threshold = ReconsDistMax*thresCoeff;
% imgComplexAbs = reshape(imgComplexAbs, [length(x_v), length(y_v), length(z_v)]);
% imgComplexAbs(imgComplexAbs < Threshold) = 0;
% imgComplexAbs(imgComplexAbs >= Threshold) = 1;
% 
% % Clustering.
% roomSize = [x_v(1), x_v(length(x_v)); y_v(1), y_v(length(y_v)); z_v(1), z_v(length(z_v))];
% voxelSize = [x_v(2)-x_v(1); y_v(2)-y_v(1); z_v(2)-z_v(1)];
% [cDistribution, clusters,centroidDist] = i4block_components(imgComplexAbs, roomSize, voxelSize);

Threshold = thresCoeff;
ReconsDistNorm = reshape(ReconsDistNorm, [length(x_v), length(y_v), length(z_v)]);
ReconsDistNorm(ReconsDistNorm < Threshold) = 0;
ReconsDistNorm(ReconsDistNorm >= Threshold) = 1;

% Clustering.
roomSize = [x_v(1), x_v(length(x_v)); y_v(1), y_v(length(y_v)); z_v(1), z_v(length(z_v))];
voxelSize = [x_v(2)-x_v(1); y_v(2)-y_v(1); z_v(2)-z_v(1)];
[cDistribution, clusters,centroidDist] = i4block_components(ReconsDistNorm, roomSize, voxelSize);

fprintf('Initial cluster number = %d\n',length(clusters.centroid));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.7; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.7;
opts.elemNumTh = 0.6; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
% clusterOut = clusterProcess(clusters,opts);
clusterOut = clusterProcess_best(clusters,opts);

centroid = clusterOut.centroid;
elemNum = clusterOut.elemNum;
for i = 1:size(centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
end


%% Plotting.
% Convert the relative/normalized brightness of the reconstructed
% reflectivity vector into a 3D matrix.
% ReconsDistNorm = reshape(ReconsDistNorm, [length(x_v), length(y_v), length(z_v)]);

% % Plot the relative/normalized brightness of the reconstructed
% % reflectivity, in a 3D grid.
% [X_V, Y_V, Z_V] = meshgrid(x_v, y_v, z_v);
% imgComplexAbs = permute(imgComplexAbs, [2 1 3]);
% h = slice(X_V, Y_V, Z_V, imgComplexAbs,x_v,y_v,z_v);
% % imgComplexAbs(imgComplexAbs == 0) = NaN;
% % h = slice(X_V, Y_V, Z_V, imgComplexAbs,x_v,y_v,z_v);
% xlabel('x (m)','FontSize',14);
% ylabel('y (m)','FontSize',14);
% zlabel('z (m)','FontSize',14);
% xlim([x_v(1), x_v(length(x_v))]);
% ylim([y_v(1), y_v(length(y_v))]);
% zlim([z_v(1), z_v(length(z_v))]);

% Plot the relative/normalized brightness of the reconstructed
% reflectivity, in a 3D grid.
[X_V, Y_V, Z_V] = meshgrid(x_v, y_v, z_v);
ReconsDistNorm = permute(ReconsDistNorm, [2 1 3]);
h = slice(X_V, Y_V, Z_V, ReconsDistNorm,x_v,y_v,z_v);
% imgComplexAbs(imgComplexAbs == 0) = NaN;
% h = slice(X_V, Y_V, Z_V, imgComplexAbs,x_v,y_v,z_v);
xlabel('x (m)','FontSize',14);
ylabel('y (m)','FontSize',14);
zlabel('z (m)','FontSize',14);
xlim([x_v(1), x_v(length(x_v))]);
ylim([y_v(1), y_v(length(y_v))]);
zlim([z_v(1), z_v(length(z_v))]);

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color');
a = alphamap('rampup',256);
% imgThresh = 180;
a(1:imgThresh) = 0;
alphamap(a);
set(gca, 'fontweight', 'bold');
% figurepalette('show')
set(gca, 'fontname', 'times', 'fontweight', 'bold', 'fontsize', 16);


figure;
ReconsDistNorm2D = sum(ReconsDistNorm, 3);
hh = pcolor(x_v, y_v, ReconsDistNorm2D);
set(hh, 'EdgeColor', 'none');
hh.FaceColor = 'interp';

myColorMap = jet(256);
myColorMap(1,:) = 1;
line = 250;
% myColorMap(2:line,1) = 0.9290;
% myColorMap(2:line,2) = 0.6940;
% myColorMap(2:line,3) = 0.1250;
myColorMap(2:line,1) = 0.9;
myColorMap(2:line,2) = 0.65;
myColorMap(2:line,3) = 0.1;
colormap(myColorMap);
% colorbar
% mesh(x_v, y_v, ReconsDistNorm2D);

set(gca, 'FontSize', 12);
xt = 0:0.5:4;
yt = 0:0.5:4;
set(gca, 'xtick', xt, 'xticklabel', xt);
set(gca, 'ytick', yt, 'yticklabel', yt);
xlim([-0.4, 4.2]);
ylim([-0.4, 4.2]);
zlim([-0.2, 2.4]);
xlabel('x / m');
ylabel('y / m');
zlabel('z / m');
set(gca, 'fontname', 'times', 'fontweight', 'bold', 'fontsize', 16);


g = figure;
for i = 1:size(clusterOut.centroid, 1)
    scatter(clusterOut.centroid(i, 1), clusterOut.centroid(i, 2), 100, 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'Marker', 'o');
    hold on
end
set(gca, 'FontSize', 12);
xt = 0:0.5:4;
yt = 0:0.5:4;
set(gca, 'xtick', xt, 'xticklabel', xt);
set(gca, 'ytick', yt, 'yticklabel', yt);
xlim([-0.4, 4.2]);
ylim([-0.4, 4.2]);
xlabel('x / m');
ylabel('y / m');
set(gca, 'fontname', 'times', 'fontweight', 'bold', 'fontsize', 16);

switch postureloc
    case 'x2y2_'
        scatter(0.6, 0.6, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x2y4_'
        scatter(0.6, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x2y6_', 'x2y6_2_'}
        scatter(0.6, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x2y8_'
        scatter(0.6, 2.4, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x2y10_'
        scatter(0.6, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x4y2_'
        scatter(1.2, 0.6, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x4y4_', 'x4y4_2_'}
        scatter(1.2, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x4y6_', 'x4y6_2_'}
        scatter(1.2, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x4y8_', 'x4y8_2_'}
        scatter(1.2, 2.4, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x4y10_'
        scatter(1.2, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y2_', 'x6y2_2_'}
        scatter(1.8, 0.6, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y4_', 'x6y4_2_'}
        scatter(1.8, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y6_', 'x6y6_2_'}
        scatter(1.8, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y8_', 'x6y8_2_'}
        scatter(1.8, 2.4, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y10_', 'x6y10_2_'}
        scatter(1.8, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x8y2_'
        scatter(2.4, 0.6, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x8y4_', 'x8y4_2_'}
        scatter(2.4, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x8y6_', 'x8y6_2_'}
        scatter(2.4, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x8y8_', 'x8y8_2_'}
        scatter(2.4, 2.4, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x8y10_'
        scatter(2.4, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x10y2_'
        scatter(3.0, 0.6, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x10y4_'
        scatter(3.0, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x10y6_', 'x10y6_2_'}
        scatter(3.0, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x10y8_'
        scatter(3.0, 2.4, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case 'x10y10_'
        scatter(3.0, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x4y4_x8y4_', 'x4y4_x8y4_2_'}
        scatter(1.2, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd'); hold on
        scatter(2.4, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x4y4_x8y8_', 'x4y4_x8y8_2_'}
        scatter(1.2, 1.2, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd'); hold on
        scatter(2.4, 2.4, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y2_x2y6_', 'x6y2_x2y6_2_'}
        scatter(1.8, 0.6, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd'); hold on
        scatter(0.6, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y2_x6y10_', 'x6y2_x6y10_2_'}
        scatter(1.8, 0.6, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd'); hold on
        scatter(1.8, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y6_x6y10_', 'x6y6_x6y10_2_'}
        scatter(1.8, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd'); hold on
        scatter(1.8, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
    case {'x6y10_x6y6_', 'x6y10_x6y6_2_'}
        scatter(1.8, 3.0, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd'); hold on
        scatter(1.8, 1.8, 100, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'Marker', 'd');
end

saveas(g, [personname, postureloc, '.emf']);
