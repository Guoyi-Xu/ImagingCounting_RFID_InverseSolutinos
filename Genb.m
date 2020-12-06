close all
clear
clc

% Parameters.
[tagPosition, rxPosition, ~] = tag_antenna_positions3D_func();
TagNum = size(tagPosition, 1);
RecvNum = size(rxPosition, 1);
FreqNumTotal = 50;
FreqNumUse = 6;
c = physconst('LightSpeed');
PowerOrder = 3;
thresCoeff = 0.8;
imgThresh = 0;
savefigure = 0;
save('ParamUsed.mat', 'TagNum', 'RecvNum', 'tagPosition', 'rxPosition', 'c', 'PowerOrder', 'thresCoeff', 'imgThresh');

% Choose the fixed frequencies to use.
freq = 1000*(902750:500:927250);
UseSavedFreq = 1;
if UseSavedFreq == 1
    load('FixedFreqUsed.mat');
else
    delete('FixedFreqUsed.mat');
    FreqInd = sort(randsample(50, FreqNumUse));
    freq = freq(FreqInd);
    save('FixedFreqUsed.mat', 'FreqInd', 'freq');
end


% The voxel coordinates.
x_v=0:0.12:3.6;
y_v=0:0.12:3.6;
z_v=0:0.3:2.4;
save('VoxelSize.mat', 'x_v', 'y_v', 'z_v');

% Threshold for standard deviation of phase.
phStdTh = 25; % In degrees. A very high value means no rejection.
% fprintf('Using phase standard deviation threshold %d\n', phStdTh);
rssiDiffTh = 10^(-0.3);
RJT = 0;
rssiDiffThAllFreq = 10^(-0.3);
RJTAllFreq = 1;
opts.calibType = 5;
OneReadRej = 0;


% Copy the raw data files to the current directory, and extract the data
% files.
year = '2020';
date = '0316';
personname = 'phantom_';
postureloc = 'x10y6_';

file_wo = ['D:\RFImagingExptDataProcess\test_', year, date, '\data_wo_', personname, postureloc, year, date, '.mat'];
file_w = ['D:\RFImagingExptDataProcess\test_', year, date, '\data_w_', personname, postureloc, year, date, '.mat'];
XLSfile_wo = ['D:\RFImagingExptDataProcess\test_', year, date, '\data_wo_', personname, postureloc, year, date, '.xlsx'];
XLSfile_w = ['D:\RFImagingExptDataProcess\test_', year, date, '\data_w_', personname, postureloc, year, date, '.xlsx'];
dir_cur = ['D:\RFImagingExptDataProcess\test_', year, date, '_Process'];

HopTable_wo = readmatrix(XLSfile_wo, 'Sheet', 3);
HopTable_wo = HopTable_wo(2:length(HopTable_wo), :);
HopTable_w = readmatrix(XLSfile_w, 'Sheet', 3);
HopTable_w = HopTable_w(2:length(HopTable_w), :);
[~, I_wo] = sort(HopTable_wo(:, 2));
[~, I_w] = sort(HopTable_w(:, 2));

%% Data Collection. 
% Load the raw data.
load(file_wo);

% Generate the data matrices.
PhaseDataCollect1 = cell(TagNum, RecvNum, FreqNumTotal);
PhaseCollectStd1 = zeros(TagNum, RecvNum, FreqNumTotal);
PhaseCollectMean1 = zeros(TagNum, RecvNum, FreqNumTotal);
RSSIDataCollect1 = cell(TagNum, RecvNum, FreqNumTotal);
RSSICollectStd1 = zeros(TagNum, RecvNum, FreqNumTotal);
RSSICollectMean1 = zeros(TagNum, RecvNum, FreqNumTotal);
CompNumDataCollect1 = cell(TagNum, RecvNum, FreqNumTotal);
CompNumCollectMean1 = zeros(TagNum, RecvNum, FreqNumTotal);
CountOneRead1 = 0;
for j = 1:RecvNum
    phase_t = phasedeglist(antennalist == j);
    rssi_t = rssiimpinjlist_d(antennalist == j);
    cha_t = chindlist(antennalist == j);
    tag_t = tagindexlist(antennalist == j);
    for k = 1:FreqNumTotal
        phase_tt = phase_t(cha_t == k);
        rssi_tt = rssi_t(cha_t == k);
        tag_tt = tag_t(cha_t == k);
        for i = 1:TagNum
            phase_ttt = phase_tt(tag_tt == i);
            phase_ttt = rad2deg(unwrap(deg2rad(phase_ttt)));
            rssi_ttt = rssi_tt(tag_tt == i);
            if ~isempty(phase_ttt)
                PhaseDataCollect1{i, j, k} = phase_ttt;
                PhaseCollectStd1(i, j, k) = std(phase_ttt);
                RSSIDataCollect1{i, j, k} = rssi_ttt;
                RSSICollectStd1(i, j, k) = std(rssi_ttt);
                CompNumDataCollect1{i, j, k} = rssi_ttt.*cos(deg2rad(phase_ttt)) + sqrt(-1).*rssi_ttt.*sin(deg2rad(phase_ttt));
                
                % Option: reject channels with only one read.
                if OneReadRej == 1
                    if PhaseCollectStd1(i, j, k) == 0
                        phase_ttt = NaN;
                        rssi_ttt = NaN;
                        CompNumDataCollect1{i, j, k} = [];
                        PhaseCollectMean1(i, j, k) = NaN;
                        RSSICollectMean1(i, j, k) = NaN;
                        CompNumCollectMean1(i, j, k) = NaN;                        
                        CountOneRead1 = CountOneRead1 + 1;
                    end
                end
                
                if PhaseCollectStd1(i, j, k) > phStdTh
                    % First completely reject the channel using phase
                    % standard deviation threshold.
                    PhaseCollectMean1(i, j, k) = NaN;
                    RSSICollectMean1(i, j, k) = NaN;
                    CompNumCollectMean1(i, j, k) = NaN;
                else
                    PhaseCollectMean1(i, j, k) = mean(phase_ttt);
                    RSSICollectMean1(i, j, k) = mean(rssi_ttt);
                    CompNumCollectMean1(i, j, k) = mean(CompNumDataCollect1{i, j, k});
                end
            else
                % Simply there's no data received for this channel.
                PhaseDataCollect1{i, j, k} = [];
                PhaseCollectStd1(i, j, k) = NaN;
                PhaseCollectMean1(i, j, k) = NaN;
                RSSIDataCollect1{i, j, k} = [];
                RSSICollectStd1(i, j, k) = NaN;
                RSSICollectMean1(i, j, k) = NaN;
                CompNumDataCollect1{i, j, k} = NaN;
                CompNumCollectMean1(i, j, k) = NaN;
            end                
        end
    end
end
fprintf('Number of channels with only one read: %3.0f \n', CountOneRead1);
fprintf('Percentage channels with only one read: %3.2f%% \n', 100*CountOneRead1/length(PhaseCollectStd1(:)));
clear chindexlist tagindexlist antennalist rssiimpinjlist rssiimpinjlist_d phasedeglist msgfreqlist

% Adjust the order in the frequency dimension.
PhaseDataCollect1 = PhaseDataCollect1(:, :, I_wo);
PhaseCollectStd1 = PhaseCollectStd1(:, :, I_wo);
PhaseCollectMean1 = PhaseCollectMean1(:, :, I_wo);
RSSIDataCollect1 = RSSIDataCollect1(:, :, I_wo);
RSSICollectStd1 = RSSICollectStd1(:, :, I_wo);
RSSICollectMean1 = RSSICollectMean1(:, :, I_wo);
CompNumDataCollect1 = CompNumDataCollect1(:, :, I_wo);
CompNumCollectMean1 = CompNumCollectMean1(:, :, I_wo);

% Pick out the data associated with the selected fixed frequencies.
PhaseDataCollect1 = PhaseDataCollect1(:, :, FreqInd);
PhaseCollectStd1 = PhaseCollectStd1(:, :, FreqInd);
PhaseCollectMean1 = PhaseCollectMean1(:, :, FreqInd);
RSSIDataCollect1 = RSSIDataCollect1(:, :, FreqInd);
RSSICollectStd1 = RSSICollectStd1(:, :, FreqInd);
RSSICollectMean1 = RSSICollectMean1(:, :, FreqInd);
CompNumDataCollect1 = CompNumDataCollect1(:, :, FreqInd);
CompNumCollectMean1 = CompNumCollectMean1(:, :, FreqInd);


% Load the raw data.
load(file_w);

% Generate the data matrices.
PhaseDataCollect2 = cell(TagNum, RecvNum, FreqNumTotal);
PhaseCollectStd2 = zeros(TagNum, RecvNum, FreqNumTotal);
PhaseCollectMean2 = zeros(TagNum, RecvNum, FreqNumTotal);
RSSIDataCollect2 = cell(TagNum, RecvNum, FreqNumTotal);
RSSICollectStd2 = zeros(TagNum, RecvNum, FreqNumTotal);
RSSICollectMean2 = zeros(TagNum, RecvNum, FreqNumTotal);
CompNumDataCollect2 = cell(TagNum, RecvNum, FreqNumTotal);
CompNumCollectMean2 = zeros(TagNum, RecvNum, FreqNumTotal);
CountOneRead2 = 0;
for j = 1:RecvNum
    phase_t = phasedeglist(antennalist == j);
    rssi_t = rssiimpinjlist_d(antennalist == j);
    cha_t = chindlist(antennalist == j);
    tag_t = tagindexlist(antennalist == j);
    for k = 1:FreqNumTotal
        phase_tt = phase_t(cha_t == k);
        rssi_tt = rssi_t(cha_t == k);
        tag_tt = tag_t(cha_t == k);
        for i = 1:TagNum
            phase_ttt = phase_tt(tag_tt == i);
            phase_ttt = rad2deg(unwrap(deg2rad(phase_ttt)));
            rssi_ttt = rssi_tt(tag_tt == i);
            if ~isempty(phase_ttt)
                PhaseDataCollect2{i, j, k} = phase_ttt;
                PhaseCollectStd2(i, j, k) = std(phase_ttt);
                RSSIDataCollect2{i, j, k} = rssi_ttt;
                RSSICollectStd2(i, j, k) = std(rssi_ttt);
                CompNumDataCollect2{i, j, k} = rssi_ttt.*cos(deg2rad(phase_ttt)) + sqrt(-1).*rssi_ttt.*sin(deg2rad(phase_ttt));
                
                % Option: reject channels with only one read.
                if OneReadRej == 1
                    if PhaseCollectStd2(i, j, k) == 0
                        phase_ttt = NaN;
                        rssi_ttt = NaN;
                        CompNumDataCollect2{i, j, k} = [];
                        PhaseCollectMean2(i, j, k) = NaN;
                        RSSICollectMean2(i, j, k) = NaN;
                        CompNumCollectMean2(i, j, k) = NaN;                        
                        CountOneRead2 = CountOneRead2 + 1;
                    end
                end
                
                if PhaseCollectStd2(i, j, k) > phStdTh
                    % First completely reject the channel using phase
                    % standard deviation threshold.
                    PhaseCollectMean2(i, j, k) = NaN;
                    RSSICollectMean2(i, j, k) = NaN;
                    CompNumCollectMean2(i, j, k) = NaN;
                else
                    PhaseCollectMean2(i, j, k) = mean(phase_ttt);
                    RSSICollectMean2(i, j, k) = mean(rssi_ttt);
                    CompNumCollectMean2(i, j, k) = mean(CompNumDataCollect2{i, j, k});
                end
            else
                % Simply there's no data received for this channel.
                PhaseDataCollect2{i, j, k} = [];
                PhaseCollectStd2(i, j, k) = NaN;
                PhaseCollectMean2(i, j, k) = NaN;
                RSSIDataCollect2{i, j, k} = [];
                RSSICollectStd2(i, j, k) = NaN;
                RSSICollectMean2(i, j, k) = NaN;
                CompNumDataCollect2{i, j, k} = NaN;
                CompNumCollectMean2(i, j, k) = NaN;
            end                
        end
    end
end
fprintf('Number of channels with only one read: %3.0f \n', CountOneRead2);
fprintf('Percentage channels with only one read: %3.2f%% \n', 100*CountOneRead2/length(PhaseCollectStd2(:)));
clear chindexlist tagindexlist antennalist rssiimpinjlist rssiimpinjlist_d phasedeglist msgfreqlist

% Adjust the order in the frequency dimension.
PhaseDataCollect2 = PhaseDataCollect2(:, :, I_w);
PhaseCollectStd2 = PhaseCollectStd2(:, :, I_w);
PhaseCollectMean2 = PhaseCollectMean2(:, :, I_w);
RSSIDataCollect2 = RSSIDataCollect2(:, :, I_w);
RSSICollectStd2 = RSSICollectStd2(:, :, I_w);
RSSICollectMean2 = RSSICollectMean2(:, :, I_w);
CompNumDataCollect2 = CompNumDataCollect2(:, :, I_w);
CompNumCollectMean2 = CompNumCollectMean2(:, :, I_w);

% Pick out the data associated with the selected fixed frequencies.
PhaseDataCollect2 = PhaseDataCollect2(:, :, FreqInd);
PhaseCollectStd2 = PhaseCollectStd2(:, :, FreqInd);
PhaseCollectMean2 = PhaseCollectMean2(:, :, FreqInd);
RSSIDataCollect2 = RSSIDataCollect2(:, :, FreqInd);
RSSICollectStd2 = RSSICollectStd2(:, :, FreqInd);
RSSICollectMean2 = RSSICollectMean2(:, :, FreqInd);
CompNumDataCollect2 = CompNumDataCollect2(:, :, FreqInd);
CompNumCollectMean2 = CompNumCollectMean2(:, :, FreqInd);


idx_wo_pre = find(RSSICollectMean1 ~= RSSICollectMean1);
idx_w_pre = find(RSSICollectMean2 ~= RSSICollectMean2);
idx_pre = [idx_wo_pre; idx_w_pre];
idx_pre = unique(idx_pre);
NumLost_pre = length(idx_pre);
PercLost_pre = NumLost_pre/length(CompNumCollectMean1(:));
fprintf('Percentage Tx-Rx-Freq combination lost: %3.2f%% \n', 100*PercLost_pre);

% Reject the channels with blocked RSSI's.
countBlocked = 0;
if RJT == 1
    for n = 1:TagNum
        for k = 1:RecvNum
            for m = 1:FreqNumUse
                if RSSICollectMean2(n, k, m)/RSSICollectMean1(n, k, m) < rssiDiffTh
%                     disp(n); disp(k); disp(m);
                    RSSICollectMean1(n, k, m) = NaN;
                    PhaseCollectMean1(n, k, m) = NaN;
                    CompNumCollectMean1(n, k, m) = NaN;
                    RSSICollectMean2(n, k, m) = NaN;
                    PhaseCollectMean2(n, k, m) = NaN;
                    CompNumCollectMean2(n, k, m) = NaN;
                    countBlocked = countBlocked + 1;
                end
            end
        end
    end
end

G_wo_s = CompNumCollectMean1;
G_w_s = CompNumCollectMean2;


% Find the zero values in the calibration matrix, and replace them with the
% maximum value in the matrix. (We can replace those zeros with any number
% we like, actually.)
idx_wo = find(RSSICollectMean1 ~= RSSICollectMean1);
idx_w = find(RSSICollectMean2 ~= RSSICollectMean2);
idx = [idx_wo; idx_w];
idx = unique(idx);
% ConstParam = max(max(max(abs(G_wo_s))));
ConstParam = 0;
% G_wo_s(idx_wo) = ConstParam;
% G_w_s(idx_wo) = ConstParam;
% G_w_s(idx_w) = G_wo_s(idx_w);
G_wo_s(idx) = ConstParam;
G_w_s(idx) = ConstParam;

% How many channels are lost.
NumLost = length(idx);
PercLost = NumLost/length(G_wo_s(:));
fprintf('#LOS possibly blocked.. %d\n',countBlocked);
fprintf('Percentage Tx-Rx-Freq combination lost: %3.2f%% \n', 100*PercLost);


%% Calibration.
% Calculate the distance from tags to receivers.
r_x = bsxfun(@minus, tagPosition(:, 1), rxPosition(:, 1)');
r_y = bsxfun(@minus, tagPosition(:, 2), rxPosition(:, 2)');
r_z = bsxfun(@minus, tagPosition(:, 3), rxPosition(:, 3)');
r = sqrt(r_x.^2 + r_y.^2 + r_z.^2);

% Perform the calibration.
G_calib = zeros(TagNum, RecvNum, FreqNumUse);
switch(opts.calibType)
    case 1
        % Differential receiving: First calculate the ratio of the data
        % received by one receiver antenna, to the data received by all the
        % other receiver antennas. Perform this for the data from both the
        % case when no object is present, and the case when object is
        % present. Then calculate the ratio of the above-calculated ratios
        % of the case when object is present to the case when no object is
        % present. The results are stored in the matrix named G_calib.
        for m = 1:FreqNumUse
            for n = 1:TagNum
                for k = 1:RecvNum
                    for k_pair = 1:RecvNum
                        if k_pair ~= k
                            if G_w_s(n, k_pair, m) ~= 0 && G_wo_s(n, k_pair, m) ~= 0 ...
                                    && G_w_s(n, k, m) ~= 0 && G_wo_s(n, k, m) ~= 0
                                % Both denominators and both numerators
                                % should not be zero.
                                G_wo_s_ratio = G_wo_s(n, k, m)/G_wo_s(n, k_pair, m);
                                G_w_s_ratio = G_w_s(n, k, m)/G_w_s(n, k_pair, m);
                                G_calib(n, k, m) = G_calib(n, k, m) + (G_w_s_ratio/G_wo_s_ratio - 1) ...
                                    *exp(-1j*(2*pi*freq(m)/c)*r(n, k));
                            else
                                G_calib(n, k, m) = G_calib(n, k, m);
                            end
                        end
                    end
                end
            end
        end
        
        % Eliminate the weighting part in Tag phase. (calibration 3)
        % For each receiver, how the tag phase is changing. (calibration 4)
        % consider RF cable and clutter as a whole. (Guoyi) (calibration 5)
        % Instead of using receiver 1 as reference, what if we use another
        % receiver as reference?
    case 5
        % This is a different calibration from 2, 3, and 4. RF cable,
        % clutter, and tag circuitry are lumped into one variable to
        % represent the phase offset introduced by these sources. We expect
        % that by directly calculating the RF cable and clutter loss part,
        % instead of calculating the relative values (like in calibration
        % 2, 3, or 4), we are able to improve the results more.
%         scale = 1e-3;
%         G_w_s = scale.*G_w_s;
%         G_wo_s = scale.*G_wo_s;
        
        NonIdealCal = zeros(TagNum, RecvNum, FreqNumUse);
        for m = 1:FreqNumUse
            kn = 2*pi*freq(m)/c;
            for k = 1:RecvNum
                for n = 1:TagNum
                    distTagRecv = norm(tagPosition(n, :) - rxPosition(k, :));
                    CalPhs = -kn*distTagRecv;
                    CalCompx = cos(CalPhs) + 1j*sin(CalPhs);
%                     CalCompx = exp(1i*CalPhs);
                    NonIdealCal(n, k, m) = CalCompx*G_wo_s(n, k, m)';
                    if abs(NonIdealCal(n, k, m)) > 0
                        NonIdealCal(n, k, m) = NonIdealCal(n, k, m)/abs(NonIdealCal(n, k, m));
                    end
                    % We HAVE TO do normalization of the undesired part of
                    % the link, the "NonIdealCal" variable. If you don't
                    % normalize it, the results will be bad.
                end
            end
        end
        
        % Background subtraction. (We could lump this into the previous for
        % loop, but for ease of understanding, we will separate the
        % previous for loop, which is in charge of getting rid of the RF
        % cable, clutter and tag circuitry loss, and this one.)
        countBlockedAllFreq = 0;
        for n = 1:TagNum
            for k = 1:RecvNum
                G_w_abs = abs(G_w_s(n, k, :));
                G_wo_abs = abs(G_wo_s(n, k, :));
                
                if RJTAllFreq == 1
                    isLosBlocked = 1;
                    for m = 1:FreqNumUse
                        if G_w_abs(m)*G_wo_abs(m) > 0
                            if G_w_abs(m)/G_wo_abs(m) < rssiDiffThAllFreq
                                isLosBLocked = 1*isLosBlocked;
                            else
                                isLosBlocked = 0*isLosBlocked;
                            end
                        end
                    end
                elseif RJTAllFreq == 0
                    isLosBlocked = 0;
                end
                
                for m = 1:FreqNumUse
                    if G_w_abs(m)*G_wo_abs(m) > 0
                        if isLosBlocked == 0
                            % Equivalent implementation.
%                             G_calib(n, k, m) = G_w_s(n, k, m) - G_wo_s(n, k, m);
%                             G_calib(n, k, m) = G_calib(n, k, m)*NonIdealCal(n, k, m);

                            % Equivalent implementation.
                            G_w_s(n, k, m) = G_w_s(n, k, m)*NonIdealCal(n, k, m);
                            G_wo_s(n, k, m) = G_wo_s(n, k, m)*NonIdealCal(n, k, m);
                            G_calib(n, k, m) = G_w_s(n, k, m) - G_wo_s(n, k, m);
                            
                            % Equivalent implementation. (But this
                            % implementation is too complicated in
                            % expression.)
%                             G_w_s(n, k, m) = G_w_s(n, k, m)*(G_wo_s(n, k, m)'/abs(G_wo_s(n, k, m)))*exp(1i*(-(2*pi*freq(m)/c)*norm(tagPosition(n, :) - rxPosition(k, :))));
%                             G_wo_s(n, k, m) = G_wo_s(n, k, m)*(G_wo_s(n, k, m)'/abs(G_wo_s(n, k, m)))*exp(1i*(-(2*pi*freq(m)/c)*norm(tagPosition(n, :) - rxPosition(k, :))));
%                             G_wo_s(n, k, m) = (abs(G_wo_s(n, k, m))^2/abs(G_wo_s(n, k, m)))*exp(1i*(-(2*pi*freq(m)/c)*norm(tagPosition(n, :) - rxPosition(k, :))));
%                             G_calib(n, k, m) = G_w_s(n, k, m) - G_wo_s(n, k, m);
                            % This above implementation is equivalent to
                            % the previous two because NonIdealCal(n, m, k)
                            % = (G_wo_s(n, k, m)'/abs(G_wo_s(n, k, m)))*
                            % exp(-1i*(2*pi*freq(m)/c)*norm(tagPosition
                            % - rxPosition(k, :))). Therefore they are
                            % exactly the same implementations.

                        else
                            countBlockedAllFreq = countBlockedAllFreq + 1;
                        end
                    end
                end
            end
        end
        fprintf('#LOS possibly blocked.. %d\n',countBlockedAllFreq);
        % Questions: The variable "NonIdealCal" represents everything else
        % except the LoS propagation delay. 1) Why the propagation delay is
        % the distance from tag to receiver, not twice the distance from
        % tag to receiver? 2) Why doesn't it matter if we normalize
        % "NonIdealCal" or not, but it does matter if we normalize "G_wo_s"
        % or not?
        
        % Thoughts: We could first do background subtraction. Then we are
        % theoretically left with only the signal reflected by the object,
        % as well as RF cable loss, background clutter associated with the
        % link of object reflection, and the tag circuitry loss (all
        % complex numbers). There could be several ways to get rid of all
        % of these undesired parts and make only the object reflection
        % remain, including differential method, and non-differential
        % methods. Let me think about this now.
    case 6
        % We combine Calibration 2 with differential receiving, in hope for
        % eliminating initial phase associated with the received data, with
        % RF cable, clutter and tag circuitry loss subtracted.
    case 7
        % We combine Calibration 5 with differential receiving, in hope for
        % eliminating initial phase associated with the received data, with
        % RF cable, clutter and tag circuitry together subtracted.
    otherwise
        fprintf('Wrong calibration option selection, select a valid method. \n');
end
b = G_calib(:);
% fprintf('Percentage Tx-Rx-Freq combination lost: %3.2f%% \n', 100*(countBlockedAllFreq+NumLost)/length(G_wo_s(:)));

save('ReconsData.mat', 'b', 'personname', 'postureloc');
