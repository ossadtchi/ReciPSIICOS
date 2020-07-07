%% Paths to the data, volume, model...
path_data = 'D:\Desktop\Ossadtchi\ITPC4Osadchii\clean_epoched';
path_forw_model = 'D:\Desktop\Ossadtchi\ITPC4Osadchii\head_model';

%% Parameters block
GainSVDTh = 0.001; % SVD Gain threshold. 0.05 results into 47 eigensensors and makes it run faster but produces less contrasting subcorr scans
% for a more reliable preformance use 0.01 to get all the sensor on board but be ready to wait;
PwrRnk = 1000; % Maximal rank of the projection
Rnk = 100; % Projection rank for the ReciPSIICOS beamformer
Rnkw = 950;
Chan_ind = []; % indexes of channels for source reconstruction (eg indexes of gradiometers)
% if empty: GRAD indexes will be found in brainstorm channel
% file
% Frequencies to filter
f_low = 35;
f_high = 45;
% Sampling rate
Fs = 1000;
%subj = 'RR015'; K0002
subj = 'RR015';
LR = 'L';
alpha = 0.1; % reg for baseline cov  matrix

%% 1. Load forward model
fprintf('\tLoading and reducing forward model\n\n')
% Load forward model
G3 = load(fullfile(path_forw_model, 'RR015_G.mat'));
% Get grid node locations
R = G3.GridLoc;
% Load channel file
chans = load(fullfile(path_forw_model, [subj '_chans.mat']));

% Set to use gradiometers only
if isempty(Chan_ind)
    Chan_ind = find(strcmp({chans.Channel.Type}, 'MEG GRAD'));
end

%% 2. Load data
fprintf('\tLoading data and calculating covariance matrix\n\n')
data = load(fullfile(path_data, [subj '_data40Hz' LR '_noSMZE.mat']));
data = data.data40HzL_noSMZE;
Nepoch = length(data.trial);
data_for_s = [];

[b,a] = butter(4, [f_low f_high]/(Fs/2));
% Filter and store data with relevant channels
for tr =1:Nepoch
    TF_data(:,:, tr)  = filtfilt(b,a, data.trial{tr}(Chan_ind, :)')';
end;
avg_data = mean(TF_data, 3);

bslt = find(data.time{1}<0);
bscov = avg_data(:, bslt)*avg_data(:, bslt)';

%% 3. Reduce forward model
[Nch, Nsites] = size(G3.Gain(Chan_ind,1:3:end));
G2d = zeros(Nch,Nsites*2);
G2d0 = zeros(Nch,Nsites*2);

G0 = G3;
range = 1:2;
for i=1:Nsites
    g = [G3.Gain(Chan_ind,1+3*(i-1)) G3.Gain(Chan_ind,2+3*(i-1)) G3.Gain(Chan_ind,3+3*(i-1))];
    [u sv v] = svd(g);
    gt = g*v(:,1:2);
    G2d(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
    G2d0(:,range) = gt;
    range = range + 2;
end;

% reduce the sensor space
[ug sg vg] = spm_svd(G2d*G2d',GainSVDTh);
UP = ug';
G2dU = UP*G2d;
G2d0U = UP*G2d0;

%% 4. Calculate data covariance

% Calculate covariance
Cov = avg_data*avg_data';
% Compute covariance in the virtual sensors space
Ca = UP*Cov'*UP';
iCa = tihinv(Ca, 0.01);

%% 5. Calculate weights
% Whitened projector
[UcorrRnk, Wps]= ProjectionAwayFromCorrelationWhitened(G2d0U, ...
    0, 1, PwrRnk, 0);
Pr = inv(Wps)*(eye(size(UcorrRnk(:,1:Rnkw),1))-UcorrRnk(:,1:Rnkw)*UcorrRnk(:,1:Rnkw)')*Wps;

% Beamformer
clear w;
w = zeros(size(G2dU));
range2d = 1:2;
for i=1:Nsites
    g = G2dU(:,range2d);
    num = g'*iCa;
    denum = inv(g'*iCa*g);
    w(:, range2d) = num'*denum;
    range2d = range2d+2;
end

% ReciPSIICOS
Upwr = ProjectorOnlyAwayFromPowerComplete(G2d0U, PwrRnk);
Cap = reshape(Upwr(:,1:Rnk)*Upwr(:,1:Rnk)'*Ca(:),size(Ca));
[e a] = eig(Cap);
Cap = e*abs(a)*e'; % fix negative eigenvalues
iCap = tihinv(Cap, 0.01);

clear w_AP;
w_AP = zeros(size(G2dU));
range2d = 1:2;
for i=1:Nsites
    g = G2dU(:,range2d);
    num = g'*iCap;
    denum = inv(g'*iCap*g);
    w_AP(:, range2d) = num'*denum;
    range2d = range2d+2;
end

% whitened ReciPSIICOS
wCap =reshape( Pr*Ca(:), size(Ca));
[e a] = eig(wCap);
wCap = e*abs(a)*e';
iwCap = tihinv(wCap, 0.01);

clear w_AP_w;
w_AP_w = zeros(size(G2dU));
range2d = 1:2;
for i=1:Nsites
    g = G2dU(:,range2d);
    num = g'*iwCap;
    denum = inv(g'*iwCap*g);
    w_AP_w(:, range2d) = num'*denum;
    range2d = range2d+2;
end

% MNE
lambd = 0.1;
Cn = eye(size(G2dU,1));
GGt = G2dU*G2dU';
Wmne = G2dU'*inv(G2dU*G2dU'+lambd*trace(GGt)/size(GGt,1)*Cn);


%% Calculate best orientations for the weights

clear w_fixed wp_fixed wmne_fixed
range2d = 1:2;
for i=1:Nsites
    % beamformer
    w_t = w(:, range2d);
    [vec eigv] = eig(w_t'*Ca*w_t);
    [~, mai] = max(diag(eigv));
    w_fixed(:, i) = (w_t*vec(:,mai))';
    
    % recipsiicos
    wp_t = w_AP(:, range2d);
    [vec eigv] = eig(wp_t'*Cap*wp_t);
    [~, mai] = max(diag(eigv));
    wp_fixed(:, i) = (wp_t*vec(:,mai))';
    
    % whitened recipsiicos
    wpw_t = w_AP_w(:, range2d);
    [vec eigv] = eig(wpw_t'*wCap*wpw_t);
    [~, mai] = max(diag(eigv));
    wpw_fixed(:, i) = (wpw_t*vec(:,mai))';
    
    % MNE
    wmne_t = Wmne(range2d,:);
    [vec eigv] = eig(wmne_t*Ca*wmne_t');
    [~, mai] = max(diag(eigv));
    wmne_fixed(:, i) = (wmne_t'*vec(:,mai))';
    range2d = range2d+2;
end

%% Reconstruct sources
red_data = UP*avg_data;
reg_BF_ser = w_fixed'*red_data;
AP_BF_ser1 = wp_fixed'*red_data;
AP_w_BF_ser1 = wpw_fixed'*red_data;
wMNE_ser1 = wmne_fixed'*red_data;

return;
%% PLOT
time = 751;

figure
subplot(2,2,1); scatter([1:Nsites], AP_BF_ser1(:,time).^2,[], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5]);
title('RP'); xlim([0,Nsites+2])
subplot(2,2,2); scatter([1:Nsites], reg_BF_ser(:,time).^2,[], [0.85 0.325 0.098], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5])
title('BF'); xlim([0,Nsites+2])
subplot(2,2,4); scatter([1:Nsites], wMNE_ser1(:,time).^2,[], [0.7725 0.3804 0.902], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5])
title('MNE'); xlim([0,Nsites+2])
subplot(2,2,3); scatter([1:Nsites], AP_w_BF_ser1(:,time).^2,[], [0.3 0.75 0.93], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5])
title('wRP'); xlim([0,Nsites+2])
