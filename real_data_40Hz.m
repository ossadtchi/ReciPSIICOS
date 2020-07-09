%% Paths to the data, volume, model...
path_data = 'D:\Desktop\Ossadtchi\ITPC4Osadchii\clean_epoched';
path_channel_file = 'D:\Desktop\Ossadtchi\ITPC4Osadchii\head_model';
path_forw_model = 'D:\Desktop\Ossadtchi\ITPC4Osadchii\head_model';
bLoad = false;

%% Parameters block
GainSVDTh = 0.01; % 0.05 results into 47 eigensensors and makes it run faster but produces less contrasting subcorr scans
% for a more reliable performance use 0.01 to get all the sensor on board but be ready to wait;
PwrRnk = 1000; % Maximal rank of the projection
Rnk = 500; % Projection rank for the ReciPSIICOS beamformer
Rnkw = 500;
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

%% 1. Load forward model and reduce it
fprintf('\tLoading and reducing forward model\n\n')
% Load forward model
G3 = load(fullfile (path_forw_model, [subj '_G.mat']));
% Get grid node locations
R = G3.GridLoc;
% Load channel file
chans = load(fullfile(path_channel_file, [subj '_chans.mat']));

% Set to use gradiometers only
if isempty(Chan_ind)
    Chan_ind = find(strcmp({chans.Channel.Type}, 'MEG GRAD'));
end

[Nch, Nsites] = size(G3.Gain(Chan_ind,1:3:end));
G2d = zeros(Nch,Nsites*2);
G2d0 = zeros(Nch,Nsites*2);

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

%% ReciPSIICOS projection

Upwr = ProjectorOnlyAwayFromPowerComplete(G2d0U, PwrRnk,0);

%% Whitened ReciPSIICOS projection

% now find span of the target space
[UcorrRnk, Wps]= ProjectionAwayFromCorrelationWhitened(G2d0U, ...
    0, 0, PwrRnk, 0);
Pr = inv(Wps)*(eye(size(UcorrRnk(:,1:Rnkw),1))-UcorrRnk(:,1:Rnkw)*UcorrRnk(:,1:Rnkw)')*Wps;

%% 2. Load data, calculate data covariance
fprintf('\tLoading data and calculating covariance matrix\n\n')
data = load(fullfile(path_data, [subj '_data40Hz' LR '_noSMZE.mat']));
data = data.data40HzL_noSMZE;
Cov = zeros(Nch, Nch);
Nepoch = length(data.trial);
data_for_s = [];

[b,a] = butter(4, [f_low f_high]/(Fs/2));
% Filter and store data with relevant channels
for tr =1:Nepoch
    TF_data(:,:, tr)  = filtfilt(b,a, data.trial{tr}(Chan_ind, :)')';
    %TF_data(:,:, tr) = hilbert(temp_d')';
end;
avg_data = mean(TF_data, 3);

% Calculate covariance
Cov = avg_data*avg_data';
% Compute covariance in the virtual sensors space
Ca = UP*Cov'*UP';

data_for_bs = zeros(size(G3.Gain,1), size(data.trial{1}, 2));
data_for_bs(Chan_ind, :) = avg_data;

%% 3. Calculate beamformer weights
fprintf('\tCalculating the inverse solution\n\n')

% Beamformer
iCa = tihinv(Ca, 0.01);
Nsrc = size(G2dU,2);
Nch = size(G2dU,1);
range2d = 1:2;
clear w;
w = zeros(size(G2dU));
for i=1:Nsrc/2
    g = G2dU(:,range2d);
    num = g'*iCa;
    denum = inv(g'*iCa*g);
    w(:, range2d) = num'*denum;
    range2d = range2d+2;
end

% now apply ReciPSIICOS beamformer
%create projection matrix
% project the covariance matrix away from zero-phase synchrony subspace
Cap = reshape(Upwr(:,1:Rnk)*Upwr(:,1:Rnk)'*Ca(:),size(Ca));
% % fix negative eigenvalues (after projection) issue
[e a] = eig(Cap);
Cap = e*abs(a)*e';
% % compute pseudoinverse
iCap = tihinv(Cap, 0.01);
% %compute source-space scan
clear w_AP;
range2d = 1:2;
w_AP = zeros(size(G2dU));
for i=1:Nsrc/2
    g = G2dU(:,range2d);
    num = g'*iCap;
    denum = inv(g'*iCap*g);
    w_AP(:, range2d) = num'*denum;
    range2d = range2d+2;
end

% Whitened ReciPSIICOS beamformer
% create projection matrix
Cap =reshape( Pr*Ca(:), size(Ca));
[e a] = eig(Cap);
Cap = e*abs(a)*e';
iCap = tihinv(Cap, 0.01);

clear w_AP_w;
w_AP_w = zeros(size(G2dU));
range2d = 1:2;
for i=1:Nsites
    g = G2dU(:,range2d);
    num = g'*iCap;
    denum = inv(g'*iCap*g);
    w_AP_w(:, range2d) = num'*denum;
    range2d = range2d+2;
end

% wMNE
lambd = 0.1;
Cn = eye(size(G2dU,1));
GGt = G2dU*G2dU';
Wmne = G2dU'*inv(G2dU*G2dU'+lambd*trace(GGt)/size(GGt,1)*Cn);

%% 4. Reconstruct timeseries
% Reduce signal space
red_data = UP*avg_data;

% Apply beamformer weights in two orientations
reg_1 =  w(:,1:2:end)'*red_data;
reg_2 =  w(:,2:2:end)'*red_data;

% Reconstruct signal as sum of squares
reg_1 = abs(hilbert(reg_1')');
reg_2 = abs(hilbert(reg_2')');
reg_BF_ser = (reg_1).^2 + (reg_2).^2;

% ReciPSIICOS
AP_1 =  w_AP(:,1:2:end)'*red_data;
AP_2 =  w_AP(:,2:2:end)'*red_data;
AP_1 = abs(hilbert(AP_1')');
AP_2 = abs(hilbert(AP_2')');
AP_BF_ser1 = (AP_1).^2 + (AP_2).^2;

% Whitened ReciPSIICOS
AP_w_1 =  w_AP_w(:,1:2:end)'*red_data;
AP_w_2 =  w_AP_w(:,2:2:end)'*red_data;
AP_w_1 = abs(hilbert(AP_w_1')');
AP_w_2 = abs(hilbert(AP_w_2')');
AP_w_BF_ser1 = (AP_w_1).^2 + (AP_w_2).^2;

% MNE 
wMNE_ser = (Wmne(1:2:end,:)*red_data).^2 + (Wmne(2:2:end,:)*red_data).^2;
wMNE_1 =  Wmne(1:2:end, :)*red_data;
wMNE_2 =  Wmne(2:2:end, :)*red_data;
wMNE_1 = abs(hilbert(wMNE_1')');
wMNE_2 = abs(hilbert(wMNE_2')');
wMNE_ser1 = (wMNE_1).^2 + (wMNE_2).^2;

return;
%% 5. Visualize the results in Brainstorm
% For the next step we uploaded resulting data to brainstorm. This way we
% can plot activations on the cortex surface
time = 751;

figure; scatter([1:Nsites], AP_BF_ser1(:,time),[], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5])
figure; scatter([1:Nsites], reg_BF_ser(:,time),[], [0.85 0.325 0.098], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5])
figure; scatter([1:Nsites], wMNE_ser1(:,time),[], [0.7725 0.3804 0.902], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5])
figure; scatter([1:Nsites], AP_w_BF_ser1(:,time),[], [0.3 0.75 0.93], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', [0.001],'MarkerEdgeAlpha', [0.5])
xticklabels([])
yticklabels([])
namess = {'AP_pow', 'BF_pow', 'MNE_pow', 'AP_w'};
for ttt = 1:4;
    figure(ttt)
    xlim([0,Nsites+2])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end
