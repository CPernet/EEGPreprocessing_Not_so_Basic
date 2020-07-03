%% OHBM2020 EEG Preprocessing: not so basic after all
%
% *This script illustrate the effect of filtering, IC labeling and
% referencing on statistical results.* 
% The data are taken from https://openneuro.org/datasets/ds002718/versions/1.0.2
% We use here sub-011 only.
% The code was writen by Cyril Pernet, reusing Nicolas Langers' code from
% <https://osf.io/z8uqx/ OSF> and snippet from LI Dong.

%% let's check the environement, files, etc

clc
clear variables
current = pwd;  %% assuming we are in \code
addpath([current filesep 'NT_tools']);
addpath([current filesep 'local_functions']);

% check all the tools we need are here
if ~exist('eeglab.m','file')
    error('eeglab is not in your matlab path')
else
    eeglab; % allows loading plugins
end

if ~exist('iclabel.m','file')
    error('IClabewl plugin is missing, please install')
end

if ~exist('pop_viewprops.m','file')
    error('ViewProps plugin is missing, please install')
end

if ~exist('rest_refer.m','file')
    error('REST plugin is missing, please install')
end

if exist('limo_eeg.m','file')
    root = fileparts(which('limo_eeg'));
    addpath([root filesep 'limo_cluster_functions']);
    addpath([root filesep 'external' filesep 'color_maps']);
    diverging_bwr = load([root filesep 'external' filesep 'color_maps' filesep 'diverging_bwr.mat']);
    diverging_bwr = diverging_bwr.dmap;
else
    error('LIMO tools plugin is missing, please install')
end

% locate the data
datafolder = [fileparts(current) filesep 'data' filesep 'sub-011' filesep 'eeg'];
if ~exist(fullfile(datafolder,'sub-011_task-FaceRecognition_eeg.set'),'file')
    warndlg('Can''t locate the data - select the data file')
    [datafile,datafolder] = uigetfile('*.set','select a subjects'' eeg file');
    if isequal(datafile,0) || isempty(datafile,0)
        disp('Selection cancelled'); return
    else
        datafile  = fullfile(datafolder, datafile);
    end
else
    datafile = fullfile(datafolder,'sub-011_task-FaceRecognition_eeg.set');
end


%% we start by loading the data and removing EKG, HEOG, VEOG
EEG = pop_loadset(datafile);
EEG = pop_select(EEG, 'nochannel',{'EEG061','EEG062','EEG063','EEG064'});

figure('Name','Raw data topography'); 
subplot(1,2,1);
topoplot(zscore(double(EEG.data(:,100))),EEG.chanlocs,'electrodes','off','colormap',diverging_bwr); 
title('Original file - data frame 100');subplot(1,2,2)
pop_spectopo(EEG, 1, [0  size(EEG.data,2)],'EEG','percent',15,...
    'freq',[6 10 22 50],'freqrange',[2 70],'electrodes','off');
set(gcf,'Colormap',diverging_bwr);
drawnow

%% bad channel indentification 
% we use here cleanraw_data.m comparing results with or witout 1Hz temporary filter

% create the copy and remove bad channels (includes a filter at 0.5Hz)
EEG_badchan = EEG;
EEG_badchan = pop_clean_rawdata(EEG_badchan,'FlatlineCriterion',5,'ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,...
    'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
    'Distance','Euclidian','WindowCriterionTolerances','off' );

% Temporary highpass filter 1Hz 
EEGf_badchan = EEG;
EEGf_badchan = pop_eegfiltnew(EEGf_badchan,1,0);
EEGf_badchan = pop_clean_rawdata(EEGf_badchan,'FlatlineCriterion',5,'ChannelCriterion',0.8,...
    'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,...
    'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
    'Distance','Euclidian','WindowCriterionTolerances','off' );

% visualize and compare results (new first, old second)
vis_artifacts(EEG_badchan, EEG)
vis_artifacts(EEGf_badchan, EEG)

% find index of the retained channels and remove from original data
tbl_channels = struct2table(EEGf_badchan.chanlocs); 
ind_retained = find(arrayfun(@(x)ismember(x.labels,tbl_channels.labels),EEG.chanlocs));
ind_bad      = setdiff(1:size(EEG.chanlocs,1),ind_retained);
EEG          = pop_select(EEG, 'nochannel',{EEG.chanlocs(ind_bad).labels});

%% <https://www.sciencedirect.com/science/article/pii/S1053811919309474 Zapline> to remove Power line artifacts
% Empirical tests showed that 7 components are required to sufficiently 
% remove the line noise - note the spatial distribution of PSD values stays
% the same at the different frequencies, except a small change posteriorly
% at the ~22Hz <https://en.wikipedia.org/wiki/Harmonic sub-harmonic> frequency

% Notch filter (better not used) would be:
% EEG = pop_eegfiltnew(EEG, 48,52,826,1,[],0); 

% test Zapline with different number of components:
line_noise_freq = 50; % 50Hz line noise
FLINE           = line_noise_freq/EEG.srate; % line frequency
p.nfft          = 1024;
ncomp           = [3 5 7];

figure('Name','Power line removal - Zapline 3-5-7 components', ...
    'units','normalized','outerposition',[0 0 1 1])
for n = 3:-1:1
    [clean, noise] = nt_zapline(EEG.data',FLINE,ncomp(n));
    EEGzap{n}      = EEG; 
    EEGzap{n}.data = clean';
    subplot(1,3,n);
    pop_spectopo(EEGzap{n}, 1, [0  size(EEGzap{n}.data,2)],'EEG','percent',15,...
        'freq', [6 10 22 line_noise_freq],'freqrange',[2 70],'electrodes','off');
    set(gcf,'Colormap',diverging_bwr);
end

% Plot the the raw vs clean/removed signal (here for 7 components):
 
Power_cleanEEG = EEGzap{3}; % this is the right data
clear EEGzap                % clear up memory

% Compute the PSD of each channel and add them (nt_dpect_plot computes
% the Power Spectral Density estimate via <https://en.wikipedia.org/wiki/Welch%27s_method#:~:text=Welch's%20method%2C%20named%20after%20Peter,a%20signal%20at%20different%20frequencies. Welch's method>) 
% Note here the normalization by sqrt(mean(EEG.data(:)).^2), avg power of
% the raw data
fprintf('proportion of non-DC power removed %g\n:', ...
    nt_wpwr(EEG.data-Power_cleanEEG.data)/nt_wpwr(nt_demean(EEG.data)));
[psd_raw ,fr] = nt_spect_plot((EEG.data/sqrt(mean(EEG.data(:).^2)))',p.nfft,[],[],1/FLINE);
[psd_zap ,fz] = nt_spect_plot((Power_cleanEEG.data/sqrt(mean(EEG.data(:).^2)))',p.nfft,[],[],1/FLINE);
[psd_diff,fd] = nt_spect_plot(((EEG.data-Power_cleanEEG.data)/sqrt(mean(EEG.data(:).^2)))',p.nfft,[],[],1/FLINE);

figure('Name','Raw vs Zapline corrected PSD');
subplot(2,1,1); semilogy(fr,abs(psd_raw)/sum(psd_raw),'LineWidth',2);
hold on; grid on; semilogy(fz,abs(psd_zap)/sum(psd_zap),'--','LineWidth',2);
legend('raw','clean'); legend boxoff; title('raw vs clean')
ylabel('relative power'); subplot(2,1,2);
semilogy(fd,abs(psd_diff)/sum(psd_diff),'k','LineWidth',2);
xlabel('frequency (relative to line)');
ylabel('relative power'); subplot(2,1,2);
title('PDS difference'); grid on

% we are using the results from ZapLine here
clear EEG; EEG = Power_cleanEEG; clear Power_cleanEEG

%% Let's look at the effect of high-pass and low-pass filters on artefacts
% We will filter using default parameters from EEGLab, compute ICA and use
% IClabel to classify artefacts.

IClabel_selection = repmat([0.8 1],7,1); % prob. [80% 100%]
IClabel_selection(1,:) = NaN;            % do not flag brain

% The recommendation for ERP is a high-pass filter between 0.01 to 0.05Hz
% and no low-pass filter

if ~exist('all_filtICA.mat','file')
    % _passband edge 0.02Hz (cut-off 0.01 hz) no low-pass_
    EEG_002hz         = pop_eegfiltnew(EEG,0.02,0);
    EEG_ica_002hz     = pop_runica(EEG_002hz, 'icatype', 'runica', 'concatcond','on',...
        'extended',1,'interrupt','on','reorder','on','options',{'pca',EEG.nbchan-1});

    % _lets' see what is happening adding a low_pass 39Hz (cut-off 43 hz)_
    EEG_002hz_39hz         = pop_eegfiltnew(EEG_002hz,0,39);
    EEG_ica_002hz_39hz     = pop_runica(EEG_002hz_39hz, 'icatype', 'runica', 'concatcond','on',...
        'extended',1,'interrupt','on','reorder','on','options',{'pca',EEG.nbchan-1});
else
    load('all_filtICA.mat'); % all the EEG* with ica (from below as well)
end

% compute classification and plot IC with labels
EEG_ica_002hz_ICL = iclabel(EEG_ica_002hz);
pop_viewprops(EEG_ica_002hz_ICL,0,1:20);
set(gcf,'color','w','Colormap',diverging_bwr);
pop_prop_extended(EEG_ica_002hz_ICL,0,1)
set(gcf,'color','w','Colormap',diverging_bwr);

%%
% if we are happy with this result, the clean data are
EEG_ica_002hz_ICL_clean = pop_icflag(EEG_ica_002hz_ICL,IClabel_selection);
EEG_ica_002hz_ICL_clean = pop_subcomp(EEG_ica_002hz_ICL_clean, ...
    find(EEG_ica_002hz_ICL_clean.reject.gcompreject), 0);

% redo with low-pass filtered data
EEG_ica_002hz_39hz_ICL = iclabel(EEG_ica_002hz_39hz);
pop_viewprops(EEG_ica_002hz_39hz_ICL,0,1:20);
set(gcf,'color','w','Colormap',diverging_bwr);
pop_prop_extended(EEG_ica_002hz_ICL,0,1)
set(gcf,'color','w','Colormap',diverging_bwr);

% _what can be said of using low-pass in ICA and labelling?_

    
%%
% _Let's test different filters_
if ~exist('all_filtICA.mat','file')
    % high-pass filter passband edge 1hz (cut-off 0.5 hz) no low-pass
    EEG_1hz = pop_eegfiltnew(EEG,1,0);
    EEG_ica_1hz  = pop_runica(EEG_1hz, 'icatype', 'runica', 'concatcond','on',...
        'extended',1,'interrupt','on','reorder','on','options',{'pca',EEG.nbchan-1});

    % add low_pass 39Hz (cut-off 43 hz)
    EEG_1hz_39hz      = pop_eegfiltnew(EEG_1hz,0,39);
    EEG_ica_1hz_39hz  = pop_runica(EEG_1hz_39hz, 'icatype', 'runica', 'concatcond','on',...
        'extended',1,'interrupt','on','reorder','on','options',{'pca',EEG.nbchan-1});

    % high_pass filter passband edge 3Hz (cut-off 2 hz) no low-pass
    EEG_2hz      = pop_eegfiltnew(EEG,3,0);
    EEG_ica_2hz  = pop_runica(EEG_2hz, 'icatype', 'runica', 'concatcond','on',...
        'extended',1,'interrupt','on','reorder','on','options',{'pca',EEG.nbchan-1});

    % add low_pass 39Hz (cut-off 43 hz)
    EEG_2hz_39hz      = pop_eegfiltnew(EEG_2hz,0,39);
    EEG_ica_2hz_39hz  = pop_runica(EEG_2hz_39hz, 'icatype', 'runica', 'concatcond','on',...
        'extended',1,'interrupt','on','reorder','on','options',{'pca',EEG.nbchan-1});

    % save all that computational time for next time we want to run the tutorial
    save all_filtICA -v7.3
end

% classify and make figures of the different filtered version
EEG_ica_1hz_ICL = iclabel(EEG_ica_1hz);
pop_viewprops(EEG_ica_1hz_ICL,0,1:20);
set(gcf,'color','w','Colormap',diverging_bwr);
EEG_ica_1hz_39hz_ICL = iclabel(EEG_ica_1hz_39hz);
pop_viewprops(EEG_ica_1hz_39hz_ICL,0,1:20);
set(gcf,'color','w','Colormap',diverging_bwr);
EEG_ica_2hz_ICL = iclabel(EEG_ica_2hz);
pop_viewprops(EEG_ica_2hz_ICL,0,1:20);
set(gcf,'color','w','Colormap',diverging_bwr);
pop_prop_extended(EEG_ica_2hz_ICL,0,1)
set(gcf,'color','w','Colormap',diverging_bwr);
EEG_ica_2hz_39hz_ICL = iclabel(EEG_ica_2hz_39hz);
pop_viewprops(EEG_ica_2hz_39hz_ICL,0,1:20);
set(gcf,'color','w','Colormap',diverging_bwr);
  
%% Temporary filtering
% *We can see that low frequency have more weights on ICA*
% The 1st component for filtered data at 0.01Hz is classified as 'others'
% with a power spectrum around 0, although it still have some 'brain' class
% associated to it. ICA/IClabel seem to work best with at least 1Hz filter,
% possibly adding a low-pass at 40Hz. At the same time, we want data 
% filtered at 0.05Hz. A solution is temporary filtering, that is filter
% the data at say 2Hz/40Hz, compute ICA and labelling. We can then remove
% the artefacts, i.e. components identified as such and backproject onto 
% the 0.01Hz we want.

% we use  for backprojection of ICA weight matrix and ICA sphere matrix
EEG_ica_2hz_39hz  = keepICA(EEG_ica_2hz_39hz); 
EEG_ica_bpr_002hz = ica_foreign_backproject(EEG_002hz,EEG_ica_2hz_39hz);

% let's see what IClabel has to say about the back projection
% the low frequency component is gone (unsurprizingly)
EEG_ica_bpr_002hz_ICL = iclabel(EEG_ica_bpr_002hz);
pop_viewprops(EEG_ica_bpr_002hz_ICL,0,1:20);
set(gcf,'color','w','Colormap',diverging_bwr);

% automatically keep Brain components 
EEG_ica_bpr_002hz_ICL_clean = pop_icflag(EEG_ica_bpr_002hz_ICL, IClabel_selection); 
EEG_ica_bpr_002hz_ICL_clean = pop_subcomp(EEG_ica_bpr_002hz_ICL_clean,...
    find(EEG_ica_bpr_002hz_ICL_clean.reject.gcompreject) , 0);

% Sanity check 
% we computed ICA for data at 0.2Hz/40Hz, back projected for the 0.01Hz
% filtered data - re-ran IClabel, flag brain components, and projected them
% onto the scalp creating a clean dataset. All components should now be
% brain

pop_prop_extended(EEG_ica_bpr_002hz_ICL_clean, 0, 1:10)

%% Re-referencing
% we have now two clean datasets, EEG_ica_002hz_ICL_clean and
% EEG_ica_bpr_002hz_ICL_clean which have the power line removed with
% ZapLine, then filtered at 0.01Hz (0.02Hz bound) and artefacts removed 
% using IClabel ; directly from the 0.01Hz data or using labeling of the
% 2Hz/40Hz filtered data. We can now re-reference, we will use the common
% average refrence and rest.

EEG_ica_002hz_ICL_clean_ca     = pop_reref(EEG_ica_002hz_ICL_clean,[],'interpchan',[]);
EEG_ica_bpr_002hz_ICL_clean_ca = pop_reref(EEG_ica_bpr_002hz_ICL_clean,[],'interpchan',[]);

% REST parameters
channel_locs = cell2mat(arrayfun(@(x) [x.X x.Y x.Z],EEG.chanlocs,'UniformOutput',false));
xyz_dipoles  = load(fullfile(fileparts(which('pop_REST_reref.m')),'corti869-3000dipoles.dat'));
% Calculate the dipole orientations.
xyz_dipOri = bsxfun ( @rdivide, xyz_dipoles, sqrt ( sum ( xyz_dipoles .^ 2, 2 ) ) );
xyz_dipOri ( 2601: 3000, 1 ) = 0;
xyz_dipOri ( 2601: 3000, 2 ) = 0;
xyz_dipOri ( 2601: 3000, 3 ) = 1;
% define headmodel
headmodel        = [];
headmodel.type   = 'concentricspheres';
headmodel.o      = [ 0.0000 0.0000 0.0000 ];
headmodel.r      = [ 0.8700,0.9200,1];
headmodel.cond   = [ 1.0000,0.0125,1];
headmodel.tissue = { 'brain' 'skull' 'scalp' };
% calculate leadfield
[G,~] = dong_calc_leadfield3(channel_locs,xyz_dipoles,xyz_dipOri,headmodel);
% compute
EEG_ica_002hz_ICL_clean_rest          = EEG_ica_002hz_ICL_clean_ca;
EEG_ica_002hz_ICL_clean_rest.data     = rest_refer(detrend(EEG_ica_002hz_ICL_clean_ca.data,'constant'),G');
EEG_ica_bpr_002hz_ICL_clean_rest      = EEG_ica_bpr_002hz_ICL_clean_ca;
EEG_ica_bpr_002hz_ICL_clean_rest.data = rest_refer(detrend(EEG_ica_bpr_002hz_ICL_clean_ca.data,'constant'),G');

%% Check ERPs and statistics

% Epoch data, remove baseline
EEG_ica_002hz_ICL_clean_ca = pop_epoch(EEG_ica_002hz_ICL_clean_ca,...
    {'famous_new','famous_second_early','famous_second_late',...
    'scrambled_new','scrambled_second_early','scrambled_second_late',...
    'unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},...
    [-0.2 1] ,'epochinfo','yes');
pop_rmbase(EEG_ica_002hz_ICL_clean_ca,[-200 0]);
EEG_ica_002hz_ICL_clean_ca = eeg_checkset(EEG_ica_002hz_ICL_clean_ca);

EEG_ica_bpr_002hz_ICL_clean_ca = pop_epoch(EEG_ica_bpr_002hz_ICL_clean_ca,...
    {'famous_new','famous_second_early','famous_second_late',...
    'scrambled_new','scrambled_second_early','scrambled_second_late',...
    'unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},...
    [-0.2 1] ,'epochinfo','yes');
pop_rmbase(EEG_ica_bpr_002hz_ICL_clean_ca,[-200 0]);
EEG_ica_bpr_002hz_ICL_clean_ca = eeg_checkset(EEG_ica_bpr_002hz_ICL_clean_ca);

EEG_ica_002hz_ICL_clean_rest = pop_epoch(EEG_ica_002hz_ICL_clean_rest,...
    {'famous_new','famous_second_early','famous_second_late',...
    'scrambled_new','scrambled_second_early','scrambled_second_late',...
    'unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},...
    [-0.2 1] ,'epochinfo','yes');
pop_rmbase(EEG_ica_002hz_ICL_clean_rest,[-200 0]);
EEG_ica_002hz_ICL_clean_rest = eeg_checkset(EEG_ica_002hz_ICL_clean_rest);

EEG_ica_bpr_002hz_ICL_clean_rest = pop_epoch(EEG_ica_bpr_002hz_ICL_clean_rest,...
    {'famous_new','famous_second_early','famous_second_late',...
    'scrambled_new','scrambled_second_early','scrambled_second_late',...
    'unfamiliar_new','unfamiliar_second_early','unfamiliar_second_late'},...
    [-0.2 1] ,'epochinfo','yes');
pop_rmbase(EEG_ica_bpr_002hz_ICL_clean_rest,[-200 0]);
EEG_ica_bpr_002hz_ICL_clean_rest = eeg_checkset(EEG_ica_bpr_002hz_ICL_clean_rest);

% Compute whole scalp statistics using LIMO tools
% 1st we find the different type of faces 
% 2nd create a design matrix X such Y = BX + e
% 3rd estimate the face effect using the GLM
conditions = arrayfun(@(x) x.eventface_type{1},EEG_ica_002hz_ICL_clean_ca.epoch,...
    'UniformOutput',false);
condition_names = unique(conditions);
DesignMatrix = [zeros(length(conditions),3) ones(length(conditions),1)];
for c = 1:3
    DesignMatrix(:,c) = cellfun(@(x) strcmp(x,condition_names{c}), conditions);
end

% let's illustrate differences on a channel
channel = 64;
cl = [1 0 0; 0 1 0; 0 0 1];
figure('Name','ERP channel 63');
for condition = 1:3
    subplot(1,2,1);
    plot(EEG_ica_002hz_ICL_clean_ca.times,...
        mean(squeeze(EEG_ica_002hz_ICL_clean_ca.data(channel,:,logical(DesignMatrix(:,condition)))),2),...
        'LineWidth',2,'Color',cl(condition,:));hold on
    plot(EEG_ica_002hz_ICL_clean_ca.times,...
        mean(squeeze(EEG_ica_002hz_ICL_clean_rest.data(channel,:,logical(DesignMatrix(:,condition)))),2),...
        '--','LineWidth',2,'Color',cl(condition,:));
    subplot(1,2,2);
    plot(EEG_ica_002hz_ICL_clean_ca.times,...
        mean(squeeze(EEG_ica_bpr_002hz_ICL_clean_ca.data(channel,:,logical(DesignMatrix(:,condition)))),2),...
        'LineWidth',2,'Color',cl(condition,:));hold on
    plot(EEG_ica_002hz_ICL_clean_ca.times,...
        mean(squeeze(EEG_ica_bpr_002hz_ICL_clean_rest.data(channel,:,logical(DesignMatrix(:,condition)))),2),...
        '--','LineWidth',2,'Color',cl(condition,:));
end
subplot(1,2,1); title('standard IClabel');grid on; box on
subplot(1,2,2); title('IClabel with temporary filtering');grid on; box on
drawnow

%% 
% run GLM and F value and p-value maps
for channel = size(EEG_ica_002hz_ICL_clean_ca.data,1):-1:1
    model = limo_glm(squeeze(EEG_ica_002hz_ICL_clean_ca.data(channel,:,:))',...
        DesignMatrix,3,0,0,'OLS','Time',[],[]);
    Fmap1(channel,:) = model.conditions.F;
    Pmap1(channel,:) = model.conditions.p;
    model = limo_glm(squeeze(EEG_ica_bpr_002hz_ICL_clean_ca.data(channel,:,:))',...
        DesignMatrix,3,0,0,'OLS','Time',[],[]);
    Fmap2(channel,:) = model.conditions.F;
    Pmap2(channel,:) = model.conditions.p;
    model = limo_glm(squeeze(EEG_ica_002hz_ICL_clean_rest.data(channel,:,:))',...
        DesignMatrix,3,0,0,'OLS','Time',[],[]);
    Fmap3(channel,:) = model.conditions.F;
    Pmap3(channel,:) = model.conditions.p;
    model = limo_glm(squeeze(EEG_ica_bpr_002hz_ICL_clean_rest.data(channel,:,:))',...
        DesignMatrix,3,0,0,'OLS','Time',[],[]);
    Fmap4(channel,:) = model.conditions.F;
    Pmap4(channel,:) = model.conditions.p;
end


% we can also check across all channels
figure('Name','Stat results (uncorrected)') 
subplot(3,3,1);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap1,1),Fmap1.*(Pmap1<.05)); 
set(gcf,'color','w','Colormap',diverging_bwr(129:end,:));
xlabel('time (ms)'); ylabel('channels');
title('Face effect IC with Filter 0.01Hz Ref: common')
subplot(3,3,2);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap2,1),Fmap2.*(Pmap2<.05)); 
set(gcf,'color','w','Colormap',diverging_bwr(129:end,:));
xlabel('time (ms)'); ylabel('channels');
title('Face effect IC with Filter 2Hz Ref: common')
subplot(3,3,3);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap2,1),Pmap1-Pmap2)
colormap(gca, diverging_bwr);
xlabel('time (ms)'); ylabel('channels');
title('p-value differences')

subplot(3,3,4);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap3,1),Fmap3.*(Pmap3<.05)); 
set(gcf,'color','w','Colormap',diverging_bwr(129:end,:));
xlabel('time (ms)'); ylabel('channels');
title('Face effect IC with Filter 0.01Hz Ref: rest')
subplot(3,3,5);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap4,1),Fmap4.*(Pmap4<.05)); 
set(gcf,'color','w','Colormap',diverging_bwr(129:end,:));
xlabel('time (ms)'); ylabel('channels');
title('Face effect IC with Filter 0.01Hz Ref: rest')
subplot(3,3,6);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap3,1),Pmap3-Pmap4)
colormap(gca, diverging_bwr);
xlabel('time (ms)'); ylabel('channels');
title('p-value differences')

subplot(3,3,7);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap1,1),Pmap1-Pmap3)
colormap(gca, diverging_bwr);
xlabel('time (ms)'); ylabel('channels');
title('p-value differences')
subplot(3,3,8);
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap2,1),Pmap2-Pmap4)
colormap(gca, diverging_bwr);
xlabel('time (ms)'); ylabel('channels');
title('p-value differences')

subplot(3,3,9);
all = zeros(size(Fmap1,1),size(Fmap1,2),4);
all(:,:,1) = Fmap1;
all(:,:,2) = Fmap2;
all(:,:,3) = Fmap3;
all(:,:,4) = Fmap4;
imagesc(EEG_ica_002hz_ICL_clean_ca.times,1:size(Fmap2,1),std(all,[],3))
colormap(gca, diverging_bwr);
xlabel('time (ms)'); ylabel('channels');
title('F-values vibration (std among methods)')

