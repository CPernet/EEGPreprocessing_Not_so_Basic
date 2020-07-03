%% OHBM Educational Course - Effects of filtering on automated ICA classifiers, test ZapLine, different high-pass filters and different baseline correction
%
% This tutorial demostrates on one subjects (subject #2):
% - ZapLine (to remove line noise) removing different number of components
% - bad channel detection (temporary filtering 1Hz): cleanraw_data
% - different high-pass 
% - filters (using pop_eegfiltnew) and applies them to ICLabel (including temporary filtering)
% - ICA (computationally expensive)
% - ICLabel with different filtered data
% - Interpolation of bad channels
% - epoch data and apply different baseline correction (200ms vs. 100ms)
% - plot ERP ('famous faces') and ERP differences ('famous faces' - 'scrambled image')
% 
%% Different high-pass filters
% high_pass 0.02Hz (cut-off 0.01 hz) no low-pass
% high_pass 0.2Hz (cut-off 0.1 hz) no low-pass
% high_pass 0.5Hz (cut-off 0.25 hz) no low-pass
% high_pass 1Hz  (cut-off 0.5 hz) no low-pass
% high_pass 2Hz  (cut-off 1 hz) no low-pass
% high_pass 3Hz  (cut-off 1 hz) no low-pass

%% added low-pass filters
% high_pass 0.2Hz (cut-off 0.1 hz) +  low_pass 39Hz (cut-off 43 hz)
% high_pass 3Hz (cut-off 2 hz) +  low_pass 39Hz (cut-off 43 hz)
%
%% computes ICA for:
% high_pass 0.2Hz (cut-off 0.1 hz) no low-pass
% high_pass 1Hz  (cut-off 0.5 hz) no low-pass
% high_pass 2Hz  (cut-off 1 hz) no low-pass
% high_pass 3Hz  (cut-off 1 hz) no low-pass
% high_pass 0.2Hz (cut-off 0.1 hz) +  low_pass 39Hz (cut-off 43 hz)
% high_pass 3Hz (cut-off 2 hz) +  low_pass 39Hz (cut-off 43 hz)
%
%
% The code is purposely written in the simplest way possible for beginners.
% written by Nicolas Langer (n.langer@psychologie.uzh.ch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% You have to specificy the path to the course folder (scripts), the path to EEGLab and the path to the data 

clc
clear all

%% set the path to the OHBM educational folder

% restoredefaultpath % set all pathes back to default

%% specify the path where you stored the OHBM educational course folder
course_folder =  
% e.g. course_folder = '~/Dropbox/AA_transfer/OHBM_educational_course/final_folder/'

% addpath('~/Dropbox/AA_transfer/ICA_test/')
addpath(genpath(course_folder))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize EEGLAB (so you can use the eeglab function from your command line):

%% cd to the path where you have saved your eeglab 2019 version
 cd 
% e.g. cd ~/Dropbox/EEG_analysis/GeneralMatlab/eeglab2019_1

eeglab
close % close the gui from EEGlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the EEG data from subject #2

%% specify the path to the subject #2 from the Wakeman and Henson data set 
path_sub2 =  
% e.g. path_sub2 = '/Volumes/methlab_data/Wakeman_Henson_EEG/sub-002/eeg/'

EEG = pop_loadset('filename','sub-002_task-FaceRecognition_eeg.set','filepath',path_sub2);


% eeglab redraw % if you like to see the subject loaded in the EEGLAB GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%% import BIDS ( to load the data for entire study)
% filepath        = 'XXXX\WakemanHenson_Faces\eeg';
% [STUDY, ALLEEG] = pop_importbids(filepath, 'bidsevent','on','bidschanloc','on', 'studyName','Face_detection');
% ALLEEG = pop_select( ALLEEG, 'nochannel',{'EEG061','EEG062','EEG063','EEG064'});
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
% % reorient
% EEG = pop_chanedit(EEG,'nosedir','+Y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% !! Change Coordinates (Problem: in the original EEG file from Wakeman & Henson the EEG coordinates were rotated).
% If you have downloaded version 1.0.2
% (https://openneuro.org/datasets/ds002718/versions/1.0.2) you don't need
% to care about this the electrodes are already rotated.

% Sanity Check: topoplot (should display a strong positivity in the frontal electrodes)
figure; topoplot(double(EEG.data(:,100)),EEG.chanlocs,'electrodes','labels'); 
title('original file')
EEG_orig = EEG;

%% in case rotation of electrodes is required
%EEG=pop_chanedit(EEG, 'nosedir','+Y');
% 
% figure; topoplot(double(EEG.data(:,100)),EEG.chanlocs,'electrodes','labels'); 
% title('rotated file +Y')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Remove the non-EEG electrodes
% EEG061	HEOG	n/a
% EEG062	VEOG	n/a
% EEG063	EKG	n/a
% EEG064	EKG	n/a

el_excl = {'EEG061' 'EEG062' 'EEG063' 'EEG064'}

ind_excl = get_el_index(EEG,el_excl)

EEG = pop_select( EEG,'nochannel',ind_excl ); % removes the specified electrodes

tbl_channels = struct2table(EEG.chanlocs); % write the electrodes into a table. We will need this for later. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ZapLine


%% Plot Frequency Spectrum without ZapLine
figure; pop_spectopo(EEG, 1, [0  size(EEG.data,2)], 'EEG' , 'percent', 15, 'freq', [6 10 22], 'freqrange',[2 70],'electrodes','off');

close all


%% test different number of components:

line_noise_freq = 50; % 50Hz line noise
FLINE=line_noise_freq/EEG.srate; % line frequency
p.nfft = 1024;

% 3 components
EEGzap3 = EEG;
NREMOVE=3; % number of components to remove (our empirical tests showed that 7 components are required to sufficiently remove the line noise), in this subject 3 are definitely not enough. 
[clean3, noise3] = nt_zapline(EEG.data',FLINE,NREMOVE);
EEGzap3.data = clean3';
figure; pop_spectopo(EEGzap3, 1, [0  size(EEGzap3.data,2)], 'EEG' , 'percent', 15, 'freq', [6 10 22 40], 'freqrange',[2 70],'electrodes','off');

% 5 components
EEGzap5 = EEG;
NREMOVE=5; % number of components to remove (our empirical tests showed that 7 components are required to sufficiently remove the line noise), in this subject 5 are definitely not enough. 
[clean5, noise5] = nt_zapline(EEG.data',FLINE,NREMOVE);
EEGzap5.data = clean5';
figure; pop_spectopo(EEGzap5, 1, [0  size(EEGzap5.data,2)], 'EEG' , 'percent', 15, 'freq', [6 10 22 40], 'freqrange',[2 70],'electrodes','off');

% 7 components
EEGzap7 = EEG;
NREMOVE=7; % number of components to remove (our empirical tests showed that 7 components are required to sufficiently remove the line noise), in this subject 5 are definitely not enough. 
[clean7, noise7] = nt_zapline(EEG.data',FLINE,NREMOVE);
EEGzap7.data = clean7';
figure; pop_spectopo(EEGzap7, 1, [0  size(EEGzap7.data,2)], 'EEG' , 'percent', 15, 'freq', [6 10 22 40], 'freqrange',[2 70],'electrodes','off');


%% to plot the the raw, clean and removed signal (here for 3 components):

plot_zapline(clean3,EEG.data,FLINE) % takes a while

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% bad channel indentification (using cleanraw_data 1Hz temporary filteron a copy of the data)

% create the copy
EEG_badchan = EEG;

% Temporary highpass filter 1Hz (clean_artifacts works better with 1Hz)
EEG_badchan = pop_eegfiltnew(EEG_badchan,1,0);

% check for bad channels
% EEG_badchan = clean_artifacts(EEG_badchan, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );

% without using Window's criterion
EEG_badchan = clean_artifacts(EEG_badchan, 'FlatlineCriterion',5,'ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off','WindowCriterion','off' );

% visualize and compare results (new first, old second)
vis_artifacts(EEG_badchan, EEG)


% find index of the retained channels

% get index of retained channels for the actual data
[ind_retained,ind_bad] = get_el_index_clean_artifacts(EEG_badchan,EEG) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Filter the data

%% notch filter (should be replaced by Zapline or CleanLine)
% EEG = pop_eegfiltnew(EEG, 48,52,826,1,[],0); 

% we are using the results from ZapLine here
EEG =  EEGzap7;

% high-pass filter with default filter parameters from EEGLab:
EEG_002hz = pop_eegfiltnew(EEG,0.02,0); % high-pass filter passband edge 0.02Hz (cut-off 0.01 hz) no low-pass
EEG_02hz = pop_eegfiltnew(EEG,0.2,0); % high-pass filter passband edge 0.2Hz (cut-off 0.1 hz) no low-pass
EEG_05hz = pop_eegfiltnew(EEG,0.5,0); % high-pass filter passband edge 0.5Hz (cut-off 0.25 hz) no low-pass
EEG_1hz = pop_eegfiltnew(EEG,1,0); % high-pass filter passband edge 0.2Hz (cut-off 0.5 hz) no low-pass
EEG_2hz = pop_eegfiltnew(EEG,2,0); % high-pass filter passband edge 2hz (cut-off 1 hz) no low-pass
EEG_3hz = pop_eegfiltnew(EEG,3,0); % high_pass filter passband edge 3Hz (cut-off 2 hz) no low-pass
 
  
% low-pass filterwith default filter parameters from EEGLab:  
EEG_02hz_39hz = pop_eegfiltnew(EEG_02hz,0,39); % add low_pass 39Hz (cut-off 43 hz)  
EEG_3hz_39hz = pop_eegfiltnew(EEG_3hz,0,39); % add low_pass 39Hz (cut-off 43 hz)
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% select only retained electrodes for ICA (exclude bad electrodes from ICA decomposition)

EEG_002hz = pop_select( EEG_002hz,'channel',ind_retained );
EEG_02hz = pop_select( EEG_02hz,'channel',ind_retained );
EEG_05hz = pop_select( EEG_05hz,'channel',ind_retained );
EEG_1hz = pop_select( EEG_1hz,'channel',ind_retained );
EEG_2hz = pop_select( EEG_2hz,'channel',ind_retained );
EEG_3hz = pop_select( EEG_3hz,'channel',ind_retained );

EEG_02hz_39hz = pop_select( EEG_02hz_39hz,'channel',ind_retained );
EEG_3hz_39hz = pop_select( EEG_3hz_39hz,'channel',ind_retained );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ICA (computational time is high. thus only executed for 0.1 Hz, 0.5 Hz, 1Hz, 2Hz, 0.1Hz+43low-pass,2Hz+43low-pass )
EEG_ica_02hz = pop_runica(EEG_02hz, 'icatype', 'runica', 'extended',1,'interrupt','on','reorder','on');
EEG_ica_1hz = pop_runica(EEG_1hz, 'icatype', 'runica', 'extended',1,'interrupt','on','reorder','on');
EEG_ica_2hz = pop_runica(EEG_2hz, 'icatype', 'runica', 'extended',1,'interrupt','on','reorder','on');
EEG_ica_3hz = pop_runica(EEG_3hz, 'icatype', 'runica', 'extended',1,'interrupt','on','reorder','on');
EEG_ica_02hz_39hz = pop_runica(EEG_02hz_39hz, 'icatype', 'runica', 'extended',1,'interrupt','on','reorder','on');
EEG_ica_3hz_39hz = pop_runica(EEG_3hz_39hz, 'icatype', 'runica', 'extended',1,'interrupt','on','reorder','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% keep all ICA components (in the .etc) - temporary filtering
% for backprojection of ICA weight matrix and ICA sphere matrix

EEG_ica_3hz_39hz  = keepICA(EEG_ica_3hz_39hz)

%% backprojection of ICA weight matrix and ICA sphere matrix
EEG_ica_bpr_002hz = ica_foreign_backproject(EEG_002hz,EEG_ica_3hz_39hz);
EEG_ica_bpr_02hz = ica_foreign_backproject(EEG_02hz,EEG_ica_3hz_39hz);
EEG_ica_bpr_05hz = ica_foreign_backproject(EEG_05hz,EEG_ica_3hz_39hz);
EEG_ica_bpr_1hz = ica_foreign_backproject(EEG_1hz,EEG_ica_3hz_39hz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  PLOT the topoplots before automated ICA classifiers (ICLabel)

pop_topoplot(EEG_ica_02hz, 0, [1:20] ,'ica 0.1hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
pop_topoplot(EEG_ica_1hz, 0, [1:20] ,'ica 0.5hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
pop_topoplot(EEG_ica_2hz, 0, [1:20] ,'ica 1hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
pop_topoplot(EEG_ica_3hz, 0, [1:20] ,'ica 2hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');

pop_topoplot(EEG_ica_02hz_39hz, 0, [1:20] ,'ica 0.1hz 43hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
pop_topoplot(EEG_ica_3hz_39hz, 0, [1:20] ,'ica 1hz 43hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');

pop_topoplot(EEG_ica_bpr_002hz, 0, [1:20] ,'ica 0.01hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
pop_topoplot(EEG_ica_bpr_02hz, 0, [1:20] ,'ica 0.1hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
pop_topoplot(EEG_ica_bpr_05hz, 0, [1:20] ,'ica 0.25hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
pop_topoplot(EEG_ica_bpr_1hz, 0, [1:20] ,'ica 0.5hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ICLabel

%% Settings (probabilty of label to be excluded)
% { 'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other' };
ic_label_crit = [NaN NaN;0.8 1;0.8 1;0.8 1;0.8 1;0.8 1;0.8 1]


EEG_ica_02hz_ICL = iclabel(EEG_ica_02hz)
pop_viewprops(EEG_ica_02hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_02hz_ICL= pop_icflag(EEG_ica_02hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_02hz = sum(EEG_ica_02hz_ICL.reject.gcompreject)
EEG_ica_02hz_ICL = pop_subcomp( EEG_ica_02hz_ICL,find(EEG_ica_02hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability


EEG_ica_1hz_ICL = iclabel(EEG_ica_1hz)
pop_viewprops(EEG_ica_1hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_1hz_ICL= pop_icflag(EEG_ica_1hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_1hz = sum(EEG_ica_1hz_ICL.reject.gcompreject)
EEG_ica_1hz_ICL = pop_subcomp( EEG_ica_1hz_ICL,find(EEG_ica_1hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability


EEG_ica_2hz_ICL = iclabel(EEG_ica_2hz)
pop_viewprops(EEG_ica_2hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_2hz_ICL= pop_icflag(EEG_ica_2hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_2hz = sum(EEG_ica_2hz_ICL.reject.gcompreject)
EEG_ica_2hz_ICL = pop_subcomp( EEG_ica_2hz_ICL,find(EEG_ica_2hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability



EEG_ica_3hz_ICL = iclabel(EEG_ica_3hz)
pop_viewprops(EEG_ica_3hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_3hz_ICL = pop_icflag(EEG_ica_3hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_3hz = sum(EEG_ica_3hz_ICL.reject.gcompreject)
EEG_ica_3hz_ICL = pop_subcomp( EEG_ica_3hz_ICL,find(EEG_ica_3hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability



EEG_ica_02hz_39hz_ICL = iclabel(EEG_ica_02hz_39hz)
pop_viewprops(EEG_ica_02hz_39hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_02hz_39hz_ICL = pop_icflag(EEG_ica_02hz_39hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_02hz_39hz = sum(EEG_ica_02hz_39hz_ICL.reject.gcompreject)
EEG_ica_02hz_39hz_ICL = pop_subcomp( EEG_ica_02hz_39hz_ICL,find(EEG_ica_02hz_39hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability


EEG_ica_3hz_39hz_ICL = iclabel(EEG_ica_3hz_39hz)
EEG_ica_3hz_39hz_ICL = pop_icflag(EEG_ica_3hz_39hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_3hz_39hz = sum(EEG_ica_3hz_39hz_ICL.reject.gcompreject)
EEG_ica_3hz_39hz_ICL = pop_subcomp( EEG_ica_3hz_39hz_ICL,find(EEG_ica_3hz_39hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability

% version with temporary filtering:

EEG_ica_bpr_002hz_ICL = iclabel(EEG_ica_bpr_002hz)
pop_viewprops(EEG_ica_bpr_002hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_bpr_002hz_ICL= pop_icflag(EEG_ica_bpr_002hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_002hz = sum(EEG_ica_bpr_002hz_ICL.reject.gcompreject)
EEG_ica_bpr_002hz_ICL = pop_subcomp( EEG_ica_bpr_002hz_ICL,find(EEG_ica_bpr_002hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability


EEG_ica_bpr_02hz_ICL = iclabel(EEG_ica_bpr_02hz)
pop_viewprops(EEG_ica_bpr_02hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_bpr_02hz_ICL= pop_icflag(EEG_ica_bpr_02hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_02hz = sum(EEG_ica_bpr_02hz_ICL.reject.gcompreject)
EEG_ica_bpr_02hz_ICL = pop_subcomp( EEG_ica_bpr_02hz_ICL,find(EEG_ica_bpr_02hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability

EEG_ica_bpr_05hz_ICL = iclabel(EEG_ica_bpr_05hz)
pop_viewprops(EEG_ica_bpr_05hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_bpr_05hz_ICL= pop_icflag(EEG_ica_bpr_05hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_05hz = sum(EEG_ica_bpr_05hz_ICL.reject.gcompreject)
EEG_ica_bpr_05hz_ICL = pop_subcomp( EEG_ica_bpr_05hz_ICL,find(EEG_ica_bpr_05hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability

EEG_ica_bpr_1hz_ICL = iclabel(EEG_ica_bpr_1hz)
pop_viewprops(EEG_ica_bpr_1hz_ICL,0,[1:20]);set(gcf,'color','w');
EEG_ica_bpr_1hz_ICL= pop_icflag(EEG_ica_bpr_1hz_ICL, ic_label_crit); % to identify all comoponents which are NOT brain with at leas 80% probability
nr_bad_comp_ICL_1hz = sum(EEG_ica_bpr_1hz_ICL.reject.gcompreject)
EEG_ica_bpr_1hz_ICL = pop_subcomp( EEG_ica_bpr_1hz_ICL,find(EEG_ica_bpr_1hz_ICL.reject.gcompreject) , 0);
%EEG = pop_icflag(EEG, [0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN;NaN NaN]); % to identify components which are brain with at leas 80% probability


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOT the retained Indepedent Components (with labels from IC Label)

EEG_ica_02hz_ICL_post = iclabel(EEG_ica_02hz_ICL)
pop_viewprops(EEG_ica_02hz_ICL_post,0,[1:20]);set(gcf,'color','w');

EEG_ica_1hz_ICL_post = iclabel(EEG_ica_1hz_ICL)
pop_viewprops(EEG_ica_1hz_ICL_post,0,[1:20]);set(gcf,'color','w');

EEG_ica_2hz_ICL_post = iclabel(EEG_ica_2hz_ICL)
pop_viewprops(EEG_ica_2hz_ICL_post,0,[1:20]);set(gcf,'color','w');

EEG_ica_3hz_ICL_post = iclabel(EEG_ica_3hz_ICL)
pop_viewprops(EEG_ica_3hz_ICL_post,0,[1:20]);set(gcf,'color','w');


EEG_ica_02hz_39hz_ICL_post = iclabel(EEG_ica_02hz_39hz_ICL)
pop_viewprops(EEG_ica_02hz_39hz_ICL_post,0,[1:20]);set(gcf,'color','w');

EEG_ica_3hz_39hz_ICL_post = iclabel(EEG_ica_3hz_39hz_ICL)
pop_viewprops(EEG_ica_3hz_39hz_ICL_post,0,[1:20]);set(gcf,'color','w');


% temporary filtered files

EEG_ica_bpr_002hz_ICL_post = iclabel(EEG_ica_bpr_002hz_ICL)
pop_viewprops(EEG_ica_bpr_002hz_ICL_post,0,[1:20]);set(gcf,'color','w');


EEG_ica_bpr_02hz_ICL_post = iclabel(EEG_ica_bpr_02hz_ICL)
pop_viewprops(EEG_ica_bpr_02hz_ICL_post,0,[1:20]);set(gcf,'color','w');

EEG_ica_bpr_05hz_ICL_post = iclabel(EEG_ica_bpr_05hz_ICL)
pop_viewprops(EEG_ica_bpr_05hz_ICL_post,0,[1:20]);set(gcf,'color','w');

EEG_ica_bpr_1hz_ICL_post = iclabel(EEG_ica_bpr_1hz_ICL)
pop_viewprops(EEG_ica_bpr_1hz_ICL_post,0,[1:20]);set(gcf,'color','w');


%% sanity check for the first 10 independent components (plot frequency spectrum, time course, topography, and the probability of all labels)
 whichComponentToPlot = [1:10] % define the number components to plot
% whichComponentToPlot = 1
pop_prop_extended(EEG_ica_3hz_39hz_ICL, 0, whichComponentToPlot, NaN)
 
  
%% An alternative way to plot the retained Indepedent Components
% pop_topoplot(EEG_ica_02hz_ICL, 0, [1:20],'after IC Label 0.1hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
% pop_topoplot(EEG_ica_1hz_ICL, 0, [1:20] ,'after IC Label 0.5hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
% pop_topoplot(EEG_ica_2hz_ICL, 0, [1:20] ,'after IC Label 1hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
% pop_topoplot(EEG_ica_3hz_ICL, 0, [1:20] ,'after IC Label 2hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
% 
% pop_topoplot(EEG_ica_02hz_39hz_ICL, 0, [1:20] ,'after IC Label 0.11hz 43hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');
% pop_topoplot(EEG_ica_2hz_39hz_ICL, 0, [1:20] ,'after IC Label 1hz 43hz',[4 5] ,0,'electrodes','off'); set(gcf,'color','w');


 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MARA (another automated classifier)

% 
% EEG_ica_2hz = pre_mara_elec_rename_wakeman_henson(EEG_ica_2hz)
% EEG = EEG_ica_2hz
% eeglab redraw 
% pop_processMARA ( ALLEEG,EEG,CURRENTSET )
% nr_bad_comp_2hz = size(EEG_ica.data,1)-size(EEG.icaact,1)
% %EEG = pop_subcomp( EEG, [1    4    5    9   13   14   16   17   19   23   24   26   28   32   33   35   36   38   43   46   48   49   50   57   60   61   62   64   65   66   68   70   71   72   74   75   77   78   79   82   84   85   86   87   88   89   91   92   93   94   96   97   98  100  101  102  103  104  106  107  108  110  111  112  113  114  115  117  118  119  120  121  122  123], 0);
% %nr_bad_comp_2hz = length([1    4    5    9   13   14   16   17   19   23   24   26   28   32   33   35   36   38   43   46   48   49   50   57   60   61   62   64   65   66   68   70   71   72   74   75   77   78   79   82   84   85   86   87   88   89   91   92   93   94   96   97   98  100  101  102  103  104  106  107  108  110  111  112  113  114  115  117  118  119  120  121  122  123])
% EEG_ica_mara_2hz = EEG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% put bad electrodes back in data and chanlocs

all_files = {
'EEG_ica_02hz_ICL',...
'EEG_ica_1hz_ICL',...
'EEG_ica_2hz_ICL',...
'EEG_ica_3hz_ICL',...
'EEG_ica_02hz_39hz_ICL',...
'EEG_ica_3hz_39hz_ICL',...
'EEG_ica_bpr_002hz_ICL',...
'EEG_ica_bpr_02hz_ICL',...
'EEG_ica_bpr_05hz_ICL',...
'EEG_ica_bpr_1hz_ICL',...
}


tbl_channels_orig = tbl_channels;

% loops through all preprocessing version (different high-pass filters)
for i = 1:size(all_files,2)
    clear EEG EEG.chanlocs
    
    eval(['EEG = ',all_files{i},';']) 


%% but back electrodes into the EEG structure
EEG.chanlocs = table2struct(tbl_channels_orig)';
EEG.nbchan = size(EEG.chanlocs,2);
tmp_nan = nan(size(tbl_channels_orig,1),size(EEG.data,2));
tmp_data = EEG.data;
EEG.data = tmp_nan;
EEG.data(sort(ind_retained),:) = tmp_data;
%eeglab redraw 


%% interpolate bad components:
EEG = eeg_interp(EEG ,ind_bad,'spherical')

%% save variable with ending _preproc
eval([all_files{i},' = EEG;'])


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract ERPs (segment and baseline correction)

% loops through all preprocessing version (different high-pass filters)
for i = 1:size(all_files,2)

    clear EEG EEG.chanlocs EEG_famous EEG_unfamiliar EEG_scrambled ERP_famous ERP_unfamiliar ERP_scrambled
    
    eval(['EEG = ',all_files{i},';']) 
        
%% Epoch data
EEG_famous = pop_epoch(EEG, {'famous_new' 'famous_second_early' 'famous_second_late'}, [-0.2, 0.8]);
EEG_unfamiliar = pop_epoch(EEG, {'unfamiliar_new' 'unfamiliar_second_early' 'unfamiliar_second_late'}, [-0.2, 0.8]);
EEG_scrambled = pop_epoch(EEG, {'scrambled_new' 'scrambled_second_early' 'scrambled_second_late'}, [-0.2, 0.8]);


%% baseline correction 200ms
EEG_famous_bc200 = pop_rmbase(EEG_famous,[-200 0]);
EEG_unfamiliar_bc200 = pop_rmbase(EEG_unfamiliar,[-200 0]);
EEG_scrambled_bc200 = pop_rmbase(EEG_scrambled,[-200 0]);

%% baseline correction 100ms 
EEG_famous_bc100 = pop_rmbase(EEG_famous,[-100 0]);
EEG_unfamiliar_bc100 = pop_rmbase(EEG_unfamiliar,[-100 0]);
EEG_scrambled_bc100 = pop_rmbase(EEG_scrambled,[-100 0]);

%% compute ERP no baseline correction
ERP_famous_bc200 = mean(EEG_famous_bc200.data(61,:,:),3);  % index 61 = electrode E065
ERP_unfamiliar_bc200 = mean(EEG_unfamiliar_bc200.data(61,:,:),3);
ERP_scrambled_bc200 = mean(EEG_scrambled_bc200.data(61,:,:),3);

%% compute ERP 100ms baseline correction
ERP_famous_bc100 = mean(EEG_famous_bc100.data(61,:,:),3); % index 61 = electrode E065
ERP_unfamiliar_bc100 = mean(EEG_unfamiliar_bc100.data(61,:,:),3);
ERP_scrambled_bc100 = mean(EEG_scrambled_bc100.data(61,:,:),3);

%% compute ERP 200ms baseline correction
ERP_famous_nobc = mean(EEG_famous.data(61,:,:),3); % index 61 = electrode E065
ERP_unfamiliar_nobc = mean(EEG_unfamiliar.data(61,:,:),3);
ERP_scrambled_nobc = mean(EEG_scrambled.data(61,:,:),3);


eval(['ERP_famous_bc200_',all_files{i}(1:end-4) ' = ERP_famous_bc200;'])
eval(['ERP_unfamiliar_bc200_',all_files{i}(1:end-4) ' = ERP_unfamiliar_bc200;'])
eval(['ERP_scrambled_bc200_',all_files{i}(1:end-4) ' = ERP_scrambled_bc200;'])

eval(['ERP_famous_bc100_',all_files{i}(1:end-4) ' = ERP_famous_bc100;'])
eval(['ERP_unfamiliar_bc100_',all_files{i}(1:end-4) ' = ERP_unfamiliar_bc100;'])
eval(['ERP_scrambled_bc100_',all_files{i}(1:end-4) ' = ERP_scrambled_bc100;'])

eval(['ERP_famous_nobc_',all_files{i}(1:end-4) ' = ERP_famous_nobc;'])
eval(['ERP_unfamiliar_nobc_',all_files{i}(1:end-4) ' = ERP_unfamiliar_nobc;'])
eval(['ERP_scrambled_nobc_',all_files{i}(1:end-4) ' = ERP_scrambled_nobc;'])

eval(['EEG_famous_',all_files{i}(1:end-4) ' = EEG_famous;'])
eval(['EEG_unfamiliar_',all_files{i}(1:end-4) ' = EEG_unfamiliar;'])
eval(['EEG_scrambled_',all_files{i}(1:end-4) ' = EEG_scrambled;'])


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calcuate ERP Differences
%famous faces - unfamiliar faces
% famous faces - scrabled

all_files_famous = who('ERP_famous*')
all_files_unfamiliar = who('ERP_unfamiliar*')
all_files_scrambled = who('ERP_scrambled*')

for i =  1:length(all_files_famous)
eval(['Diff_fam_unfam',all_files_famous{i}(11:end), ' = ERP_famous',all_files_famous{i}(11:end),'- ERP_unfamiliar',all_files_unfamiliar{i}(15:end),';'])
eval(['Diff_fam_scram',all_files_famous{i}(11:end), ' = ERP_famous',all_files_famous{i}(11:end),'- ERP_scrambled',all_files_scrambled{i}(14:end),';'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOTTING: 


% course_folder =  '~/Dropbox/AA_transfer/OHBM_educational_course/final_folder/'
% cd(course_folder)
% load('single_subjects_preprocessed.mat')


%% Plot all different ERPs for different high-pass filters and ICLabel versions
% 
% here ERP for "famous faces" baseline correction 100ms

% vector for plotting (needs to be adjusted with different baseline length)
times = -196:4:800;

figure
hold on
plot(times,ERP_famous_bc100_EEG_ica_02hz,'r','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_1hz,'-.r','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_2hz,'--r','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_3hz,':r','lineWidth',2)

plot(times,ERP_famous_bc100_EEG_ica_3hz_39hz,'b','lineWidth',2)

plot(times,ERP_famous_bc100_EEG_ica_bpr_002hz,'--g','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_bpr_02hz,'g','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_bpr_05hz,':g','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_bpr_1hz,'-.g','lineWidth',2)


legend({'0.1Hz ICLabel','0.5Hz ICLabel','1Hz ICLabel','2Hz ICLabel','2Hz & 43 Hz ICLabel','temp 2Hz-43Hz ICLabel passed to 0.01 Hz','temp 2Hz-43Hz ICLabel passed to 0.1 Hz','temp 2Hz-43Hz ICLabel passed to 0.25 Hz','temp 1Hz-49Hz ICLabel passed to 1 Hz','temp 2Hz-39Hz ICLabel passed to 0.5 Hz', 'no ICA unfiltered'},'Location','northeastoutside')
set(gcf,'color','w');
title('ERP famous faces')



%% Plot all Differences between "famous faces" and "scrambled" different high-pass filters and ICLabel versions
% 
% here ERP for "famous faces"  is substracted by "scrambled images"


figure
hold on
plot(times,Diff_fam_scram_nobc_EEG_ica_02hz,'r','lineWidth',2)
plot(times,Diff_fam_scram_nobc_EEG_ica_1hz,'-.r','lineWidth',2)
plot(times,Diff_fam_scram_nobc_EEG_ica_2hz,'--r','lineWidth',2)
plot(times,Diff_fam_scram_nobc_EEG_ica_3hz,':r','lineWidth',2)

plot(times,Diff_fam_scram_nobc_EEG_ica_3hz_39hz,'b','lineWidth',2)

plot(times,Diff_fam_scram_nobc_EEG_ica_bpr_002hz,'--g','lineWidth',2)
plot(times,Diff_fam_scram_nobc_EEG_ica_bpr_02hz,'g','lineWidth',2)
plot(times,Diff_fam_scram_nobc_EEG_ica_bpr_05hz,':g','lineWidth',2)
plot(times,Diff_fam_scram_nobc_EEG_ica_bpr_1hz,'-.g','lineWidth',2)


legend({'0.1Hz ICLabel','0.5Hz ICLabel','1Hz ICLabel','2Hz ICLabel','2Hz & 43 Hz ICLabel','temp 2Hz-43Hz ICLabel passed to 0.01 Hz','temp 2Hz-43Hz ICLabel passed to 0.1 Hz','temp 2Hz-43Hz ICLabel passed to 0.25 Hz','temp 1Hz-49Hz ICLabel passed to 1 Hz','temp 2Hz-39Hz ICLabel passed to 0.5 Hz', 'no ICA unfiltered'},'Location','northeastoutside')
set(gcf,'color','w');
title('DIFF famous faces - scrambled')



%% Effect of different baseline corrections:


%% Effect of Baseline Correction and High-Pass Filter interaction: 

% dashed line = 200ms baseline
% solid line = 100ms baseline
% almost no difference if 200 or 100ms baseline

%% Grand Average "famous faces"
figure
hold on
plot(times,ERP_famous_nobc_EEG_ica_bpr_002hz,'r','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_bpr_002hz,':r','lineWidth',2)
plot(times,ERP_famous_bc200_EEG_ica_bpr_002hz,'--r','lineWidth',2)
plot(times,ERP_famous_nobc_EEG_ica_bpr_02hz,'Color',[0 255 255]/255,'lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_bpr_02hz,'Color',[0 255 255]/255,'lineWidth',2,'LineStyle',':')
plot(times,ERP_famous_bc200_EEG_ica_bpr_02hz,'Color',[0 255 255]/255,'lineWidth',2,'LineStyle','--')
plot(times,ERP_famous_bc100_EEG_ica_3hz,'b','lineWidth',2)
plot(times,ERP_famous_bc100_EEG_ica_3hz,':b','lineWidth',2)
plot(times,ERP_famous_bc200_EEG_ica_3hz,'--b','lineWidth',2)
legend({'no bc 0.01 hz', '100ms bc 0.01 hz','200ms bc 0.01 hz','no bc 0.1 hz', '100ms bc 0.1 hz','200ms bc 0.1 hz','no bc 2 hz', '100ms bc 2 hz','200ms bc 2 hz'})
set(gcf,'color','w');
title('Interaction between baseline corrections and high-pass filtering (Grand Average)')



%% ERP Differences "famous faces" - "scrambled images"

figure
hold on
plot(times,Diff_fam_scram_nobc_EEG_ica_bpr_002hz,'r','lineWidth',2)
plot(times,Diff_fam_scram_bc100_EEG_ica_bpr_002hz,':r','lineWidth',2)
plot(times,Diff_fam_scram_bc200_EEG_ica_bpr_002hz,'--r','lineWidth',2)
plot(times,Diff_fam_scram_nobc_EEG_ica_bpr_02hz,'Color',[0 255 255]/255,'lineWidth',2)
plot(times,Diff_fam_scram_bc100_EEG_ica_bpr_02hz,'Color',[0 255 255]/255,'lineWidth',2,'LineStyle',':')
plot(times,Diff_fam_scram_bc200_EEG_ica_bpr_02hz,'Color',[0 255 255]/255,'lineWidth',2,'LineStyle','--')
plot(times,Diff_fam_scram_bc100_EEG_ica_3hz,'b','lineWidth',2)
plot(times,Diff_fam_scram_bc100_EEG_ica_3hz,':b','lineWidth',2)
plot(times,Diff_fam_scram_bc200_EEG_ica_3hz,'--b','lineWidth',2)
legend({'no bc 0.01 hz', '100ms bc 0.01 hz','200ms bc 0.01 hz','no bc 0.1 hz', '100ms bc 0.1 hz','200ms bc 0.1 hz','no bc 2 hz', '100ms bc 2 hz','200ms bc 2 hz'})
set(gcf,'color','w');
title('Interaction between baseline corrections and high-pass filtering (ERP differences)')







% 
% 
% 
% %% PLOT all different preprocessings
% 
% all_files_famous = who('ERP_famous_bc100_EEG_*')
% 
% % vector for plotting (needs to be adjusted with different baseline length)
% times = -196:4:800;
% 
% % for color grading
% startColor = [255];
% scaling_f = 255/length(all_files_famous)-0.3
% clear legend_names
% figure
% for ii = 1:length(all_files_famous)
%     plot(times, eval(all_files_famous{ii}), 'Color', [startColor/startColor,(startColor - (ii*scaling_f))/startColor,0] );
%     hold on;
%     legend_names{ii} = all_files_famous{ii}(18:end);
% end
% legend(legend_names)
% 
