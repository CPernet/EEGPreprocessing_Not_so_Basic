 
function EEG = ica_foreign_backproject(EEG,EEG_temp);

% EEG = EEG file to backproject (e.g. EEG_02hz)
% EEG_temp = EEG file from which ICA components come from (e.g. EEG_ica_2hz_39hz)


    wts = EEG_temp.icaweights;
    sph = EEG_temp.icasphere;
    

% Remove any existing ICA solutions from your original dataset
    EEG.icaact      = [];
    EEG.icasphere   = [];
    EEG.icaweights  = [];
    EEG.icachansind = [];
    EEG.icawinv     = [];
    
    EEG.icasphere   = sph;
    EEG.icaweights  = wts;
    EEG.icachansind = EEG_temp.icachansind;
    EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv
    
end
