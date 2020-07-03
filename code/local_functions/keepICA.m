function EEG = keepICA(EEG)

    EEG.etc.beforeICremove.icaact = EEG.icaact;
    EEG.etc.beforeICremove.icawinv = EEG.icawinv;
    EEG.etc.beforeICremove.icasphere = EEG.icasphere;
    EEG.etc.beforeICremove.icaweights = EEG.icaweights;
    EEG.etc.beforeICremove.chanlocs = EEG.chanlocs;
    
end