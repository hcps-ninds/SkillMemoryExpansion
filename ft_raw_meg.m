function ft_raw_meg(subjectID)
addpath(['./fieldtrip_latest/']);
ft_defaults; 

% subjectID = 'BRWDDIXO';

meg_dir = dir(['../../MEGdata/' subjectID '_Day1_targ_*.ds']);

% Loading MEG Data
cfg = [];
cfg.dataset = [meg_dir.folder '/' meg_dir.name];
cfg.continuous = 'yes';
meg_data = ft_preprocessing(cfg);

cfg = [];
cfg.continuous = 'yes';
cfg.channel = 'UPPT001';
meg_trigger = ft_preprocessing(cfg, meg_data);

% Applying third-order gradient
cfg = [];
cfg.gradient = 'G3BR';
meg_data = ft_denoise_synthetic(cfg, meg_data);

if ismember(subjectID, {'FLDSUUNE', 'HDYJKAXR', 'KLORMPJU', 'OVMYOBTX'})
    cfg = [];
    meg_label = meg_data.hdr.label;
    rm_idx = ismember(meg_label, 'MLT16');
    meg_label(rm_idx) = []; 
    cfg.channel = meg_label;
    meg_data = ft_selectdata(cfg, meg_data);
end

% High-pass and notch filtering
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = .1;
cfg.hpfiltord = 4;
cfg.dftfilter = 'yes';
cfg.dftreplace = 'neighbour';
cfg.dftfreq = [60 120];
cfg.channel = 'meg';
meg_data = ft_preprocessing(cfg, meg_data);

save(['./raw_meg/' subjectID '.mat'], 'meg_data', 'meg_trigger', '-v7.3')