%% Loading NWB data
matnwb_path = '/Users/mac/Projects/matnwb';
addpath(genpath(matnwb_path));
generateCore();

data_path = '/Volumes/rc/lab/D/DuvelleE/mvdmlab_ED_fs/Gillespie_et_al_2021_dataset/sub-despereaux_ses-despereaux-07_behavior+ecephys_draft.nwb';
nwb = nwbRead(data_path);

% Get LFP data
LFP = nwb.acquisition.values{1}.data;

% Get position data
pos = nwb.processing.get('behavior').nwbdatainterface.get('position').spatialseries.get('series_0').data;