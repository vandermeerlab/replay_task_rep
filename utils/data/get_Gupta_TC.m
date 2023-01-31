%%
data_paths = getAdrDataPath([]);
for p_i = 1:length(data_paths)
    cd(data_paths{p_i});
    load_adrlab_data();
end