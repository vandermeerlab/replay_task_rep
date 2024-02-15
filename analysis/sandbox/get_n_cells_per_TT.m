%%
data_paths = getGuptaDataPath([]);
for p_i = 1:length(data_paths)
    cd(data_paths{p_i});
    cfg_spikes = {};
    cfg_spikes.load_questionable_cells = 1;
    S = LoadSpikes(cfg_spikes);

    tt_nums = max(unique(S.usr.tt_num));
    cell_counts{p_i} = zeros(1, tt_nums);
    for i = 1:tt_nums
        cell_counts{p_i}(i) = sum(S.usr.tt_num == i);
    end
end