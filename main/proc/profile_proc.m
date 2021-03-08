octa_ffnames = dir('F:\img\DM_175003\OCTA\2020_08_12_OS\Raw\*.octa');
Gc = [];
all_times = cell(numel(octa_ffnames), 1);
for ii=1:numel(octa_ffnames)
    octa_ffname = fullfile(octa_ffnames(ii).folder, octa_ffnames(ii).name);
    [~, all_times{ii}, Gc] = lv_to_ml_proc_oct_gpu(octa_ffname, Gc, true);
end