
% Resolve paths relative to this script's location (works on any machine)
script_dir = fileparts(mfilename('fullpath'));
data_dir_scal     = fullfile(script_dir, 'benchmark-output', 'scaling');
data_dir_cond     = fullfile(script_dir, 'benchmark-output', 'conditioning');
data_dir_amg      = fullfile(script_dir, 'benchmark-output', 'amg');
data_dir_100_b100 = fullfile(script_dir, 'benchmark-output', 'amg_100x100_b100');
data_dir_100_b256 = fullfile(script_dir, 'benchmark-output', 'amg_100x100_b256');
data_dir_100_b512 = fullfile(script_dir, 'benchmark-output', 'amg_100x100_b512');
data_dir_real     = fullfile(script_dir, 'benchmark-output', 'real');

%% Scal benchmarks with k = 10^3, 10^6, 10^9
%plot_scal_results(data_dir_scal, '20260224_102848_scaling_results.csv', '20260224_102848_scaling_breakdown.csv', 'worst_ortho');
%plot_scal_results(data_dir_scal, '20260224_104655_scaling_results.csv', '20260224_104655_scaling_breakdown.csv', 'worst_ortho');
%plot_scal_results(data_dir_scal, '20260224_105319_scaling_results.csv', '20260224_105319_scaling_breakdown.csv', 'worst_ortho');

%% AMG 64x64 well-conditioned (κ ≈ 54), b_sz=64
%plot_gsvd_results(data_dir_amg, '20260227_114514_gsvd_results.csv', '20260227_114514_gsvd_breakdown.csv', 'worst_ortho');

%% AMG 64x64 ill-conditioned random V (κ ≈ 1.5e7), b_sz=64
%plot_gsvd_results(data_dir_amg, '20260227_114525_gsvd_results.csv', '20260227_114525_gsvd_breakdown.csv', 'worst_ortho');

%% AMG 100x100 well-conditioned (κ ≈ 84), block size sweep (skip_apps)
%plot_gsvd_results(data_dir_100_b100, '20260227_152336_gsvd_results.csv', '20260227_152336_gsvd_breakdown.csv', 'worst_ortho');
%plot_gsvd_results(data_dir_100_b256, '20260227_153049_gsvd_results.csv', '20260227_153049_gsvd_breakdown.csv', 'worst_ortho');
%plot_gsvd_results(data_dir_100_b512, '20260227_153355_gsvd_results.csv', '20260227_153355_gsvd_breakdown.csv', 'worst_ortho');

%% FEM_Problem_1 (m=73631, n=4421), 128 threads, skip_apps, block size sweep
%plot_gsvd_results(data_dir_real, '20260227_165034_gsvd_results.csv', '20260227_165034_gsvd_breakdown.csv', 'worst_ortho');  % b=256
%plot_gsvd_results(data_dir_real, '20260227_174043_gsvd_results.csv', '20260227_174043_gsvd_breakdown.csv', 'worst_ortho');  % b=512
%plot_gsvd_results(data_dir_real, '20260227_183047_gsvd_results.csv', '20260227_183047_gsvd_breakdown.csv', 'worst_ortho');  % b=1024

%% FEM_Problem_2 (m=75860, n=4812), 128 threads, block size sweep
plot_gsvd_results(data_dir_real, '20260305_140836_gsvd_results.csv', '20260305_140836_gsvd_breakdown.csv', 'worst_ortho');  % 128t, b=256
plot_gsvd_results(data_dir_real, '20260305_142015_gsvd_results.csv', '20260305_142015_gsvd_breakdown.csv', 'worst_ortho');  % 128t, b=512
plot_gsvd_results(data_dir_real, '20260305_143356_gsvd_results.csv', '20260305_143356_gsvd_breakdown.csv', 'worst_ortho');  % 128t, b=1024

%% FEM_Problem_2 (m=75860, n=4812), thread scaling (b=256)
plot_gsvd_results(data_dir_real, '20260305_160014_gsvd_results.csv', '20260305_160014_gsvd_breakdown.csv', 'worst_ortho');  % 8t, b=256
plot_gsvd_results(data_dir_real, '20260305_161745_gsvd_results.csv', '20260305_161745_gsvd_breakdown.csv', 'worst_ortho');  % 16t, b=256
plot_gsvd_results(data_dir_real, '20260305_162725_gsvd_results.csv', '20260305_162725_gsvd_breakdown.csv', 'worst_ortho');  % 32t, b=256
plot_gsvd_results(data_dir_real, '20260305_163520_gsvd_results.csv', '20260305_163520_gsvd_breakdown.csv', 'worst_ortho');  % 64t, b=256
