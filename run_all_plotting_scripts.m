
data_dir_scal = 'C:\Users\mymel\Documents\GitHub\QLCQRRT_benchmark\benchmark-output\scaling\';
data_dir_cond = 'C:\Users\mymel\Documents\GitHub\QLCQRRT_benchmark\benchmark-output\conditioning\';
data_dir_amg  = 'C:\Users\mymel\Documents\GitHub\QLCQRRT_benchmark\benchmark-output\amg\';

%% Scal benchmarks with k = 10^3, 10^6, 10^9
%plot_scal_results(data_dir_scal, '20260224_102848_scaling_results.csv',
%'20260224_102848_scaling_breakdown.csv', 'worst_ortho');
%plot_scal_results(data_dir_scal, '20260224_104655_scaling_results.csv',
%'20260224_104655_scaling_breakdown.csv', 'worst_ortho');
%plot_scal_results(data_dir_scal, '20260224_105319_scaling_results.csv',
%'20260224_105319_scaling_breakdown.csv', 'worst_ortho');

%% Case 1: Well-conditioned AMG base (κ ≈ 54)
plot_gsvd_results(data_dir_amg, '20260226_154140_gsvd_results.csv', '20260226_154140_gsvd_svals.csv', 'worst_ortho');

%% Case 2: Ill-conditioned random V (κ ≈ 1.5e7)
%plot_gsvd_results(data_dir_amg, '20260226_154156_gsvd_results.csv', '20260226_154156_gsvd_svals.csv', 'worst_ortho');