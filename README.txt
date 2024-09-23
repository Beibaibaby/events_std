Run
"export QT_QPA_PLATFORM='offscreen' "
to aviod gui warning


julia run_6.jl --T 2000 --dir_name_in "/gpfs/data/doiron-lab/draco/results_corr/exp_test" --c_noise 0.01 --sigma_noise 1.0 --event_thre 0.8 --large_peak_mean 200 --peak_ratio 4000 --ie_sign false --ee_sign false