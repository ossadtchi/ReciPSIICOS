# ReciPSIICOS

Data source: https://drive.google.com/drive/folders/1kZjO7CgZNYVGtrciPKmdTiAvkHW4pUsb?usp=sharing

Simulations:

1.1 save_simulations_snr.m -- Generate simulations for different methods and different SNRs, the resulting activation maps saved at Z_total_full.mat, and source locations in picked_src.mat (available in GDrive).
	
1.2 plot_metrics_simulations_snr.m -- calculate point spreading and bias for precomputed Z_total_full and picked_src and draw them on a graph.

2.1 save_simulations_Gerror.m -- Generate simualtions for different forward model inaccuracies, using matrices C_re_005, C_re_01, C_re_02 (GDrive), save the results into Z_total_G (GDrive).
	
2.2 plot_metrics_simulations_Gerror.m -- calculate point spreading and bias for precomputed Z_total_G and picked_src and draw them on a graph.

3.1 save_simulations_rank.m -- Generate simulations for different methods and different SNRs, the resulting activation maps saved at Z_total_rank.mat, and source locations in picked_src_rank.mat (available in GDrive).
	
3.2 plot_simulations_rank.m -- calculate point spreading and bias for precomputed Z_total_rank and picked_src_rank and draw them on a graph.

4.1 plot_simulations_threesrc.m -- visualize simulations with three sources (input data on GDrive).

Real data:

5.1 MMN analysis: real_data_mmn.m -- analysis of auditory evoked MMN, using mmn.mat and G3_mmn.mat (GDrive).
5.2 40 Hz audio stim: ReciPSIICOS_real_data.m , today 09072020, will upload data by 11072020
