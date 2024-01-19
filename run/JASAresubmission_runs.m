clear all
savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';
source = 1; % 1 - TIMIT wav files; 2 - synthesized wav files, 3 - synthesized with f0, 4 - synthesized F3 reflected, 5 - synthesized with f0 F3 reflected
trackBW = 1;
runVTRReplicate
save_error

% 
clear all
savedir = '..\results\EKS_WS_Praat_trackBW0_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';
source = 1; % 1 - TIMIT wav files; 2 - synthesized wav files, 3 - synthesized with f0, 4 - synthesized F3 reflected, 5 - synthesized with f0 F3 reflected
trackBW = 0;
runVTRReplicate
save_error

% 
clear all
savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';
source = 4; % 1 - TIMIT wav files; 2 - synthesized wav files, 3 - synthesized with f0, 4 - synthesized F3 reflected, 5 - synthesized with f0 F3 reflected
trackBW = 1;
runVTRReplicate
save_error_noclass

clear all
savedir = '..\results\EKS_WS_Praat_trackBW1_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';
source = 5; % 1 - TIMIT wav files; 2 - synthesized wav files, 3 - synthesized with f0, 4 - synthesized F3 reflected, 5 - synthesized with f0 F3 reflected
trackBW = 1;
runVTRReplicate
save_error

clear all
savedir = '..\results\EKS_WS_Praat_trackBW0_ARMAcep15_ar12ma0_win20_pe07_fs7000_timitVAD_useCorr1_Jacoblin_Q5e4bw1e4_BWfix';
source = 5; % 1 - TIMIT wav files; 2 - synthesized wav files, 3 - synthesized with f0, 4 - synthesized F3 reflected, 5 - synthesized with f0 F3 reflected
trackBW = 0;
runVTRReplicate
save_error

% make PDFs of plots
clear all
run_make_VTR_Results_All
