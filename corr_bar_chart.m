clear;clc;

modl = {'CCSM4', 'IPSLcm5aMR', 'MIROC5', 'HadGEM2A', 'GFDLCM3', 'MPIesmMR', 'MRI-CGCM3', 'ACCESS1-3', 'NorESM1', 'Ensemble'};

%%%%%%%%%%%%%%% EL NINO %%%%%%%%%%%%%%%%
load CORR_SH_ExptA_ElNino; a = rLRTM;
load CORR_SH_ExptB_ElNino; b = rLRTM;
load CORR_SH_ExptC_ElNino; c = rLRTM;
d = rmod;
figure(1)
subplot(2,1,1)
bar([a b c d])
ax=gca;
ax.FontSize = 8;
ylim([-0.1 0.9])
ylabel('Correlation')
xticklabels(modl)
xtickangle(45)
grid on
legend('Expt-A', 'Expt-B', 'Expt-C', 'Model ens', 'orientation', 'horizontal')
title ('El Nino')
%%%%%%%%%%%%%%% LA NINA %%%%%%%%%%%%%%%%
load CORR_SH_ExptA_LaNina; a = rLRTM;
load CORR_SH_ExptB_LaNina; b = rLRTM;
load CORR_SH_ExptC_LaNina; c = rLRTM;
d = rmod;
figure(1)
subplot(2,1,2)
bar([a b c d])
ax=gca;
ax.FontSize = 8;
ylim([-0.1 0.9])
ylabel('Correlation')
xticklabels(modl)
xtickangle(45)
grid on
legend('Expt-A', 'Expt-B', 'Expt-C', 'Model ens', 'orientation', 'horizontal')
title ('La Nina')
