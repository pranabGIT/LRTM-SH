%{
Details of Expt-A can be found here:

https://docs.google.com/document/d/1kU337S6c0JNcKNqAqmFwEGcnsaxZ41Rgk9MpLqnXuIs/edit

%}

clc; clear all;
addpath /home/pranab/Documents/MATLAB-AUX/m_map/

%%%%% YEARS/ El Nino or La Nina??
yr = 1980:2004; % Dec - centred == 1980 means DJF/D1980/J81/F81

% ElNino
yr1 = [1982, 1987, 1991, 1997, 2002];
    
% list for La Nina
%yr1 = [1988, 1995, 1998, 1999, 2000]; 

%%
% calculate NCEP DJF ENSO composite

load NcepZG_sig_mat_DJF_1980_2004SH_NoFilt.mat
rst = 'NCEP';

za = sig_mat; [mt, nt, ot] = size(za);
% clear lon lat; % to be used to interpolate stpavg data to this grid

zaDJF = zeros(mt/90, nt, ot);

k = 1;
for i = 1:90:length(za)
    zaDJF(k,:,:) = nanmean(za(i:i+89,:,:),1);
    k = k + 1;
end

k2 = 1;
for j5 = 1:length(yr1)
    ensInd(k2) = find (yr == yr1(j5));
    k2 = k2+1;
end

ghano_ENSO = squeeze(nanmean(zaDJF(ensInd,:,:),1));
ghano_ENSO = ghano_ENSO';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKING LON 0==>360 /// This step is needed because sig_mat of the models
% are on 0-->360 RANGE
lon(lon<0) = 180 + lon(lon>0);
pln = find(lon>180); 
ln1 = lon(pln);
pln2 = find(lon<=180);
ln2 = lon(pln2);
lon = [ln2 ln1];
ghln = ghano_ENSO(pln,:); ghln2 = ghano_ENSO(pln2,:);
ghano_ENSO360 = cat(1, ghln2, ghln);

    % TARGET GRID %
latR = lat; lonR = lon;
[Xq360,Yq360] = meshgrid(latR, lonR); % TARGET GRID for interpolation

% for correlation with stp ==> LRTM [0 --> +360 lon range]
lt = find (latR<=-20);
aera360 = ghano_ENSO360(:, lt);

%% Plot Reanalysis ENSO COMPOSITE
figure(1)
subplot(2,2,1)

m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonR,latR,ghano_ENSO360',[-200,-100:10:100,200], 'LineColor','none')
m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
title (['EL Nino Comp - ', rst])

% fnm = ['ERA_ElNino_SH.png'];
% print ('-r300', fnm, '-dpng')

%% loading step average for models

stp_ens = zeros (nt, ot);
ghENSO_ens = zeros (nt, ot);
% idc = input ('Which model is it?? 1 for CCSM4 :: 2 for IPSL :: 3 for MIROC5 :: 4 for HadGEM2A :: 5 for GFDLCM3 :: 6 for MPIesmMR :: 7 for MRI-CGCM3 :: 8 for ACCESS1-3 :: 9 for NorESM1 ::');
k = 0;
for idc = 1:9
    % STP AVG COMPOSITE
    if idc == 1
            modl = 'CCSM4';
        elseif idc == 2
            modl = 'IPSLcm5aMR';
        elseif idc == 3
            modl = 'MIROC5';
        elseif idc == 4
            modl = 'HadGEM2A';
        elseif idc == 5
            modl = 'GFDLCM3';
        elseif idc == 6
            modl = 'MPIesmMR';
        elseif idc == 7
            modl = 'MRI-CGCM3';
        elseif idc == 8
            modl = 'ACCESS1-3';
        elseif idc == 9
            modl = 'NorESM1';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Loading step response to ENSO for the model ///// has 0 --> +360 lon range
    stpsv = ['stp_',modl, '_DJF_1980_2004_ExptA.mat'];
    load (stpsv)
    [X,Y] = meshgrid(lat, lon);
 
    a1 = interp2(X,Y,stpavg',Xq360,Yq360); % interpolation from X,Y to Xq,Yq
    
    stp_ens = stp_ens + a1';
   
    lt = find (latR<=-20);
    bstp = a1(:,lt);
    % checking correlation of individual models
    ['corr for ',modl, ' LRTM v NCEP']
    
    nancorr2(bstp,aera360)
    
    
%#####################################################################################################
%#####################################################################################################

    % MODEL COMPOSITE ///// has 0 --> +360 lon range
    % loading sig_mat for the model to get ENSO composite for the model
    fl = ['sig_mat_DJF_1980_2004_',modl,'_amip.mat'];
    load (fl)
    clear lon lat
   
    [mc, nc, oc] = size(sig_mat);
    ghENSO_mod = zeros (nc, oc);

    for i = ensInd
        a = nanmean(sig_mat((i-1)*90+1:(i-1)*90+1+90, : ,:),1);
        % Interpolate from model grid to target NorESM1 grid (taken as the target grid for interpolation)
        
        ghENSO_mod = ghENSO_mod+squeeze(a);
    end
    
    ghENSO_mod1 = ghENSO_mod/length(ensInd);
    ghENSO_mod2 = squeeze(ghENSO_mod1);

    a1 = interp2(X,Y,ghENSO_mod2',Xq360,Yq360); % interpolation from X,Y to Xq,Yq
    ghENSO_ens = ghENSO_ens + a1';
  
    % for individual model composite vs ERA correlation
    cmod = a1; cmod = cmod(:, lt); 
    ['corr for ',modl, ' ENSO comp v NCEP']
    nancorr2(cmod,aera360)
    
    k = k+1;
end

stp_ens = stp_ens/k;
ghENSO_ens = ghENSO_ens/k;

['shape of Model ensemble stpavg :: ']
size(stp_ens)


['shape of Model ensemble stpavg :: ']
size(ghENSO_ens)


%%
% Plot STP ENS
figure(1)
subplot(2,2,2)
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonR,latR,stp_ens,[-200,-100:10:100,200], 'LineColor','none')
m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
title ('step avg. ensemble - ENSO')

%% Plot MODEL GH ENSO COMPOSITE
figure(1)
subplot(2,2,3)
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonR,latR,ghENSO_ens,[-200,-100:10:100,200], 'LineColor','none')
% m_contourf(lonm,latm,ghENSO_ens,[-200,-100:5:100,200], 'LineColor','none')

m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
title ('model ghano ensemble - ENSO')

% 
% This following one is to produce the colormap and delete the spatial plt
% manually
subplot(2,2,4)
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonR,latR,ghENSO_ens,[-200,-100:10:100,200], 'LineColor','none')
% m_contourf(lonm,latm,ghENSO_ens,[-200,-100:5:100,200], 'LineColor','none')

m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
colorbar

%% Spatial correlation between stpavg and ghano_ENSO(ERA)
era = ghano_ENSO';

% Removing NaN values in these metrices
% as = find(isnan(stp_ens));
% stp_ens(as) = []; ghENSO_ens(as) = []; era(as) = [];

lt = find(latR<=-20);
stp_ens = stp_ens(lt,:);
era = era (lt,:);
ghENSO_ens = ghENSO_ens(lt, :);

rspat1 = nancorr2(stp_ens,era)
rspat2 = nancorr2(ghENSO_ens, era)