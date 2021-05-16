clc; clear all;

addpath /home/pranab/Documents/MATLAB-AUX/m_map/

%%%%% YEARS/ El Nino or La Nina??

yr = 1980:2004; % Dec - centred == 1980 means DJF/D1980/J81/F81

% ElNino
yr1 = [1982, 1987, 1991, 1997, 2002];
    
% list for La Nina
%yr1 = [1988, 1995, 1998, 1999, 2000]; 


%%%%% Using NorESM1 grid as the basic grid fr interpolation %%%%%%%%%%%%%%%%%%%%%%
load NorESM1_DJF_1980_2004; [mt, nt] = size(stpavg); latm = lat; lonm = lon; clear stpavg
[Xq,Yq] = meshgrid(lat, lon); clear lon lat; % to be used to interpolate stpavg data to this grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% load ERA5 DJF data
ncfile = 'ERA5_DJF_1980_2005_gh250.nc';
ncinfo(ncfile)
ncdisp(ncfile)
gh250 = ncread(ncfile,'z')* 0.23678 + 100537.0379; % SCALE FACTOR and OFFSET included

lon = ncread(ncfile,'longitude') ;
% MAKING LON 0==>360
lon(lon<0) = 360+lon(lon<0);
pln = find(lon>=180); pln2 = find(lon<180);
ln1 = lon(pln); ln2 = lon(pln2);
lon = [ln2;ln1];

ghln = gh250(pln,:,:); ghln2 = gh250(pln2,:,:);
gh250ln = cat(1, ghln, ghln2);

lat = ncread(ncfile,'latitude') ;

% Finding DJF clim and anomaly
DJF_clim = zeros(1440, 441, 3);
DJF_clim(:,:,1) = nanmean(gh250ln(:,:,1:3:end),3);
DJF_clim(:,:,2) = nanmean(gh250ln(:,:,2:3:end),3);
DJF_clim(:,:,3) = nanmean(gh250ln(:,:,3:3:end),3);
summClim = repmat(DJF_clim, [1,1,26]);
gh250ano = gh250ln - summClim;
[m,n,o] = size(gh250ano);

ghanoDJF = zeros(1440, 441, 25); % 26 ==> no.s of years
k = 1;
for j = 3:3:o-1 % starts from Dec-1980
    ghanoDJF(:,:,k) = nanmean(gh250ano(:,:,j:j+2),3);
    k = k+1;
end

k2 = 1;
for j5 = 1:length(yr1)
    ensInd(k2) = find (yr == yr1(j5));
    k2 = k2+1;
end

ghano_ElNino = nanmean(ghanoDJF(:,:,ensInd),3);

lone = lon; late = lat;
% Interpolate from ERA grid to target NorESM1 grid (taken as the target grid for interpolation)
[X,Y] = meshgrid(late, lone);
ghano_ENSO = interp2(X,Y,ghano_ElNino,Xq,Yq);


% for correlation
lt = find (latm<=-20);
aera = ghano_ENSO(:, lt);
clear ghano_ElNino ghanoDJF DJF_clim summClim gh250ln gh250 ghln* gh250ano gh250
%% Plot ERA ENSO COMPOSITE
figure(1)
subplot(2,2,1)

m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonm,latm,ghano_ENSO',[-200,-100:10:100,200], 'LineColor','none')
m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
title ('EL Nino Comp - ERA')

% fnm = ['ERA_ElNino_SH.png'];
% print ('-r300', fnm, '-dpng')

%% loading step average for models


stp_ens = zeros (mt, nt);
ghENSO_ens = zeros (mt, nt);
% idc = input ('Which model is it?? 1 for CCSM4 :: 2 for IPSL :: 3 for MIROC5 :: 4 for HadGEM2A :: 5 for GFDLCM3 :: 6 for MPIesmMR :: 7 for MRI-CGCM3 :: 8 for ACCESS1-3 :: 9 for NorESM1 ::');
k = 0;
for idc = 8:8
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

    % Loading step response to ENSO for the model
    stpsv = [modl, '_DJF_1980_2004.mat'];
    load (stpsv)
    lons = lon; lats = lat;
    [X,Y] = meshgrid(lats, lons);

    % Interpolate from model grid to target NorESM1 grid (taken as the target grid for interpolation)
    stpavg1 = interp2(X,Y,stpavg',Xq,Yq);
%     'size of a1-rtansposed ', size(stpavg1')

    stp_ens = stp_ens+stpavg1';
    
    % for individual model LRTM vs ERA correlation
    bstp = stpavg1; bstp = bstp(:, lt); 
    ['corr for ',modl, ' LRTM v ERA']
    nancorr2(bstp,aera)
    
    
%#####################################################################################################
%#####################################################################################################

    % MODEL COMPOSITE
    % loading sig_mat for the model to get ENSO composite for the model
    fl = ['sig_mat_DJF_1980_2004_',modl,'_amip.mat'];
    load (fl)

    [mc, nc, oc] = size(sig_mat);
    ghENSO_mod = zeros (nc, oc);

    for i = ensInd
        a = nanmean(sig_mat((i-1)*90+1:(i-1)*90+1+90, : ,:),1);
        % Interpolate from model grid to target NorESM1 grid (taken as the target grid for interpolation)
        
        ghENSO_mod = ghENSO_mod+squeeze(a);
    end
    
    ghENSO_mod1 = ghENSO_mod/length(ensInd);
    ghENSO_mod2 = squeeze(ghENSO_mod1);

    a1 = interp2(X,Y,ghENSO_mod2',Xq,Yq);
%     'size of a1-rtansposed ', size(a1')
    ghENSO_ens = ghENSO_ens + a1';
  
    
    % for individual model composite vs ERA correlation
    cmod = a1; cmod = cmod(:, lt); 
    ['corr for ',modl, ' ENSO comp v ERA']
    nancorr2(cmod,aera)
    
    k = k+1;
end

stp_ens = stp_ens/k;
ghENSO_ens = ghENSO_ens/k;

% ['shape of Model ensemble stpavg :: ']
% size(stp_ens)


% ['shape of Model ensemble stpavg :: ']
% size(ghENSO_ens)


%%
% Plot STP ENS
figure(1)
subplot(2,2,2)
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonm,latm,stp_ens,[-200,-100:10:100,200], 'LineColor','none')
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
m_contourf(lonm,latm,ghENSO_ens,[-200,-100:10:100,200], 'LineColor','none')
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
m_contourf(lonm,latm,ghENSO_ens,[-200,-100:10:100,200], 'LineColor','none')
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

lt = find(latm<=-20);
stp_ens = stp_ens(lt,:);
era = era (lt,:);
ghENSO_ens = ghENSO_ens(lt, :);

rspat1 = nancorr2(stp_ens,era)
rspat2 = nancorr2(ghENSO_ens, era)