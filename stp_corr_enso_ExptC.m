%{

Details of Expt-C can be found here:

https://docs.google.com/document/d/1kU337S6c0JNcKNqAqmFwEGcnsaxZ41Rgk9MpLqnXuIs/edit

%}

clc; clear; close;
addpath /home/pranab/Documents/MATLAB-AUX/m_map/
addpath /home/pranab/Documents/MATLAB-AUX/ClimateDataToolbox/cdt/

%%%%% YEARS/ El Nino or La Nina??
yr = 1980:2005; rst = 'NCEP:1980-2005';
ien = input ('Response to - (1) El Nino precipitation anom & (2) La Nina precipitation anom :: ');
if ien == 1
    trm = 'ElNino';
    yr1 = [1982, 1987, 1991, 1997, 2002];
%     yr1 = [1986, 1987, 1991, 1997, 2002]; % leaving aside 1982 - the
%     anomalous ENSO year

%     yr1 = [2002, 2004, 2006, 2009, 2015]; % JClim paper ENSO
elseif ien == 2
    trm = 'LaNina';
    yr1 = [1988, 1995, 1998, 1999, 2000];
end

% calculate NCEP DJF ENSO composite
load NcepZG_sig_mat_DJF_1980_2005NH_NoFilt_Balaji



k = 1;
for i = 1:90:length(sig_mat)
    hgtDJF(k,:,:) = nanmean(sig_mat(i:i+89,:,:),1);
    k = k+1;
end

% extract mean anom for DJF during ENSO years
k2 = 1;
for j5 = 1:length(yr1)
    ensInd(k2) = find (yr == yr1(j5));
    k2 = k2+1;
end

yr(ensInd)

ghano_ENSO = squeeze(nanmean(hgtDJF(ensInd,:,:),1));
ghano_ENSO = ghano_ENSO';

% removing the lon=180 which appears twice as -180/+180 grid is converted to 0/360 grid
p180 = find(lon == 180);
ghano_ENSO(p180(1),:)=[]; lon(p180(1))=[]; clear p180

%%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
[mg, ng] = size(ghano_ENSO);
for i = 1:ng
    ghano_ENSO(:, i) = ghano_ENSO (:, i) - nanmean(ghano_ENSO(:, i));
end


% [X,Y] = meshgrid(lat, lon);
% % following interp step is important to avoid two 180 deg lons appearing in lon
% lonT = 0:2.5:360; latT = lat; [X360,Y360] = meshgrid(latT, lonT); % TARGET GRID for interpolation
% a1 = interp2(X,Y,ghano_ENSO,X360,Y360); % interpolation from X,Y to Xq,Yq
% clear ghano_ENSO; ghano_ENSO = a1; clear a1;




%%
%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(2,2,1)
[ln,lt] = meshgrid(lon,lat);
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(ln,lt,ghano_ENSO',[-200,-100:10:100,200], 'LineColor','none')
% m_contourf(lat,lon,ghano_ENSO,[-200,-100:10:100,200], 'LineColor','none')

m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
title ([trm,' Comp - ', rst])

    % TARGET GRID %
latR = lat; lonR = lon;
[Xq360,Yq360] = meshgrid(latR, lonR); % TARGET GRID for interpolation

% for correlation with stp ==> LRTM [0 --> +360 lon range]
ghano_ENSO360 = ghano_ENSO;
lt = find (latR<=-20);
aera360 = ghano_ENSO360(:, lt);

%% loading step average for models

stp_ens = zeros (mg, ng);
ghENSO_ens = zeros (mg, ng);
% idc = input ('Which model is it?? 1 for CCSM4 :: 2 for IPSL :: 3 for MIROC5 :: 4 for HadGEM2A :: 5 for GFDLCM3 :: 6 for MPIesmMR :: 7 for MRI-CGCM3 :: 8 for ACCESS1-3 :: 9 for NorESM1 ::');
k = 0;
rmod = []; rLRTM = [];
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading step response to ENSO for the model ///// has -180 --> +180 lon range
    stpsv = ['stp_',modl, '_DJF_1980_2004_ExptC_',trm]
    load (stpsv)
    p180 = find(lon == 180);
    stpavg(:, p180(1))=[]; lon(p180(1))=[]; clear p180
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % uncomment when using NON-balaji sig_mat (lon: -180E <--> 180E)
    % MAKING LON 0==>360 /// This step is needed because sig_mat of the models
    % are on 0-->360 RANGE
%     lon(lon<0) = 180 + lon(lon>0);
%     pln = find(lon>180); 
%     ln1 = lon(pln);
%     pln2 = find(lon<=180);
%     ln2 = lon(pln2);
%     lon = [ln2 ln1];
%     ghln = stpavg(:,pln); ghln2 = stpavg(:,pln2);
%     stp360 = cat(2, ghln2, ghln);
%     stp360 = stp360';

    stp360 = stpavg'; % uncomment when using balaji sig_mat (lon: 0-360E)
%     %%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
    [mg, ng] = size(stp360);
%     for i = 1:ng
%         stp360(:, i) = stp360 (:, i) - nanmean(stp360(:, i));
%     end

    stp_ens = stp_ens + stp360;
 
    %%%%%% CORR %%%%%%
    lonm = lon; latm = lat;
    lt = find (latm<=-20);
    bstp = stp360(:,lt);
    
    % %%%%%%%%%% AREA WEIGHTED correlation %%%%%%%%%%%%
    % % grid area
    [Xg,Yg] = meshgrid(latR(lt),lonR);
    A = cdtarea(Xg, Yg, 'km2');
    % normalization of area weights #1
    normA = A - min(A(:));
    normA = normA ./ max(normA(:));

    % checking correlation of individual models
%     ['corr for ',modl, ' LRTM v NCEP']
    
    r1 = nancorr2(bstp.*normA,aera360.*normA);
    rLRTM = [rLRTM; r1];
%#####################################################################################################
%#####################################################################################################

    % MODEL COMPOSITE ///// has 0 --> +360 lon range
    % loading sig_mat for the model to get ENSO composite for the model
    fl = ['sig_mat_DJF_1980_2004_',modl,'_amip.mat'];
    load (fl)
    [X,Y] = meshgrid(lat, lon); clear lon lat % model grid
  
    [mc, nc, oc] = size(sig_mat);
    ghENSO_mod = zeros (nc, oc);
    for i = ensInd
        a = nanmean(sig_mat((i-1)*90+1:(i-1)*90+1+90, : ,:),1);
        ghENSO_mod = ghENSO_mod+squeeze(a);
    end
    
    ghENSO_mod1 = ghENSO_mod/length(ensInd);
    ghENSO_mod2 = squeeze(ghENSO_mod1);

      %%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
    gzm = ghENSO_mod2';
    [ma, na] = size(gzm);
    for i = 1:na
        gzm(:, i) = gzm (:, i) - nanmean(gzm(:, i));
    end
    ghENSO_mod3 = gzm;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interpolate from model grid to target NorESM1 grid (taken as the target grid for interpolation)
    a1 = interp2(X,Y,ghENSO_mod3,Xq360,Yq360); % interpolation from X,Y to Xq,Yq
%     'size of a1-rtansposed ', size(a1')
    ghENSO_ens = ghENSO_ens + a1;
  
    
    % for individual model composite vs ERA correlation
    cmod = a1; cmod = cmod(:, lt); 
%     ['corr for ',modl, ' ENSO comp v NCEP']
    r2 = nancorr2(cmod.*normA,aera360.*normA);
    rmod = [rmod; r2];
    
    k = k+1;
end

stp_ens = stp_ens/k;
ghENSO_ens = ghENSO_ens/k;

% ['shape of Model ensemble stpavg :: ']
% size(stp_ens)
% 
% 
% ['shape of Model ensemble stpavg :: ']
% size(ghENSO_ens)

%% Spatial correlation between stpavg and ghano_ENSO(ERA)
lt = find(latR<=-20);
lat = latR(lt); lon = lonR;
lrtm_ens = stp_ens(:,lt);
modl_ens = ghENSO_ens(:, lt);


% %%%%%%%%%% AREA WEIGHTED correlation %%%%%%%%%%%%
% % grid area
[Xg,Yg] = meshgrid(lat, lon);
A = cdtarea(Xg, Yg, 'km2');

% normalization of area weights #1
normA = A - min(A(:));
normA = normA ./ max(normA(:));

% 
% % normalization of area weights #2
% NormRows = sqrt(sum(A.*A,2));
% normA = bsxfun(@rdivide,abs(A),NormRows);

lrtm_ens = lrtm_ens.*normA;
modl_ens = modl_ens.*normA;
aera360 = aera360.*normA;

whos lon lat lrtm* modl_ens aera360 normA

rspat1 = nancorr2(lrtm_ens,aera360);
rLRTM = [rLRTM; rspat1];
rspat2 = nancorr2(modl_ens, aera360);
rmod = [rmod; rspat2];
%%
% Plot STP ENS
figure(1)
subplot(2,2,2)
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
% m_contourf(Xq360,Yq360,stp_ens,[-200,-100:10:100,200], 'LineColor','none')
m_contourf(lonm,latm,stp_ens',[-200,-100:10:100,200], 'LineColor','none')

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
% m_contourf(Xq360,Yq360,ghENSO_ens,[-200,-100:10:100,200], 'LineColor','none')
m_contourf(lonm,latm,ghENSO_ens',[-200,-100:5:100,200], 'LineColor','none')

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
m_contourf(lonm,latm,ghENSO_ens',[-200,-100:10:100,200], 'LineColor','none')
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
title ('Expt - C')

fnm = ['ExptC_CORR_SH_',trm,'.png'];
% print ('-r300', fnm, '-dpng')

[rmod rLRTM]
fln = ['CORR_SH_ExptC_',trm,'.mat']
save(fln, 'rmod', 'rLRTM')
