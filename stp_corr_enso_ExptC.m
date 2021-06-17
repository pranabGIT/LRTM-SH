%{
Details of Expt-A can be found here:

https://docs.google.com/document/d/1kU337S6c0JNcKNqAqmFwEGcnsaxZ41Rgk9MpLqnXuIs/edit

%}
clc; clear
addpath /home/pranab/Documents/MATLAB-AUX/m_map/

%%%%% YEARS/ El Nino or La Nina??
ien = input ('Response to - (1) El Nino precipitation anom & (2) La Nina precipitation anom :: ');
if ien == 1
    trm = 'ElNino';
    yr1 = [1982, 1987, 1991, 1997, 2002];
%     yr1 = [2002, 2004, 2006, 2009, 2015]; % JClim paper ENSO
elseif ien == 2
    trm = 'LaNina';
    yr1 = [1988, 1995, 1998, 1999, 2000]; 
end

% calculate NCEP DJF ENSO composite

% ncfile = '/media/pranab/Backup Plus/Backup_12_05_2021/NCEP_gh/hgt_NCEP_NCAR_1980JAN_2020DEC_monthly.nc'; rst = 'NCEP:1980-2019'; yr = 1980:2019; % Dec - centred == 1980 means DJF/D1980/J81/F81
ncfile = '/media/pranab/Backup Plus/Backup_12_05_2021/NCEP_gh/hgt250_NCEP_NCAR_1980JAN_2005DEC_monthly.nc'; rst = 'NCEP:1980-2005'; yr = 1980:2005;

hgt = ncread(ncfile,'HGT'); hgt = squeeze(hgt); lon = ncread(ncfile,'LONN71_73'); lat = ncread(ncfile,'LAT');
lt = find(lat<=0); hgt = hgt(:,lt,:); lat = lat(lt); [m,n,o] = size(hgt);

 % Calculate climatology
b = zeros (m, n, 12);
for i = 1:12
    b(:,:,i) = nanmean(hgt(:,:,i:i+12:end), 3);
end
hgtclim = repmat(b,[1,1,o/12]);
% Calculate anomaly
hgtano = hgt-hgtclim; clear hgtclim hgt

% detrend anomaly
hgtd = [];
for i = 1:m
    for j = 1:n
        hgtd(i,j,:) = detrend(squeeze(hgtano(i,j,:)),1);
    end
end
hgtano = hgtd;

% mean anomaly during DJF of each year
hgtDJF = zeros(o/12-1, m, n);
k=1;
for i = 12:12:o-1
    hgtDJF(k,:,:) = nanmean(hgtano(:,:,i:i+2),3);
    k=k+1;
end
% extract mean anom for DJF during ENSO years
k2 = 1;
for j5 = 1:length(yr1)
    ensInd(k2) = find (yr == yr1(j5));
    k2 = k2+1;
end

ghano_ENSO = squeeze(nanmean(hgtDJF(ensInd,:,:),1));
% ghano_ENSO = ghano_ENSO';

%%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
[mg, ng] = size(ghano_ENSO);
for i = 1:ng
    ghano_ENSO(:, i) = ghano_ENSO (:, i) - nanmean(ghano_ENSO(:, i));
end

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

% fnm = ['ERA_ElNino_SH.png'];
% print ('-r300', fnm, '-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKING LON 0==>360 /// This step is needed because sig_mat of the models
% are on 0-->360 RANGE
lon(lon<0) = 180 + lon(lon>0);
pln = find(lon>180); 
ln1 = lon(pln);
pln2 = find(lon<=180);
ln2 = lon(pln2);
lon = [ln2; ln1];
ghln = ghano_ENSO(pln,:); ghln2 = ghano_ENSO(pln2,:);
ghano_ENSO360 = cat(1, ghln2, ghln);

    % TARGET GRID %
latR = lat; lonR = lon;
[Xq360,Yq360] = meshgrid(latR, lonR); % TARGET GRID for interpolation

% for correlation with stp ==> LRTM [0 --> +360 lon range]
lt = find (latR<=-20);
aera360 = ghano_ENSO360(:, lt);
for i=1:1 % FOR LOOP JUST FOR MINIMISING THIS PART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % clc; clear
% % % % % 
% % % % % addpath /home/pranab/Documents/MATLAB-AUX/m_map/
% % % % % 
% % % % % %%%%% YEARS/ El Nino or La Nina??
% % % % % yr = 1980:2003; % Dec - centred == 1980 means DJF/D1980/J81/F81
% % % % % 
% % % % % ien = input ('Response to - (1) El Nino precipitation anom & (2) La Nina precipitation anom :: ');
% % % % % if ien == 1
% % % % %     trm = 'ElNino';
% % % % %     % ElNino
% % % % %     yr1 = [1982, 1987, 1991, 1997, 2002];
% % % % % 
% % % % % elseif ien == 2
% % % % %     trm = 'LaNina';
% % % % %     % list for La Nina
% % % % %     yr1 = [1988, 1995, 1998, 1999, 2000]; 
% % % % % 
% % % % % end
% % % % % 
% % % % % %%
% % % % % % calculate NCEP DJF ENSO composite
% % % % % 
% % % % % load NcepZG_sig_mat_DJF_1980_2004SH_NoFilt.mat
% % % % % rst = 'NCEP';
% % % % % 
% % % % % za = sig_mat; [mt, nt, ot] = size(za);
% % % % % % clear lon lat; % to be used to interpolate stpavg data to this grid
% % % % % 
% % % % % zaDJF = zeros(mt/90, nt, ot);
% % % % % 
% % % % % k = 1;
% % % % % for i = 1:90:length(za)
% % % % %     zaDJF(k,:,:) = nanmean(za(i:i+89,:,:),1);
% % % % %     i, i+89
% % % % %     k = k + 1;
% % % % % end
% % % % % 
% % % % % k2 = 1;
% % % % % for j5 = 1:length(yr1)
% % % % %     ensInd(k2) = find (yr == yr1(j5));
% % % % %     k2 = k2+1;
% % % % % end
% % % % % 
% % % % % ghano_ENSO = squeeze(nanmean(zaDJF(ensInd,:,:),1));
% % % % % ghano_ENSO = ghano_ENSO';
% % % % % 
% % % % % %%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
% % % % % [mg, ng] = size(ghano_ENSO);
% % % % % for i = 1:ng
% % % % %     ghano_ENSO(:, i) = ghano_ENSO (:, i) - nanmean(ghano_ENSO(:, i));
% % % % % %     ghano_ENSO(:, i) = ghano_ENSO (:, i) - nanmean(nanmean(ghano_ENSO,1),2);
% % % % % end
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % MAKING LON 0==>360 /// This step is needed because sig_mat of the models
% % % % % % are on 0-->360 RANGE
% % % % % lon(lon<0) = 180 + lon(lon>0);
% % % % % pln = find(lon>180); 
% % % % % ln1 = lon(pln);
% % % % % pln2 = find(lon<=180);
% % % % % ln2 = lon(pln2);
% % % % % lon = [ln2 ln1];
% % % % % ghln = ghano_ENSO(pln,:); ghln2 = ghano_ENSO(pln2,:);
% % % % % ghano_ENSO360 = cat(1, ghln2, ghln);
% % % % % 
% % % % %     % TARGET GRID %
% % % % % latR = lat; lonR = lon;
% % % % % [Xq360,Yq360] = meshgrid(latR, lonR); % TARGET GRID for interpolation
% % % % % 
% % % % % % for correlation with stp ==> LRTM [0 --> +360 lon range]
% % % % % lt = find (latR<=-20);
% % % % % aera360 = ghano_ENSO360(:, lt);
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %%%%%%%%%% Plot Reanalysis ENSO COMPOSITE %%%%%%%%%%
% % % % % figure(1)
% % % % % subplot(2,2,1)
% % % % % 
% % % % % m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
% % % % % m_contourf(lonR,latR,ghano_ENSO360',[-200,-100:10:100,200], 'LineColor','none')
% % % % % m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
% % % % % m_coast('color','r');
% % % % % % caxis([-200 200])
% % % % % caxis([-100 100])
% % % % % 
% % % % % % edited to remove the confusing/whitish part of rainbow clorbar
% % % % % c = colormap(jet(28));
% % % % % c(11:18,:) = [];
% % % % % colormap(c)
% % % % % title ([trm,' Comp - ', rst])
% % % % % 
% % % % % % fnm = ['ERA_ElNino_SH.png'];
% % % % % % print ('-r300', fnm, '-dpng')
end
%% loading step average for models

stp_ens = zeros (m, n);
ghENSO_ens = zeros (m, n);
% idc = input ('Which model is it?? 1 for CCSM4 :: 2 for IPSL :: 3 for MIROC5 :: 4 for HadGEM2A :: 5 for GFDLCM3 :: 6 for MPIesmMR :: 7 for MRI-CGCM3 :: 8 for ACCESS1-3 :: 9 for NorESM1 ::');
k = 0;
rLRTM = []; rmod = [];
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
    stpsv = ['stp_',modl, '_DJF_1980_2004_ExptC_',trm];
    load (stpsv)
    [X,Y] = meshgrid(lat, lon);
 
    a1 = interp2(X,Y,stpavg',Xq360,Yq360); % interpolation from X,Y to Xq,Yq

    %%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
    [mg, ng] = size(a1);
    for i = 1:ng
        a1(:, i) = a1 (:, i) - nanmean(a1(:, i));
    end

    stp_ens = stp_ens + a1;
   
    lt = find (latR<=-20);
    bstp = a1(:,lt);
    
    % %%%%%%%%%% AREA WEIGHTED correlation %%%%%%%%%%%%
    % % grid area
    [Xg,Yg] = meshgrid(latR(lt),lonR);
    A = cdtarea(Xg, Yg, 'km2');
    % normalization of area weights #1
    normA = A - min(A(:));
    normA = normA ./ max(normA(:));
    
    % checking correlation of individual models
%     ['corr for ',modl, ' LRTM v NCEP'];
    
   r1 =  nancorr2(bstp.*normA,aera360.*normA);
   rLRTM = [rLRTM; r1];    
    
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

      %%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
    gzm = ghENSO_mod2';
    
    [ma, na] = size(gzm);

    for i = 1:na
        gzm(:, i) = gzm (:, i) - nanmean(gzm(:, i));
    end

    ghENSO_mod3 = gzm;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    a1 = interp2(X,Y,ghENSO_mod3,Xq360,Yq360); % interpolation from X,Y to Xq,Yq
    
    
    ghENSO_ens = ghENSO_ens + a1;
  
    % Correlation for individual model composite vs ERA correlation
    cmod = a1; cmod = cmod(:, lt); 
%     ['corr for ',modl, ' ENSO comp v NCEP']
    r2 = nancorr2(cmod.*normA,aera360.*normA);
    rmod = [rmod; r2];
    
    k = k+1;
end

stp_ens = stp_ens/k;
ghENSO_ens = ghENSO_ens/k;

% ['shape of LRTM ensemble stpavg :: ']
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
% normB = bsxfun(@rdivide,abs(A),NormRows);

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
m_contourf(lonR,latR,stp_ens',[-200,-100:10:100,200], 'LineColor','none')
m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
title (['step avg. ensemble - ', trm])

%% Plot MODEL GH ENSO COMPOSITE
figure(1)
subplot(2,2,3)
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonR,latR,ghENSO_ens',[-200,-100:10:100,200], 'LineColor','none')
% m_contourf(lonm,latm,ghENSO_ens,[-200,-100:5:100,200], 'LineColor','none')

m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',8,'color',[0.5 0.4 0.7])%,'color','b');
m_coast('color','r');
% caxis([-200 200])
caxis([-100 100])

% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
c(11:18,:) = [];
colormap(c)
title (['model ghano ensemble - ', trm])

% 
% This following one is to produce the colormap and delete the spatial plt
% manually
subplot(2,2,4)
m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);
m_contourf(lonR,latR,ghENSO_ens',[-200,-100:10:100,200], 'LineColor','none')
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
print ('-r300', fnm, '-dpng')

[rmod rLRTM]

fln = ['CORR_SH_ExptC_',trm,'.mat']
save(fln, 'rmod', 'rLRTM')
