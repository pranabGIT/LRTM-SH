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

% ncfile = '/media/pranab/Backup Plus/Backup_12_05_2021/hgt_NCEP_NCAR_1980JAN_2020DEC_monthly.nc'; rst = 'NCEP:1980-2019'; yr = 1980:2019; % Dec - centred == 1980 means DJF/D1980/J81/F81
ncfile = '/media/pranab/Backup Plus/Backup_12_05_2021/hgt_NCEP_NCAR_1980JAN_2005DEC_monthly.nc'; rst = 'NCEP:1980-2005'; yr = 1980:2005;

hgt = ncread(ncfile,'HGT'); hgt = squeeze(hgt); lon = ncread(ncfile,'LONN71_73'); lat = ncread(ncfile,'LAT');
lt = find(lat<=0); hgt = hgt(:,lt,:); lat = lat(lt); [m,n,o] = size(hgt);
%%
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
%%
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

