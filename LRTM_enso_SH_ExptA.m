
%!/usr/bin/python

%{
% The code takes in:

Expt-A

1. forcing data (DJF + 40 lag days) : : projection of CMAP ENSO prec anom on
daily model prec anomaly

2. sigmal matrix (DJF all grid points over NH) : geopot. height anom @250
hPa from resp. models

Next, applies the LRT to produce Impulse and Step responses for different
lag days at each grid point over the NH

finds the significance by comparing against SRF generated by random AR1
instead of forcing prec ano

%}

clear;
close all


%  add m_map path
% addpath G:\MY-WORK-LIVE\RealProj_Jclim\codes-data_MATLAB\MATLAB_Aux_files\m_map\
% addpath G:\MY-WORK-LIVE\RealProj_Jclim\codes-data_MATLAB\sig_mat_files
addpath /home/pranab/Documents/MATLAB-AUX/m_map/

% addpath '~/Downloads/MATLAB_Aux_files/'
% addpath '~/Downloads/MATLAB_Aux_files/m_map/'
% % add data path
% addpath G:\RealProj_Jclim\codes-data\
% addpath G:\RealProj_Jclim\codes-data\frc_mat_files\
% addpath G:\RealProj_Jclim\codes-data\sig_mat_files\
% addpath G:\RealProj_Jclim\codes-data\MATLAB_Aux_files\
% addpath G:\RealProj_Jclim\codes-data\MATLAB_Aux_files\m_map\

klgnd=1;
lgt = 41;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ien = input ('Response to - (1) El Nino precipitation anom & (2) La Nina precipitation anom :: ');
if ien == 1
    trm = 'ElNino';
elseif ien == 2
    trm = 'LaNina';
end

idc = input ('Which model is it?? 1 for CCSM4 :: 2 for IPSL :: 3 for MIROC5 :: 4 for HadGEM2A :: 5 for GFDLCM3 :: 6 for MPIesmMR :: 7 for MRI-CGCM3 :: 8for ACCESS1-3 :: 9 for NorESM1 ::');
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

% Forcing data
fln = ['Proj_cmap',trm,'_frc_mat_lag40_DJF_1980_2004_',modl,'_amip.mat'];
load (fln)

% Signal data
fln = ['sig_mat_DJF_1980_2004_',modl,'_amip.mat'];
load (fln)

        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time lag or time lag Range::

itm1 = 31; itm2 = 41;

clim1 = -100; clim2 = 100;

% apply detrend again, summer trend could remain
frc_mat = detrend(frc_mat,1); 
[m,n] = size(frc_mat);
[ms,ns,os] = size(sig_mat);


%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply LRT
% forcing : independent variables (prec at diff time lags)
x1=frc_mat(:,n-lgt+1:n);
x2 = fliplr(x1); % flipping left to right so as to have max lag on the right-most column of forcing data (x here)
x=[ones(ms,1) x2];

stp = zeros(lgt, ns, os);
for iln = 1:os
    for ilt = 1:ns
        % signal : dependent variable
        y=sig_mat(:,ilt,iln);
    %     y=y';
            % apply detrend again, as python detrending does not work (Don't understand why!!) 
            % -> may be since we removed trend from the whole data, summer trend remains
        y = detrend(y);

            % whos x,y
        [bo, bint, r, rint, stats] = regress(y,x); % r --> residuals
%         whos r
        
%         b = bo(2:end); % first element of b is the intercept
%         b = bo(2:end)*std(x1(:,1)); % impulse response scaled by std of MC prec ano

        b = bo(2:end)*std(x1(:,1)); % impulse response scaled by a factor repreesentative of st dev in obs data

        stp(:,ilt,iln) = cumsum(b);

    end
end
% display 'St. dev. of trop prec. : ', std(x1(:,1))

% sct = squeeze(stp(itm,:,:)); % for specific time lag
stpavg = squeeze(nanmean(stp(itm1:itm2,:,:),1)); % given time lag range

% %%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
% 
% [mg, ng] = size(stpavg);
% 
% for i = 1:ng
%     stpavgA(:, i) = stpavg (:, i) - nanmean(stpavg(:, i));
% end
% 
% stpavg = stpavgA;


stpsv = ['stp_',modl, '_DJF_1980_2004_ExptA_',trm];


 %%%%%%%%%%% Remove the zonal mean to eliminate SAM %%%%%%%%%%%
gzm = stpavg;

[ma, na] = size(gzm);

for i = 1:na
    gzm(:, i) = gzm (:, i) - nanmean(gzm(:, i));
end

stpavg = gzm;


% save(stpsv, 'stpavg', 'lon', 'lat')
%% Calculate Geostrophic wind
% eta=stpavg';
% %
% % Constants
% g=9.8; % gravity constant
% R=6400000; % Earth radius
% omega=2*pi/(24*60*60); % Earth rotation angle velocity
%  
% % Set grid
% %x=linspace(129.5,130.5,50);
% %y=linspace(36,37,50);
% [x y]=meshgrid(lon,lat);
% x=x';y=y';
%  
% % Set Coriolis force coefficients
% f=2*omega*sind(y);
% 
% u=zeros(size(eta));
% v=zeros(size(eta));
% 
% % Calculate geostrophic current using numerical method
% for i=2:size(x,1)-1
%     for j=2:size(y,2)-1
%         dx=(x(i+1,j)-x(i-1,j))*(R*cosd(y(i,j))*pi/180);
%         dy=(y(i,j+1)-y(i,j-1))*(R*pi/180);
%          
%         v(i,j)=g/f(i,j)*(eta(i+1,j)-eta(i-1,j))/dx;
%         u(i,j)=-g/f(i,j)*(eta(i,j+1)-eta(i,j-1))/dy;
%         
%         if f(i,j) == 0
%             u(i,j) = NaN; v(i,j) = NaN;
%         end
%         
%         if abs(y(i,j)) <= 15
%             u(i,j) = NaN; v(i,j) = NaN;
%         end
%         
%     end
% end


%% PLOTS %%%%%%%%%%%%%%%%%%
figure

m_proj('stereographic','lat',-90,'long',0,'radius',91);%91);

m_contourf(lon,lat,stpavg,[-200,-100:10:100,200], 'LineColor','none')

m_grid('xtick',12,'tickdir','out','ytick',[0 30 60 90],'linest','--', 'fontsize',14,'color',[0.5 0.4 0.7])%,'color','b');
% m_coast('patch',[.7 .7 .7],'edgecolor','r');
% m_coast('color',[0.5 0.4 0.7]);
m_coast('color','r');


% caxis
% caxis([-120 120])
caxis([-100 100])
% caxis([-80 80])

% edited to remove white part of jet colorbar
% c = colormap(jet(20));
% c(9:12,:) = [];

% blue-white-red
% c = colormap(bluewhitered(16));
% c = colormap(jet(16));



% edited to remove the confusing/whitish part of rainbow clorbar
c = colormap(jet(28));
% if (idat == 3 && idc == 4)
%     c = colormap(bluewhitered(28));
% end
c(11:18,:) = [];

colormap(c)

% add colorbar
%cH = colorbar;
%set(cH,'FontSize',12, 'Fontweight','b');

% hold on
% m_quiver(x(1:2:end, 1:2:end),y(1:2:end, 1:2:end),u(1:2:end, 1:2:end),v(1:2:end, 1:2:end),0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title
% if idat == 1
%     title(['SRF TRMM ', yrst(1:4),'-',yrst(6:end), ' lags: ', num2str(itm1-1),'-', num2str(itm2-1), ' days @ Ph: ', num2str(iph)])
% elseif idat == 2
%     title(['SRF ', mdnm(1:end-18),' lags: ', num2str(itm1-1),'-', num2str(itm2-1), ' days @ Ph: ', num2str(iph)])
% end
title ([modl,' ',trm])

%fnm = ['SRF_AVG_30_40_',trm,'_SH_',yrst,'.png'];
fnm = ['SRF_AVG_30_40_',trm,'_SH_',modl,'_',trm,'.png'];
% print ('-r300', fnm, '-dpng')
