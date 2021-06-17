nccreate('test_files_1.nc','lat','Dimensions',{'lat' length(lat)});
ncwriteatt('test_files_1.nc', 'lat', 'bounds', 'lat_bnds');
ncwriteatt('test_files_1.nc', 'lat', 'standard_name', 'latitude');
ncwriteatt('test_files_1.nc', 'lat', 'long_name', 'latitude');
ncwriteatt('test_files_1.nc', 'lat', 'units', 'degrees north');
ncwriteatt('test_files_1.nc', 'lat', '_CoordinateAxisType', 'Lat');

nccreate('test_files_1.nc','lon','Dimensions',{'lon' length(lon)});
ncwriteatt('test_files_1.nc', 'lon', 'bounds', 'lon_bnds');
ncwriteatt('test_files_1.nc', 'lon', 'standard_name', 'longitude');
ncwriteatt('test_files_1.nc', 'lon', 'long_name', 'longitude');
ncwriteatt('test_files_1.nc', 'lon', 'units', 'degrees north');
ncwriteatt('test_files_1.nc', 'lon', '_CoordinateAxisType', 'Lon');

% nccreate('test_files.nc','time','Dimensions',{'time' 9862});
% ncwriteatt('test_files.nc', 'time', 'long_name', 'Time variable');
% ncwriteatt('test_files.nc', 'time', 'units', 'days since 1979-01-01 00:00:00');
% ncwriteatt('test_files.nc', 'time', '_CoordinateAxisType', 'Time');
% nccreate('test_files.nc','lev','Dimensions',{'lev' 4});
% ncwriteatt('test_files.nc', 'lev', 'standard_name', 'pressure');
% ncwriteatt('test_files.nc', 'lev', 'long_name', 'pressure');
% ncwriteatt('test_files.nc', 'lev', 'units', 'Pa');
% ncwriteatt('test_files.nc', 'lev', '_CoordinateAxisType', 'Z');

nccreate('test_files_1.nc','zg','datatype','single','Dimensions',{'lon' length(lon) 'lat' length(lat)});
ncwriteatt('test_files_1.nc', 'zg', 'standard_name', 'geopotential_height');
ncwriteatt('test_files_1.nc', 'zg', 'long_name', 'Geopotential height');
ncwriteatt('test_files_1.nc', 'zg', 'units', 'm');
% ncwriteatt('test_files_1.nc', 'ncepenso', 'missing_value', '-1e+04');
% ncdisp('test_files.nc');

% nccreate('test_files_1.nc','rspat1','datatype','single','Dimensions',{'lon' 145 'lat' 37});
% ncwriteatt('test_files_1.nc', 'rspat1', 'standard_name', 'lrtm enso');
% ncwriteatt('test_files_1.nc', 'rspat1', 'long_name', 'ncep enso');
% ncwriteatt('test_files_1.nc', 'rspat1', 'units', 'm');
% ncwriteatt('test_files_1.nc', 'rspat1', 'missing_value', '-1e+04');
% 
% nccreate('test_files_1.nc','rspat2','datatype','single','Dimensions',{'lon' 145 'lat' 37});
% ncwriteatt('test_files_1.nc', 'rspat2', 'standard_name', 'model enso');
% ncwriteatt('test_files_1.nc', 'rspat2', 'long_name', 'ncep enso');
% ncwriteatt('test_files_1.nc', 'rspat2', 'units', 'm');
% ncwriteatt('test_files_1.nc', 'rspat2', 'missing_value', '-1e+04');





%%%%%%%%%%%%%%%% WRITING %%%%%%%%%%%%%%%%%%%%
ncwrite('test_files_1.nc','lat',lat);
ncwrite('test_files_1.nc','lon',lon);

ncwrite('test_files_1.nc','zg',aera360);

% ncwrite('test_files_1.nc','lrtm_ens',lrtm_ens);
% ncwrite('test_files_1.nc','lrtm_CCSM4',lrtm_CCSM4);
% ncwrite('test_files_1.nc','lrtm_ACCESS1_3',lrtm_ACCESS1_3);
% ncwrite('test_files_1.nc','lrtm_GFDLCM3',lrtm_GFDLCM3);
% ncwrite('test_files_1.nc','lrtm_HadGEM2A',lrtm_HadGEM2A);
% ncwrite('test_files_1.nc','lrtm_IPSLcm5aMR',lrtm_IPSLcm5aMR);
% ncwrite('test_files_1.nc','lrtm_MIROC5',lrtm_MIROC5);
% ncwrite('test_files_1.nc','lrtm_MPIesmMR',lrtm_MPIesmMR);
% ncwrite('test_files_1.nc','lrtm_MRICGCM3',lrtm_MRICGCM3);
% ncwrite('test_files_1.nc','lrtm_NorESM1',lrtm_NorESM1);
% 
% ncwrite('test_files_1.nc','modl_ens',modl_ens);


ncdisp('test_files_1.nc');
