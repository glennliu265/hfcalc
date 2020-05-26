%% User Input

% Name of variable
varname = 'HMXL';
outpath = '/home/glliu/01_Data/Project12860/';

%% Startup
allstart = datetime('now');
fprintf(['Now Coarsening: ',varname])


if strcmp(varname,'SST') == 1
    ncpath = '/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/ocn/proc/tseries/monthly/SST/';
elseif strcmp(varname,'SSS') == 1
    ncpath = '/stormtrack/data0/yokwon/CESM1LE/SSS/';
elseif strcmp(varname,'HMXL') == 1
    ncpath = '/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/ocn/proc/tseries/monthly/HMXL/';
end

addpath(ncpath)
ncnames = dir([ncpath,'b.e11.B20TRC5CNBDRD.f09_g16.*',varname,'*.nc']);

% Remove OIC
%noOIC = structfun(@(x) ~any(contains(x,'OIC')),ncnames)
ncnames(43:end) = [];

% Create 5x5 Lat/Lon Grid
lat5 = [-90:5:90];
lon5 = [0:5:360];
[gx,gy] = meshgrid(lon5,lat5);

%% Loop

% Loop for each ncfile
for n = 1:length(ncnames)
    
    hfst = datetime('now');
    % Get ncfile name
    nc = [ncpath,'/',ncnames(n).name];
    
    % Read in Data
    var = ncread(nc,varname);
    tlat = ncread(nc,'TLAT');       % Array of t-grid latitudes
    tlon = ncread(nc,'TLONG');

    % Move time dimension to front to linearly index position
    rvar = permute(squeeze(var),[3,1,2]);
    
    % Create array to store ensemble average
    if n == 1
        varsum = zeros(1032,length(lon5),length(lat5));
    end

    % Remap data to 5x5 Grid (Average everything within)...
    varnew = NaN(size(rvar,1),length(lon5),length(lat5));
    icnt = 0;   
    % Loop for longitude
    for o = 1:length(lon5)-1
        % Bounding Lons
        x0 = lon5(o);
        x1 = lon5(o+1);
        
        % Loop for latitude
        for a = 1:length(lat5)-1
            icnt = icnt + 1;
%             if mod(icnt,1000) == 0 || icnt == 2590
%                 fprintf(['Processed ',num2str(icnt),...
%                     ' of ',num2str((length(lon5)-1)*(length(lat5)-1)-2),' points.\n'])
%             end

            % Bounding Lats
            y0 = lat5(a);
            y1 = lat5(a+1);

            % Find coordinates that satisfy the value
            kxy = find((tlon > x0) & (tlon < x1) & (tlat > y0) & (tlat < y1));
            tval = mean(rvar(:,kxy),2);

            % Assign NaN to empty regions
            try
                varnew(:,o,a) = tval;
            catch
                varnew(:,o,a) = NaN;
            end     
        end
    end
    
    % Add to variable
    if n==1
        varsum = varsum + varnew(841:end,:,:);
    else
        varsum = varsum + varnew;
    end
    
    % Save variable to matfile
    matname = [outpath,'HTR_5deg_',varname,'_ensnum',num2str(n,'%03d'),'.mat'];
    save(matname,'varnew','lat5','lon5');
    
     % Timekeeping
     hfend = datetime('now');
     elapsed = char(hfend-hfst);
     fprintf('\n Output file for ensnum %s in %s seconds\n',num2str(n,'%03d'),elapsed)
end

% Take Ensemble Mean (Post Calculation) 
var_ensmean = varsum./length(ncnames);

% Output Ensmean and variables
matname = [outpath,'HTR_5deg_',varname,'_ensnumavg.mat'];
save(matname,'var_ensmean','lat5','lon5');

% Timekeeping
allend = datetime('now');
elapsed = char(allend-allstart);
fprintf('\n Completed script in %s seconds', elapsed)
