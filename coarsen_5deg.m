%% User Input

% Name of variable (choices are below)
% ocn [HMXL, SST, SSS, IFRAC]
% atm [TS,ICEFRAC,LANDFRAC,SHFLX,LHFLX,FSNS,FLNS]
vars = {'SHFLX','LHFLX','FSNS','FLNS'};
outpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/5deg/';
domain = 'atm'                  ; %['ocn' or 'atm']

% Set path to data, based on domain
% dat_atm = '/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/';
% dat_ocn = '/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/ocn/proc/tseries/monthly/';
datpath = strcat('/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/',domain,'/proc/tseries/monthly/');

% Toggles
areaweight = 1; % Option to weight each cell by its area when averaging (ocn variables)

%% Startup
allstart = datetime('now');

for v = 1:length(vars)
    varname = char(vars(v));
    
    fprintf('Now Coarsening %s | Domain: %s | Area-Weight: %s \n',varname,domain,num2str(areaweight))
    fprintf('\t Output Location: %s',outpath)
    addpath('/home/glliu/')
    startup

    % --------------------------
    % Construct NC file location
    % --------------------------
    % Salinity is stored at a different location
    if strcmp(varname,'SSS') == 1
        ncpath = '/stormtrack/data0/yokwon/CESM1LE/SSS/';
    else
        ncpath = strcat(datpath,varname,'/');
    end

    % Add ncpath and pull names from directory using glob (Historical Exp)
    addpath(ncpath)
    ncnames = dir([ncpath,'b.e11.B20TRC5CNBDRD.f09_g16.*',varname,'*.nc']);

    % Remove OIC
    ncnames(43:end) = [];

    % ------------------------
    % Create 5x5 Lat/Lon Grid
    % ------------------------
    lat5 = [-90:5:90];
    lon5 = [-0:5:360];

    % Make new lat/lon indices (at midpoints)
    LAT = [-87.5:5:87.5];
    LON = [2.5:5:357.5];

    % ------------------------
    % Setup years
    % ------------------------
    yrs   = 1850:2005;

    % Get Index of start year for first ensemble member
    tstart = (find(yrs==1920)-1)*12+1;

    % Get total size of time dimension
    tsize = (yrs(end) - 1920 + 1) * 12;



    %% Loop
    % Preallocate [time x lon x lat x ens], values represent midpoint numbers
    varnew = NaN(tsize,length(lon5)-1,length(lat5)-1,length(ncnames));
    ANEW = NaN(length(lon5)-1,length(lat5)-1);

    % Loop for each ncfile
    for n = 1:length(ncnames)

        hfst = datetime('now');

        % Get ncfile name
        nc = [ncpath,ncnames(n).name];

        % Read in Data
        var = ncread(nc,varname);

        if strcmp(domain,'ocn') == 1
            tlat = ncread(nc,'TLAT');       % Array of t-grid latitudes
            tlon = ncread(nc,'TLONG');
        else
            tlat = ncread(nc,'lat'); 
            tlon =  ncread(nc,'lon'); 
        end

        % Multiply variables by the area of the cell
        if areaweight == 1

            % Ocean Grid, use area variable
            if strcmp(domain,'ocn') == 1

                tarea = ncread(nc,'TAREA');   % Area of T cells (cm2)
                var   = squeeze(var) .* tarea(:,:,ones(1,size(var,3)));         % Variable times Area of cell

            % Atm Grid, cosine of latitude
            elseif strcmp(domain,'atm') == 1

                % Create meshgrid and scale by cosine of latitude
                [X,Y] = meshgrid(tlon,tlat);
                %AREA  = sqrt(cos(Y/180*pi)); % Old method from Young-Oh
                tarea = cdtarea(Y,X);


                % Apply to variable
                var   = var .* permute(tarea(:,:,ones(1,size(var,3))),[2,1,3]); 
            end

        end

        % Reduce dimensions if ensemble member 1
        if n == 1
            var = var(:,:,tstart:end);
        end

        % Move time dimension to front to linearly index position
        rvar = permute(squeeze(var),[3,1,2]);

        % Remap data to 5x5 Grid (Average everything within)...  
        % Loop for longitude
        for o = 1:length(lon5)-1
            % Bounding Lons
            x0 = lon5(o);
            x1 = lon5(o+1);

            % Loop for latitude
            for a = 1:length(lat5)-1

                % Bounding Lats
                y0 = lat5(a);
                y1 = lat5(a+1);

                if strcmp(domain,'ocn') == 1
                    % Find coordinates that satisfy the value
                    kxy = find((tlon >= x0) & (tlon <= x1) & (tlat >= y0) & (tlat <= y1));

                    % Take mean of indexed coordinates
                    tval = rvar(:,kxy);
                    tmean = mean(tval,2);

                    %tval = mean(rvar(:,kxy),2);
                elseif strcmp(domain,'atm') == 1
                    kx = find((tlon >= x0) & (tlon <= x1));
                    ky = find((tlat >= y0) & (tlat <= y1));

                    % Get values, combine spatial index and take mean
                    tval = rvar(:,kx,ky);
                    tmean = mean(tval(:,:),2);

                end




                % Divide by total area for weighting (ocean grid)
                if areaweight == 1
                    if strcmp(domain,'ocn') == 1

                        % Sum area of all points
                        totalarea = sum(tarea(kxy));

                        % Divide values to get weighted values
                        tval = tval ./ totalarea;

                        % Sum to get mean
                        tmean = sum(tval,2)

                    end

                    if strcmp(domain,'atm') == 1

                        % Sum area of all points
                        totalarea = sum(sum(tarea(ky,kx),2),1);

                        % Divide values to get weighted values
                        tval = tval ./ totalarea;  

                        % Sum for weighted area average
                        tmean = sum(tval(:,:),2);
                        ANEW(o,a) = totalarea;

                    end

                else


                end

                % Assign NaN to empty regions
                try
                    varnew(:,o,a,n) = tmean;
                catch
                    varnew(:,o,a,n) = NaN;
                end     
            end
        end

         % Timekeeping
         hfend = datetime('now');
         elapsed = char(hfend-hfst);
         totelapsed = char(hfend-allstart);
         fprintf('\n Processed for ensnum %s (%s | %s)',num2str(n,'%03d'),elapsed,totelapsed)
    end



    % Save variable to matfile
    matname = [outpath,'HTR_5deg_',varname,'_ensnum.mat'];
    save(matname,'varnew','LAT','LON');

    % Timekeeping
    allend = datetime('now');
    elapsed = char(allend-allstart);
    fprintf('\n Completed %s in %s seconds',varname,elapsed)
end
% Timekeeping
allend = datetime('now');
elapsed = char(allend-allstart);
fprintf('\n Completed script in %s seconds', elapsed)