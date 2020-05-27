% Land Ice Mask Maker (CESM1)
% Make land/ice masks for CESM1 Data. This relies 
% upon the build_e11name.m function and follows the naming 
% convention of the CESM1.1 data on stormtrack
%
% Inputs
% 1) inarray    - input array, in lon x lat x time
% 2) rtype      - "PIC", "HTR" or "RCP85"
% 3) ensnum     - Number of the ensemble member
%
% Outputs
% 1) mat file with variable "mask_landice" in outloc for each ensemble member.
%% User Edits
rtype   = 'HTR'                         ;% Experiment (HTR, CTL, PIC)
mnum    = [1:35,101:107]                ;% Ensemble member numbers
outloc  = '/home/glliu/01_Data/masks/'  ;% Output Location
deg5    = 1                             ;% Set to 1 for 5 deg option

%% Script Start
% Timekeeping
allstart = datetime('now');
fprintf('\nStarting landicemaskmaker (%s)\n',datetime('now'))
fprintf('\tRtype: %s | deg5: %i \n',rtype,deg5)
fprintf('\tOutput Path: %s \n',outloc)

% Add path
addpath('/home/glliu/')
startup

if deg5 == 1
    
    % Load matfiles with smoothed data
    datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/5deg/';
    lndmat = 'HTR_5deg_LANDFRAC.mat';
    icemat = 'HTR_5deg_ICEFRAC.mat';
    load([datpath,lndmat])
    lndall = varnew;
    load([datpath,icemat])
    iceall  = varnew;
    
else
    
    % Get names
    lndNames = build_e11name(rtype,'LANDFRAC');
    iceNames = build_e11name(rtype,'ICEFRAC') ;
end

% Loop for each ensemble member
for idx = 1:length(mnum)
    %% Setup
    lstart = datetime('now');
    enum = mnum(idx);
    
    if deg5 == 1
        
        % Index for ensemble member
        lnd = lndall(:,:,:,idx);
        ice = iceall(:,:,:,idx);
        
        % Permute to [Lon Lat Time]
        lnd = permute(lnd,[2,3,1]);
        ice = permute(ice,[2,3,1]);
        
    else
        % Get appropriate netCDF names
        icename = iceNames(idx);
        lndname = lndNames(idx);

        % Add paths (ICEFRAC and LANDFRAC data on stormtrack)
        addpath(['/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/LANDFRAC/']);
        addpath(['/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/ICEFRAC/']);

        % Read in ice and land data
        lnd = ncread(lndname,'LANDFRAC');
        ice = ncread(icename,'ICEFRAC');
    end


    %% Make land mask -----
    % Get average landfracs
    meanLF = nanmean(lnd,3);
    maskland = meanLF;

    % Set all land portions (avg 30%+ landfrac to NaN)
    maskland(maskland > 0.3) = NaN;

    % Set all ocean (non-land) portions to 1.
    maskland(~isnan(maskland)) = 1;

    %% Make ice mask -----

    % Set all points, set to NaN if ice frac ever exceeds 5%
    maskice = ice;
    maskice(maskice > 0.05) = NaN;
    maskice = mean(maskice,3);
    maskice(~isnan(maskice)) = 1;

    %% Combine land_ice mask and save
    mask_landice = maskland.*maskice;
    if deg5 == 1
        % Add "_5deg"
        outfile = [outloc,'landicemask_ensnum',num2str(enum,'%03d'),'_5deg.mat'];
    else
        outfile = [outloc,'landicemask_ensnum',num2str(enum,'%03d'),'.mat'];
    end
    save(outfile,'mask_landice');
    
    % Timekeeping
    lend = datetime('now');
    lelapsed = lend - lstart;
    fprintf('Mask created for ensnum %s! (%s seconds)\n',num2str(enum,'%03d'),lelapsed)
end

 % Timekeeping
allend = datetime('now');
elapsed = allend - allstart;
fprintf('All masks created! (%s seconds)',elapsed)
fprintf('\nFinished running landicemaskmaker (%s)\n',datetime('now'))
