%% ---------------------
%  hf4_sumfluxes
%  ---------------------
% Sum fluxes together
%% User Inputs

mnum   = [1:35,101:107] ;% Ensemble Members List
lags   = [1:3]          ;% Number of lags to include
monwin = 1              ;% Month Window
net    = 1              ;% Set to just output net heatflux. See note below-
% Note: Net outputs heat flux damping, test results, and correlation coefficients
% automatically, regardless of the toggle options below!

% Paths to data and output
datpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/03_hf3out/',num2str(monwin),'mon/'];
outpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/04_hf4out/',num2str(monwin),'mon/'];

% Enso Removal?
ensorem     = 1; % Set to 1 if ENSO was removed
savedamping = 0; % Set to 1 to save damping
savesigtest = 0; % Set to 1 to save significance testing
removepts   = 0; % 1 for remove sst fails, 2 for remove flx-sst crosscor fail

%% Script Start

% -----------------------------------------
% Preallocate (LON x LAT x ENS x MON x LAG)
% -----------------------------------------
pasize = [288,192,42,12,length(lags)];

if savedamping == 1
    dshflx = NaN(pasize);
    dlhflx = NaN(pasize);
    dfsns  = NaN(pasize);
    dflns  = NaN(pasize);
end
if savesigtest == 1
    tlhflx = NaN(pasize);
    tshflx = NaN(pasize);
    tfsns  = NaN(pasize);
    tflns  = NaN(pasize);
    tsst   = NaN(pasize);
end
if net == 1
    dnhflx = NaN(pasize);
    tnhflx = NaN(pasize);
    tsst   = NaN(pasize);
    rnhflx = NaN(pasize);
    rsst   = NaN(pasize);
end

% Loop by ensemble member to combine heat fluxes
for n = 1:length(mnum)
    ensnum = mnum(n);

    % Load and sum damping variables
    if savedamping == 1
        
        % Construct matfile name (ex: hfdamping_ENS001_ensorem1_monwin1.mat)
        matname = [datpath,'hfdamping_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
         
        % Load in variables
        load(matname)
        
        % Directly Read in damping variables
        dlhflx(:,:,n,:,:) = LHFLX_damping(:,:,:,1:3);
        dshflx(:,:,n,:,:) = SHFLX_damping(:,:,:,1:3);
        dfsns(:,:,n,:,:)  = FSNS_damping(:,:,:,1:3);
        dflns(:,:,n,:,:)  = FLNS_damping(:,:,:,1:3); 
    end
        
    %% Sum and save significance test results
    if savesigtest == 1
        
        % Construct matfile name (ex: corr_ENS001_ensorem1_monwin1.mat)
        matname = [datpath,'corr_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        load(matname)
        
        % Read in significance test results(accidental, extra dim at the
        % end of each test variable, corresponding to ensemble number
        tlhflx(:,:,n,:,:) = LHFLX_test(:,:,:,1:3,n);
        tshflx(:,:,n,:,:) = SHFLX_test(:,:,:,1:3,n);
        tfsns(:,:,n,:,:)  = FSNS_test(:,:,:,1:3,n);
        tflns(:,:,n,:,:)  = FLNS_test(:,:,:,1:3,n);
        tsst(:,:,n,:,:)   = SST_test(:,:,:,1:3);
    end
    
    % For net heat flux, automatically save the net heat fluxes,
    % significance testing results, and the correlation coefficients
    if net == 1
        
        % Load Net Heat Flux Damping (ex. nhfdamping_ENS001_ensorem1_monwin1.mat)
        matname =  [datpath,'nhfdamping_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        load(matname)
        
        % Load Test Results and CorrCoeffs (ex. ncorr_ENS001_ensorem1_monwin1.mat)
        matname = [datpath,'ncorr_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        load(matname)
        
        
        dnhflx(:,:,n,:,:) = NHFLX_damping(:,:,:,1:3); % Net Heat Flux Damping
        tnhflx(:,:,n,:,:) = NHFLX_test(:,:,:,1:3);    % Significance Test Results (Cross-Correlation)
        rnhflx(:,:,n,:,:) = NHFLX_rho(:,:,:,1:3);     % Cross-Correlation Coefficients
        tsst(:,:,n,:,:)   = SST_test(:,:,:,1:3);      % Significance Test Results (Autocorrelation)
        rsst(:,:,n,:,:)   = SST_rho(:,:,:,1:3);       % Autocorrelation Coefficients
    end
    
   fprintf('\nRead in data for ENS%s',num2str(ensnum,'%03d'))
end


    
% Remove points that failed sigtest if option is set
if removepts ~=0 && savedamping == 1
    matname = [datpath,'corr_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
    % temporarily force it to read in 1 month daya
%     datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/03_hf3out/1mon/';
%     matname = [datpath,'corr_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin3.mat'];
%     
    load(matname)
           
    % Remove points that failed SST autocorrelation
    dlhflx = dlhflx .* tsst;
    a = dlhflx .* tsst;
    dshflx = dshflx .* tsst;
    dfsns  = dfsns .* tsst;
    dflns  = dflns .* tsst;

    if removepts == 2
        dlhflx = dlhflx .* tlhflx;
        dshflx = dshflx .* tshflx;
        dfsns  = dfsns .* tfsns;
        dflns  = dflns .* tflns; 
    end                
end
%     if length(lags) > 1
%         % Load in lags 2 and 3
%         matname = [datpath,'/hfdamping_ENS',num2str(ensnum,'%03d')];
%         load(matname)
% 
%         dlhflx(:,:,n,:,2:3) = LHFLX_damping(:,:,:,2:3);
%         dshflx(:,:,n,:,2:3) = SHFLX_damping(:,:,:,2:3);
%         dfsns(:,:,n,:,2:3)  = FSNS_damping(:,:,:,2:3);
%         dflns(:,:,n,:,2:3)  = FLNS_damping(:,:,:,2:3);
%     end
    
if savedamping == 1
    % Add together heat fluxes for cumulative values
    drhflx = -1*dfsns + dflns;
    dthflx = dlhflx + dshflx;
    dnhflx = drhflx + dthflx;
end

% Save variables
if length(lags) == 1
    outpath = strcat(outpath,'lag1/');
end

if ensorem == 0
    out1 = [outpath,'damping_netfluxes_noensorem_rmpts',num2str(removepts),'.mat'];
    out2 = [outpath,'damping_hfluxes_noensorem_rmpts',num2str(removepts),'.mat'];
    out3 = [outpath,'sigtest_hfluxes_noensorem_rmpts',num2str(removepts),'.mat'];
else
    out1 = [outpath,'damping_netfluxes_rmpts',num2str(removepts),'.mat'];
    out2 = [outpath,'damping_hfluxes_rmpts',num2str(removepts),'.mat'];
    out3 = [outpath,'sigtest_hfluxes_rmpts',num2str(removepts),'.mat'];
end

if savedamping == 1
    save(out1,'dnhflx','drhflx','dthflx')
    save(out2,'dlhflx','dshflx','dfsns','dflns')
end

if savesigtest == 1
    %save(out3,'tsst')
    save(out3,'tlhflx','tshflx','tfsns','tflns','tsst')   
end

% Save net damping file to "outpath/alldamping_nhflx_monwin1.mat"
if net == 1
    save([outpath,'alldamping_nhflx_monwin',num2str(monwin),'.mat'],'dnhflx','tnhflx','rnhflx','tsst','rsst')
end


