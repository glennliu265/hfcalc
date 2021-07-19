
%HF Calc Ens Average

% Calculate Ensemble Average 

%% ------------------------------
%  User Input
%  ------------------------------
%fluxes = {'RHFLX','THFLX','NHFLX','SHFLX','LHFLX','FSNS','FLNS'};
fluxes   = {'nhflx','thflx','rhflx'}  ; %Flux to calculate ensemble average

lonremap  = 1          ; %Set 1 to remap LON from 360 -> -180 (Change [bbox] acc.)
lag       = [1:2]      ; %List of lags to average over/sum over
monwin    = 3          ; %Months considered (1 or 3)
mode      = 4          ; % (1) No Sig Testing, (2) SST testing (3) Flx testing (4) Both
ensavgf   = 1          ; %Set to 1 if you want to take the ensemble average first
lag1_elim = 0          ;% Set to 1 to eliminate rest of lags if lag 1 fails
insigNaN  = 0          ;% Set to 1 to change insignificant feedbacks to NaN
deg5      = 0          ;% Set to 1 to use smoothed data
ensorem   = 1          ;% Indicate if ENSO was removed

% Shouldn't change these, but I'll put these here for now...
monlist   = [1:12]     ; %List of months 
ensnum    = [1:42]     ; %List of ensemble members 

% -------------------
% SET SIGNIFICANCE TESTING OPTIONS
% -------------------
% T-Test parameters
p        = 0.05         ; % p-value that will be used
tails    = 2             ; % 1 or 2-tailed T-test
dof_man  = 82            ; % Manual DOF value


% Set Paths
projpath = '/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/';

% Add Paths

addpath('/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/01_Data'); % Path to common data
addpath('/Users/gliu/') % Path to scripts

if deg5 == 1 
    datpath  = [projpath,'01_Data/5deg/'];
    load('CESM1_LATLON_5deg.mat');
    lonsize = 72;
    latsize = 36;
else
    datpath  = [projpath,'01_Data/'];
    load('CESM1_LATLON.mat')
    lonsize = 288;
    latsize = 192;
end

% Set output location to same folder where data was stored
outpath = datpath;

% Path to data
addpath(datpath) 

% Load SST autocorrelation coefficients
load(strcat(datpath,['SST_rho_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat']))


startup
%% ------------------------------------------------------------------------
%  Script Start
% Determine effective DOF
if dof_man == 0
    if monwin == 1
        n_eff = 84;
    elseif monwin == 3
        n_eff = 82;
    end
% Manually set dof
else
    n_eff = dof_man;
end

% Calculate Correlation Threshold (fixed dof)
ptilde  = (1-p/tails);
critval = tinv(ptilde,n_eff);
corrthres = sqrt(1/((n_eff/critval^2)+1));

% Loop for each element of fluxes
for i = 1:length(fluxes)
    
    % Load variable and pull flux
    fluxtype = char(fluxes(i));
    fluxlab = upper(fluxtype);
    
    % Load damping variable
    load(strcat(datpath,[fluxlab,'_damping_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat']))
    v = damping;
    
    % Load correlation coefficients
    load(strcat(datpath,[fluxlab,'_rho_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat']))
    rflx = rho;
    
    % -----------------------------------
    % Masking
    % -----------------------------------
    % Make masks  
    if mode == 1
        % No mask (identity matrix)
        amask = ones(size(v));
        
    else
        
        % Create masks with threshold
        msst = rsst   > corrthres;
        mflx = rflx   > corrthres;

        
        % Lag 1-based elimination
        if lag1_elim == 1
            
            % Eliminate points that failed lag 1 for...
            % Autocorrelation
            msst(:,:,:,:,2) = msst(:,:,:,:,2) .* msst(:,:,:,:,1);
            msst(:,:,:,:,3) = msst(:,:,:,:,3) .* msst(:,:,:,:,1);
            
            % Cross Correlation
            mflx(:,:,:,:,2) = mflx(:,:,:,:,2) .* mflx(:,:,:,:,1);
            mflx(:,:,:,:,3) = mflx(:,:,:,:,3) .* mflx(:,:,:,:,1);            

        end
        
        % Selectively apply testing
        switch mode
            % Apply autocorrelation only
            case 2
                amask = msst;
            % Apply cross correlation only
            case 3
                amask = mflx;
            % Apply Both
            case 4
                amask = mflx .* msst;      
        end
        % End mask creation
    end
    
    % Turn 0 to NaNs if option is set
    if insigNaN == 1
        amask = double(amask);
        amask(amask==0) = NaN;
    end
    
    % Apply Mask
    v = v .* amask;
    
    % ----------------------------------------------
    % Further Adjustments
    % ----------------------------------------------
    
    % Check size and adjust variable for plotting
    chksize = size(v);
    refsize = [lonsize, latsize, 12, 42, 3];
    if ~isequal(chksize,refsize)
        fprintf('Permuting variable to match script requirements')
        corrsize = zeros(1,length(refsize));
            for ii = 1:length(refsize)
                corrsize(1,ii) = find(chksize==refsize(ii));
            end
        v = permute(v,corrsize);
    end
    
    % Remap LON to -180-180 if toggle is on
    if lonremap == 1
        if any(LON> 181)        
            [LON1,v] = lon360to180(LON,v);
        else
            LON1 = LON;
        end
    else
        LON1 = LON;
    end
    
    % Limit to specified lag [lon x lat x 12 x 42 x lag]
    ensavg = squeeze(v(:,:,:,:,lag));
    
    
    
    
    % Take ensemble average (Keep Monthly Data Intact)
    if ensavgf == 1
        % Take ensemble average first, then lag
        while length(size(ensavg)) > 3 
                ensavg = squeeze(nanmean(ensavg,4));
        end
    else
        fprintf("Ensavgf set to 0\n")
        ensavg = ensavg(:,:,:,:);
        ensavg = nanmean(ensavg,4);
    end
    
    % Save Variable (ex:
    % ensavg_NHFLXdamping_monwin1_sig010_dof082_mode4.mat)
    outname = [outpath,'ensavg_',fluxtype,'damping'...
        '_monwin',num2str(monwin),...
        '_sig',num2str(p*100,'%03d'),...
        '_dof',num2str(dof_man,'%03d'),...
        '_mode',num2str(mode),...
        '.mat'];
    save(outname,'ensavg','LON1','LAT','fluxtype')
    fprintf('\nSaved Ensavg for %s\n',fluxtype)
end