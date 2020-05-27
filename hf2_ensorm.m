% -------------------------------------------------
% h2_ensorm
% -------------------------------------------------
%
% This script is split into two parts
% 1. Compute PC1 and PC2 from TS (product of hf1)
% 2. Regress PC1/PC2 onto a specific variable and remove the ENSO-related
% component of that variable
%
% Input
% 1. vars    - List of Fluxes/Variables {TS,'LHFLX','SHFLX','FSNS','FLNS'}
%              varrm [288 x 192 x 12 x 86]
% 2. ensolag - Lag in months between Variable and PC before the regression
% 3. monwin  - Window of months to consider (1 or 3)
% 4. pcrem   - PCs to remove [1,2,3]
% 5. EOFcorr - Set to 1 to correct the sign of the EOF based on predef
% boxes
%
% Output
% 1. outfile - Contains EOF1-3, PC1-3, VarPerc1-3
%              ex. ENSO_EOF_ens.mat
%              (1) EOF1-3       [42 x 137 x  12 x 107] lon|lat|mon|ens
%              (2) PC1-3        [ 1 x  86 x  12 x 107]  - |yrs|mon|ens
%              (3) VarPerc1-3   [ 1 x  12 x 107 x    ]  - |mon|ens
%
% 2. outname - Contains variable with ENSO Component removed, by ENSmember
%              ex. LHFLX_noenso_ENS070.mat;
%              (1) var_noenso   [12 x 84 x 288 x 192    ] mon|yrs|lat|lon
%              (2) all_ensocomp [12 x 84 x 288 x 192 x 2] mon|yrs|lat|lon|pcn
%
%  
%

%% User Input ------------
vars    = {'TS'};
ensolag = 1     ;% Lag between variable and ENSO Component Removed
monwin  = 3     ;
pcrem   = [1,2] ;% PCs to remove
EOFcorr = 1     ;% Option to correct signs using coords in Lat/Lon Boxes Section
deg5    = 1     ;% Use smoothed 5 degree data
printflip = 0   ;% Toggle to print warning when EOF was flipped due to criteria

% Path to data (monthly folder before the variable addition)
if deg5 == 1
    % Just add +5deg to everything
    datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/01_hf1out/5deg/';
    outpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/02_hf2out/5deg/';
    latlon  = '/home/glliu/01_Data/CESM1_LATLON_5deg.mat';
    outfile  =[outpath,'ENSO_EOF_ens_5deg.mat'];
else
    datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/01_hf1out/';
    outpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/02_hf2out/';
    latlon  = '/home/glliu/01_Data/CESM1_LATLON.mat';
    outfile  =[outpath,'ENSO_EOF_ens.mat'];
end

% Other sets
mnum  = [1:35,101:107];% Ensemble Members List


% Toggles
calc_PC = 0; % Calculate PC from TS data if set to 1
step3   = 1; % ENSO Removal

%% --------------------------------
%  Lat/Lon Box for EOF Calculation
%  --------------------------------
LAT1 = -20;
LAT2 = 20;
LON1 = 120;
LON2 = 290;

% EOF1 sign check Lat/Lon Box (Currently using nino3.4)
cLAT1 = -5;
cLAT2 = 5;
cLON1 = 190;
cLON2 = 240;

% EOF2 sign check Lat/Lon Box (Frankignoul et al 2011: 
% extent of positive contours offshore S. Am > 0.1)
dLAT1 = -5;
dLAT2 = 5;
dLON1 = 250;
dLON2 = 280;

% EOF3 sign check Lat/Lon Box (Frankignoul et al 2011"
% narrow band of positive values around equator > 0.1)
eLAT1 = -5;
eLAT2 = 5;
eLON1 = 200;
eLON2 = 260;

%% Script Start
allstart = datetime('now');
fprintf('Now running hf2_ensorm (%s)',allstart)

% Load startup
addpath('/home/glliu/')
startup

% Load Lat/lon
load(latlon)

% Set windowsize
if monwin == 3
    winsize = 1;
elseif monwin == 1
    winsize = 0;
end

%%  ------------------------
%   EOF/PC Calculation
%   ------------------------

if calc_PC == 1
    cpstart = datetime('now');
    fprintf('Calculating PC from SST (%s)',cpstart)
    
    % ---------------------
    % Preallocate Variables
    % ---------------------
    % Dimensions [1 x yr x mon x ensnum]
    PC1 = NaN(1,86,12,length(mnum)); PC2 = NaN(size(PC1)); PC3 = NaN(size(PC1));
    
    VarPerc1 = NaN(1,12,length(mnum)); VarPerc2 = size(VarPerc1); VarPerc3 = size(VarPerc1);
    
    EOF1=[];EOF2=[]; EOF3=[];
    
    % -------------------------------------
    % Load Lat/Lon and compute area weights
    % -------------------------------------
    [X,Y]=meshgrid(LON,LAT); 
    AREA=cos(Y/180*pi); % Note this is 192 x 288
    
    % -------------------------------------
    % Limit to area of interest
    % -------------------------------------
    % Find indices for area of interest
    X00=LON';
    Y00=LAT'; 
    ky=find((Y00>=LAT1)&(Y00<=LAT2));
    if LON1 > LON2 % Crossing Antimeridian
        kx=[find((X00>=LON1)&(X00<=180)),find((X00>=-180)&(X00<=LON2))];
        kx = sort(kx,2,'ascend');
    else
        kx=find((X00>=LON1)&(X00<=LON2));
    end
    
    % Restrict to area of interest
    AREA=AREA(ky,kx);
    Y0  = Y00(ky);   % Need size later for eof 
    X0  = X00(kx);   % Need size later for eof
    X   = X(ky,kx);
    Y   = Y(ky,kx);
    
    
    % -------------------------------------
    % Compute PC1/PC2 via yo_eof
    % -------------------------------------
    for n = 1:length(mnum)

        % Load SST data
        ensnum  = mnum(n);
        matname = [datpath,'/TS_ens',num2str(ensnum,'%03d'),'_proc.mat'];
        load(matname)
    
        % Find NaN Values
        varrm      = permute(varrm,[3,4,2,1]);   % move time dims to front (M X Y LAT X LON)
        varrm      = varrm(:,:,ky,kx);           % Take only selected region
        varr_tmean = mean(mean(varrm,1),2);      % Take time mean
        inan       = find(isfinite(varr_tmean)); % Find NaN Values
        varr_nonan = varrm(:,:,inan);            % Select non-NaN values

        % Preallocate variable to fill
        var_fill   = NaN(1,size(varrm,3)*size(varrm,4));
        
        % Calculate EOFs by month
        for m = 1:12

            % Select Month
            imon = m;
            var_mon = varr_nonan(m,:,:);
            Xin = X(:);
            Yin = Y(:);

            % Apply Area Weight
            WGT=sqrt(AREA(inan)'); % Apply NaN removal to WGT (inan is indexed 288x192)
            var_wgt = var_mon .* WGT(:,ones(1,size(var_mon,2)),:);

            % Run EOF
            fprintf('\n * Running EOF for mon %i ...',m);
            estart = datetime('now');
            [EV,PC,q]=yo_eof(squeeze(var_wgt)',3);

            % Normalize PC
            stdPC=std(PC,0,2); 
            PC=PC./stdPC(:,ones(1,size(PC,2)));
            clear var_wgt;
            elapsed = datetime('now') - estart;
            fprintf('  EOF done (%s)...',elapsed);

            % Regress PCs to calculate EOF pattern
            PC1a=PC(1,:); PC2a=PC(2,:); PC3a=PC(3,:);
            EV1a=PC1a*squeeze(var_mon)/length(PC1a);
            EV2a=PC2a*squeeze(var_mon)/length(PC2a);
            EV3a=PC3a*squeeze(var_mon)/length(PC3a);

            % Reindex EOFs into original lat/lon configurations pre NaN removal
            ev1a = var_fill; ev2a = var_fill; ev3a = var_fill;
            ev1a(:,inan) = EV1a;
            ev2a(:,inan) = EV2a;
            ev3a(:,inan) = EV3a;

            % Correct EOF signs if option is toggled on
            if EOFcorr == 1
                % Keep nino3.4 region positive for EOF1
                % Find indices
                kxy = find((Xin >= cLON1) & (Xin <= cLON2)...
                    & (Yin >= cLAT1) & (Yin <= cLAT2));
                % Change sign accordingly
                if nansum(ev1a(kxy)) < 0
                    if printflip == 1
                        fprintf('\nFlipping EOF 1 for mon %s since sum = %s',num2str(imon,'%02d'),num2str(nansum(ev1a(kxy))))
                    end
                    ev1a = ev1a*-1;
                    PC1a = PC1a*-1;
                end

                % Keep west coast of SA positive for EOF2
                kxy = find((Xin >= dLON1) & (Xin <= dLON2)...
                    & (Yin >= dLAT1) & (Yin <= dLAT2));
                if nansum(ev2a(kxy)) < 0
                    if printflip == 1
                        fprintf('\nFlipping EOF 2 for mon %s since sum = %s',num2str(imon,'%02d'),num2str(nansum(ev2a(kxy))))
                    end
                    ev2a = ev2a*-1;
                    PC2a = PC2a*-1;
                end

                % Keep narrow equatorial band positive for EOF4
                kxy = find((Xin >= eLON1) & (Xin <= eLON2)...
                    & (Yin >= eLAT1) & (Yin <= eLAT2));
                if nansum(ev3a(kxy)) < 0
                    if printflip == 1
                        fprintf('\n**Flipping EOF 3 for mon %s since sum = %s',num2str(imon,'%02d'),num2str(nansum(ev3a(kxy))))
                    end
                    ev3a = ev3a*-1;
                    PC3a = PC3a*-1;
                end
            end

            % store the results 
            VarPerc1(:,imon,n)=q(1);
            VarPerc2(:,imon,n)=q(2);
            VarPerc3(:,imon,n)=q(3);
            PC1(:,:,imon,n)=PC1a;
            PC2(:,:,imon,n)=PC2a;
            PC3(:,:,imon,n)=PC3a;
            EOF1(:,:,imon,n)=reshape(ev1a,length(Y0),length(X0));
            EOF2(:,:,imon,n)=reshape(ev2a,length(Y0),length(X0));
            EOF3(:,:,imon,n)=reshape(ev3a,length(Y0),length(X0));
            
            iloopend = datetime('now');
            elapsed = iloopend-cpstart;
            fprintf('\nCompleted ens %s at %s (%s)',num2str(ensnum,'%03d'),iloopend,elapsed)
        end
    end
    save(outfile,'EOF1','EOF2','EOF3','PC1','PC2','PC3','VarPerc1','VarPerc2','VarPerc3');
    cpend = datetime('now');
    elapsed = cpend-cpstart;
    fprintf('PC has been computed in %s (%s)',elapsed,cpend)
else
    load(outfile)
end

%%  ------------------------
%   ENSO Removal
%   ------------------------

for v = 1:length(vars)
    vstart = datetime('now');
    
    % Get Variable Name
    vname = char(vars(v));
    
    % Loop by each ensemble member
    for n = 1:length(mnum)
        
        % Load in variable
        ensnum  = mnum(n);
        matname = [datpath,'/',vname,'_ens',num2str(ensnum,'%03d'),'_proc.mat'];
        load(matname)
        
        % Permute to move time dim first
        varrm = permute(varrm,[3,4,1,2]);
        
        % Reduce yr dimension of variable to prepare for removal
        % Drop first and last month for 3-month window
        if winsize > 0
            var_noenso = varrm(:,2:end-1,:,:);
        else
            var_noenso = varrm(:,2:end,:,:);
        end
        finsize = size(var_noenso);
        
        % Preallocate enso component variable, with # of PCs to remove
        all_ensocomp = zeros([size(var_noenso),length(pcrem)]);

        % Loop by PC
        for pcn = 1:length(pcrem)   
            
            % Get PC (1 x yr x mon x ensnum)
            eval(['PC = PC',num2str(pcn),'(:,:,:,',num2str(n),');']);
            
            % Permute PC for input into laglead indexer
            % [1 x 86 x 12] -> [12 x 86]
            PC = squeeze(permute(PC,[3,2,1]));
            
            % Crossing Flags           
            tflag   = 0;         
            tshift  = 0;
            
            % If windowsize > 0, reduce leading end to allow for year
            % crossing acommodation
            if winsize > 0
                lshift  = -1;
            else
                lshift  = 0;
            end
            
            % Total time series length
            totyr  = size(varrm,2);
            
            % Loop by month
            for m = 1:12
                mstart = datetime('now');
                fprintf('\nRemoving Var:%s ENS:%s PC:%i Mon:%i...',vname,...
                    num2str(ensnum,'%03d'),pcn,m)
                % Get Lag Month
                lm = mod(m - ensolag,12);
                if lm == 0
                    lm = 12;
                end
                
                % When leading edge of the base variable's window crosses 
                % into next year, shift the end of the period forward
                % (When varidx -> ND|J)
                if m + winsize > 12
                    lshift = lshift + 1;
                end
                
                % When the trailing end of the lagged variable crosses into
                % the central period, shift the start of the period forward
                
                if tflag == 1
                    tshift = tshift + 1;
                    tflag = 0;
                end
                
                % Flag for tshift based on lag month and window size
                % (When pcidx -> D|JF)
                if winsize > 0
                    if lm - winsize == 0 && tflag == 0
                        tflag = 1;                 
                    end
                else
                    if lm == 12
                        tflag = 1;
                    end
                end
                
                % Produce time period
                base_ts = 1+tshift:totyr+lshift;
                
                % Index according to lag
                [varidx,basflag] = laglead_indexer(varrm,m,winsize,base_ts);
                [pcidx,lagflag] = laglead_indexer(PC,lm,winsize,base_ts);
                
                % Reduce yr dim (correct time periods) based on flag
                % If both months are NOT in the interior
                if lagflag+basflag ~= 0 || tshift == 0
                    % Take the interior period if month has crossed the initial
                    % threshold (trailing end of window no longer in the 
                    % previous year) (pcidx -> D|JF ; varidx -> |JFM)
                    if basflag == 0
                        varidx = varidx(:,2:end,:,:);
                    end
                
                    % Take the interior period once the leading edge of the
                    % window of the base variable has crossed into the next
                    % year (pcidx -> OND| ; varidx -> ND|J)
                    if lagflag == 0
                        pcidx = pcidx(:,1:end-1);
                    end
                end
                
                % Combine time dimensions to create 2D matrix
                % var2d: [252 x 55296], pc2d: [252 x 1]
                pc2d  = pcidx(:);
                vsize = size(varidx);
                var2d = reshape(varidx,prod(vsize(1:2)),prod(vsize(3:4)));
                
                % Regress field onto PC for beta (proportionality constant)
                % [ 1 x 55296 ]
                beta  = pc2d'*var2d./length(pc2d);
                
                % Repeat PC to match spdim of beta, and beta to PC
                % pc_rep [252 x 55296], beta_rep [252 x 55296]
                pc_rep    = pc2d(:,ones(1,size(beta,2)));
                beta_rep  = beta(ones(1,prod(vsize(1:2))),:);
                
                % Multiply beta by PC to obtain the ENSO Related component
                % enso_comp [252 x 55296]
                enso_comp = pc2d .* beta_rep;
                
                % Reshape and take mean along first dimension (month)
                enso_comp = reshape(enso_comp,[monwin,finsize(2:end)]); 
                enso_comp = nanmean(enso_comp,1);
                
                % Remove ENSO component and store into variable
                var_noenso(m,:,:,:) = var_noenso(m,:,:,:) - enso_comp;
                all_ensocomp(m,:,:,:,pcn) = all_ensocomp(m,:,:,:,pcn) + enso_comp;
                
                elapsed = datetime('now') - mstart;
                fprintf('Complete! (%s | %s)',elapsed,datetime('now')-allstart)
            end               
        end
        
        % Save variables after pc and month loops
        sstart  = datetime('now');
        outname = [outpath,vname,'_noenso_ENS',num2str(ensnum,'%03d'),'.mat'];
        save(outname,'all_ensocomp','var_noenso')
        ssend    = datetime('now');
        fprintf('\nSaved variable for ENS %s! (%s | %s)',...
            num2str(ensnum,'%03d'),ssend-sstart,ssend-allstart)
        clear var_noenso
        clear all_ensocomp
    end
    vend = datetime('now');
    elapsed = vend-vstart;
    fprintf('Completed %s in %s (%s)',vname,elapsed,vend-allstart)
    
end

                
                
                

                
                
                


        




