% meanback_infer.m
% Inputs (example):
% tracks: table or struct array with fields .id, .frame, .x, .y
% dt: frame time (s)
% T: temperature (K)
% sigma_loc: localization std dev (m)
% min_len: minimum frames per track (e.g. 20)
% -------------------------------------------------------------------------
function results = meanback_infer(trajectories, dt, T)

kB = 1.380649e-23;

ntraj = size(trajectories,2);
res = struct();
res.ntraj = ntraj;
res.k_per_dim = nan(ntraj,1);
res.k_radial = nan(ntraj,1);
res.tau_rel = nan(ntraj,1);
res.D_short = nan(ntraj,1);

% compute per particle quantities
for i=1:ntraj
    tr = trajectories(i).positions;
    % detrend position by subtracting mean (removes static offset)
    x = tr(:,1);
    y = tr(:,2);
    varx = var(x,1) ; % use population var (1/N)
    vary = var(y,1) ;
    if varx <= 0 || vary <= 0
        % skip or set to NaN
        continue;
    end
    % k per dimension
    kx = kB*T / varx;
    ky = kB*T / vary;
    res.k_per_dim(i) = 0.5*(kx+ky); % average of x and y
    res.k_radial(i) = kB*T / (varx + vary); % effective radial stiffness
    
    % compute autocorrelation C(tau) upto maxlag
    N = size(tr,1);
    maxlag = min(floor(N/2), 50);
    C = zeros(maxlag+1,1); cnt=zeros(maxlag+1,1);
    % Handle tau=0 case separately
    C(1) = mean(x.^2 + y.^2);
    cnt(1) = length(x);
    
    for tau = 1:maxlag
        % Check that the end index is valid
        if (length(x) - tau) >= 1
            X1 = x(1:end-tau); 
            X2 = x(1+tau:end);
            Y1 = y(1:end-tau); 
            Y2 = y(1+tau:end);
            C(tau+1) = mean(X1.*X2 + Y1.*Y2); 
            cnt(tau+1) = length(X1);
        else
            % Stop the loop if the array slice becomes invalid
            break;
        end
    end

    steps = diff([x y]);                % (N-1)x2 displacements
    maxlag_mbr = min(floor((size(steps,1))/2), 100);
    Cstep = zeros(maxlag_mbr,1);

    for tau = 1:maxlag_mbr
        % correlation of steps separated by tau frames
        s1 = steps(1:end-tau,:);
        s2 = steps(1+tau:end,:);
        Cstep(tau) = mean(sum(s1.*s2,2));  % dot product per step
    end

    % normalise by zero-lag (variance of steps)
    C0step = mean(sum(steps.^2,2));
    Rstep = Cstep / C0step;  % normalised MBR curve

    % store for later plotting
    res.MBR{i} = Rstep;                 % cell array, one per trajectory
    res.MBR_t{i} = (1:maxlag_mbr)'*dt;  % time axis per trajectory

    % normalize and fit exponential
    tvec = (0:maxlag)' * dt;
    C0 = C(1);
    if C0 <= 0
        continue;
    end
    R = C / C0;
    % fit R(t) = exp(-t/tau_rel) (ignore t=0 for fit)
    try
        ft = fit(tvec(2:end), R(2:end), 'exp1'); % exp(a*x) where a = -1/tau
        a = ft.a; b = ft.b; % exp(a*t)+b ; but exp1 = a*exp(b*x)? depends on fit
        % safer to do custom fit:
        p = polyfit(tvec(2:end), log(R(2:end)),1);
        slope = p(1);
        tau_rel = -1/slope; % because log(R) ~ -t/tau
        res.tau_rel(i) = tau_rel;
    catch
        res.tau_rel(i) = NaN;
    end
    
    % short time D from MSD slope using first few lags
    % compute MSD for lags 1:5
    maxlag_msd = min(5, floor(N/2));
    msd = zeros(maxlag_msd,1);
    lags = (1:maxlag_msd)';
    for m=1:maxlag_msd
        diffs = tr(1+m:end,:) - tr(1:end-m,:);
        msd(m) = mean(sum(diffs.^2,2));
    end
    % subtract localization offset (4*sigma^2)
    if any(msd<=0)
        res.D_short(i) = NaN;
    else
        % linear fit msd = 4 D t
        tt = lags*dt;
        p = polyfit(tt, msd, 1);
        slope = p(1); % slope = 4D
        res.D_short(i) = slope / 4;

    end
end

% compute k, gamma, D from averages or per particle if available
% use per-particle tau and k to compute gamma and D_est
res.gamma = res.tau_rel .* res.k_per_dim;
res.D_from_gamma = (kB*T) ./ res.gamma;
res.D_trap = (kB*T)./ res.gamma;


% summary statistics
results = res;
end