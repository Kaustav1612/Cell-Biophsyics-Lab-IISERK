function [omega, Gstar, Gp, Gpp] = GSER(MSD,nframes,new_stats,max_paricles,t, MSD, a, T)

% Computes complex modulus from MSD(t)
%
% Inputs:
%   t    : time vector (s)
%   MSD  : mean squared displacement at time t (m^2)
%   a    : probe radius (m)
%   T    : temperature (K)
%
% Outputs:
%   omega : angular frequency vector (rad/s)
%   Gstar : complex modulus G*(omega) (Pa)
%   Gp    : storage modulus G'(omega) (Pa)
%   Gpp   : loss modulus G''(omega) (Pa)

kB = 1.380649e-23; % Boltzmann constant J/K
for k =1 :size(new_stats)
   a = [a,new_stats]
end
% ---- step 1: differentiate MSD to get D(t)
dMSDdt = gradient(MSD, t); % derivative
D_t = dMSDdt / 6;          % time-dependent diffusion coefficient (m^2/s)

% ---- step 2: frequency vector for FFT
N = length(t);
dt = mean(diff(t));
fs = 1/dt;
omega = (0:N-1)' * 2*pi*fs/N; % angular frequency vector

% ---- step 3: Fourier transform of derivative of MSD
FT_dMSD = fft(dMSDdt); 

% ---- step 4: complex diffusion coefficient D*(omega)
D_star = FT_dMSD ./ (6); % (m^2/s) in frequency domain

% ---- step 5: compute G*(omega)
Gstar = (kB*T) ./ (6*pi*a*D_star); % Pa

% Extract storage (real) and loss (imag) moduli
Gp = real(Gstar);
Gpp = imag(Gstar);

% optional: only keep first half of spectrum (positive frequencies)
n_half = floor(N/2);
omega = omega(1:n_half);
Gp = Gp(1:n_half);
Gpp = Gpp(1:n_half);
Gstar = Gstar(1:n_half);

end

