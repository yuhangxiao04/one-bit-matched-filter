%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMCW (Phase-Modulated Continuous-Wave) radar simulation
%
%   • Waveform    : 1-GHz BPSK code, 1024 chips → 1.024 µs per code
%   • Array       : 3-Tx (time division) × 4-Rx  → 12-element virtual ULA
%   • Processing  : one-bit matched filter → Doppler FFT →   MUSIC
%   • Track model : analytic (x,y) curve; **one plotted dot every dt = 0.2 s**
%                   The radar itself keeps working continuously:
%                     – each coherent CPI = N_codes/PRF = 6.4 ms
%                     – ~31 CPIs occur inside every 0.2 s “visual frame”
%                   (Set dt = N_codes/PRF to emulate a fully duty-cycled
%                    PMCW radar without touching any other code line.)
%   • Power law   : amplitude ∝ 1/R² → power ∝ 1/R⁴  (radar equation)
%   • Quantizer   : complex 1-bit (sign on I & Q)
%   • Noise       : AWGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all


%% -------------------- Radar / waveform parameters --------------------
c            = 3e8;                               % speed of light (m/s)
f0           = 77e9;                              % carrier frequency (Hz)
Rc           = 1e9;                               % chip-rate / code rate (chips per second)
codeLen      = 1024;                              % number of chips per code word
T_code       = codeLen/Rc;                        % code duration = 1.024 µs
PRF          = 1e4;                               % super-PRI (code-burst repetition frequency, after 3 Tx)
N_codes      = 64;                                % number of code words per coherent CPI
radar_pos    = [0,0];                             % radar origin (x,y) in metres
 
fs           = Rc;                                % sampling rate: one sample per chip
chips        = 0:codeLen-1;                       % fast-time index within one code
codeSeq      = 2*randi([0 1],1,codeLen)-1;        % random BPSK (+1 / -1) spreading code

lambda       = c/f0;
range_res    = c/(2*Rc);

%% --------------------  Track generator: produces N points (x,y) and matching velocity -----------------

dt         = 0.1;                    % time between exported track points (visual sampling);
                                     % radar actually runs ~31 CPI = 6.4 ms each inside this interval
T_total    = 6.4;                    % total duration (s)
t_vec      = 0:dt:T_total;           % time vector
Npts       = numel(t_vec);

% --- Shape design ------------------------------------------------------
% The track rises in y while bending gently right, then back left – no loops
y_traj = 20 + 4.5*t_vec + 0.5*t_vec.^2;                 % always increasing
x_traj = 10 + 6*sin(0.5*t_vec) + 4.5*t_vec.*exp(-0.6*t_vec);

% --- Velocity (first-order finite difference) --------------------------
vx = [diff(x_traj)/dt  0];          % append a 0 so vectors are same length
vy = [diff(y_traj)/dt  0];
spd = hypot(vx,vy);                 % speed magnitude (optional)

% --- Package outputs ---------------------------------------------------
pos = [x_traj.'  , y_traj.'];     % N×2   [x  y]
vel = [vx.'      , vy.'    ];     % N×2   [vx vy]

R_ref        = norm(pos(1,:) - radar_pos);  % Reference range: distance from radar to the first trajectory point               
SNR_dB       = 0;                           % Desired receive-SNR (dB) when the target is at the reference range R_ref

%% -------------------- Antenna geometry –– 3 Tx × 4 Rx → 12-element virtual ULA -----------------
Ntx = 3;              % transmit antennas (time-division)
Nrx = 4;              % simultaneous receive antennas
M   = Ntx*Nrx;        % virtual array length (should stay 12)

d_rx = lambda/2;                    % physical Rx spacing
rx_pos = (0:Nrx-1)*d_rx;            % Rx x-coordinates

d_tx = Nrx*d_rx;                    % Tx spacing → contiguous virtual grid
tx_pos = (0:Ntx-1)*d_tx;            % Tx x-coordinates

est_x = zeros(1,length(t_vec));
est_y = zeros(1,length(t_vec));

%% -------------------- Target scenario (single target) -----------------
for kk=1:length(t_vec)
    
true_pos         = pos(kk,:);
true_vel         = vel(kk,:);

% Complex reflection coefficient:
%   • |reflectivity| is set so that the received SNR equals SNR_dB at R =  R_ref, with random phase in each scan
reflectivity     = 10^(SNR_dB/20)*exp(1j*2*pi*rand); 

%% -------------------- Pre-allocate 4-D cube  (fast × slow × Rx × Tx)
rx_cube = zeros(codeLen, N_codes, Nrx, Ntx);

%% -------------------- Simulation loop over repeated codes -------------
for n = 1:N_codes
    t0 = (n-1)/PRF;                         % super-PRI origin
    true_pos = true_pos + true_vel/PRF;     % linear motion

    R                = norm(true_pos - radar_pos);
    theta_true       = atan2(true_pos(1),true_pos(2));

    radial_vel = dot(true_vel,(true_pos-radar_pos)/R);
    fd = 2*radial_vel/lambda;               % Doppler (Hz)

    % --- free-space amplitude for this CPI (amplitude ∝ 1/R²) --------------
    amp_n = reflectivity * (R_ref/R)^2;        

    % iterate over Tx bursts (time-division)
    for k_tx = 1:Ntx
        t_tx = (k_tx-1)*T_code;             % Tx-k firing delay

        for m_rx = 1:Nrx
            %% ----- geometry for this Tx/Rx pair -----------------------
            virt_x   = tx_pos(k_tx) + rx_pos(m_rx);
            R_m      = R + virt_x*sin(theta_true);
            tau_m    = 2*R_m/c;

            %% ----- ideal baseband echo  (code + Doppler) --------------
            t_fast = chips/fs;
            sig = amp_n * codeSeq .* exp(1j*2*pi*fd*(t_fast + t0 + t_tx));

            %% ----- integer & fractional chip delay --------------------
            shift_samples = round(tau_m*fs);
            frac_delay    = tau_m - shift_samples/fs;

            sig = circshift(sig,[0, shift_samples]);
            kf  = 0:codeLen-1;
            SIG = fft(sig);
            SIG = SIG .* exp(-1j*2*pi*kf*frac_delay/codeLen);
            sig = ifft(SIG);

            %% ----- steering phase BEFORE noise & 1-bit  ---------------
            geom_phase = exp(1j*2*pi*virt_x*sin(theta_true)/lambda);
            sig = sig .* geom_phase;

            %% ----- AWGN + 1-bit quantization --------------------------
            sig = sig + sqrt(0.5)*(randn(size(sig))+1j*randn(size(sig)));
            sig = csign(sig);

            %% ----- store into 4-D cube --------------------------------
            rx_cube(:,n,m_rx,k_tx) = sig.';
        end
    end
end

%% -------------------- Range processing: one-bit matched filter ----------------
CODE_FFT  = conj(fft(codeSeq));
range_mat = zeros(codeLen, N_codes, Nrx, Ntx);
for k_tx = 1:Ntx
  for m_rx = 1:Nrx
    RANGE = ifft( fft(rx_cube(:,:,m_rx,k_tx),[],1) .* CODE_FFT.', [], 1 );
    range_mat(:,:,m_rx,k_tx) = RANGE;
  end
end
range_axis = (0:codeLen-1)*range_res;

%% -------------------- Doppler processing ------------------------------
N_fft_doppler = N_codes;
dopp_cube = fftshift( fft(range_mat, N_fft_doppler, 2), 2 );

doppler_axis  = ( -N_fft_doppler/2 : N_fft_doppler/2-1 ) * PRF/N_fft_doppler;
velocity_axis = doppler_axis * lambda/2;

%% -------------------- Coarse peak search ------------------------------
abs_map = squeeze(sum(sum(abs(dopp_cube).^2,3),4));   % integrate Rx+Tx
[~, idx_max]   = max(abs_map(:));
[r_idx, d_idx] = ind2sub(size(abs_map), idx_max);

est_range    = range_axis(r_idx);
fd_hat       = doppler_axis(d_idx);
est_velocity = fd_hat*lambda/2;
fprintf('Estimated Range   : %.2f m\n', est_range);
fprintf('Estimated Velocity: %.2f m/s\n', est_velocity);

%% -------------------- Build deskewed 12-element virtual snapshot -----------------
snapshot = zeros(M,1);
cnt = 1;
for k_tx = 1:Ntx
    deskew = exp(-1j*2*pi*fd_hat*(k_tx-1)*T_code);  % Doppler deskew
    for m_rx = 1:Nrx
        snapshot(cnt) = dopp_cube(r_idx,d_idx,m_rx,k_tx) * deskew;
        cnt = cnt + 1;
    end
end

%% -------------------- Spatial-smoothing MUSIC -------------------------
P = 4;  L = M-P+1;  X = zeros(P,L);
for k = 1:L, X(:,k) = snapshot(k:P+k-1); end
Rx = (X*X')/L;

[Evec, Eval] = eig(Rx);
[~, sortID]  = sort(diag(Eval), 'descend');
Un           = Evec(:,sortID(2:end));

angle_grid   = -90:0.025:90; 
MUSIC_spec   = zeros(size(angle_grid));
for ii = 1:length(angle_grid)
    theta = angle_grid(ii)*pi/180;
    steering = exp(1j*(0:P-1).' * pi * sin(theta));
    MUSIC_spec(ii) = 1./abs( steering'*(Un*Un')*steering );
end
[~, ang_idx]  = max(MUSIC_spec);
angle_est_deg = angle_grid(ang_idx);

fprintf('Estimated Angle  : %.2f deg\n', angle_est_deg);
fprintf('True Angle  : %.2f deg\n', theta_true*180/pi);
%% -------------------- Cartesian position ------------------------------
angle_est_rad = angle_est_deg*pi/180;
est_x(kk) = est_range * sin(angle_est_rad);
est_y(kk) = est_range * cos(angle_est_rad);


end

%% -------------------- Sanity plot ------------------------------------------------
figure;  plot(pos(:,1),pos(:,2),'r.');hold on
plot(est_x(:),est_y(:),'b.');
axis equal
xlabel('Horizontal Distance (m)'); ylabel('Vertical Distance (m)'); grid on
legend('Ground truth','Estimated');
title('Simulated Trajectory');



%% -------------------- Utility : complex sign (1-bit) ------------------
function y = csign(x)
    y = sign(real(x)) + 1j*sign(imag(x));
end

