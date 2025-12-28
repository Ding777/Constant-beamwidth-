%% superdirective_pipeline.m

close all;
clc;
clear all

% ---------------- USER PARAMETERS (tweak if needed) ----------------
angle_az   = -138 ;        % degrees (target azimuth)
angle_el   = -18;        % degrees (target elevation)
Nfft       = 320;       % FFT length (must be same in frame loop)
fs         = 16000;     % sampling rate
c          = 340;       % speed of sound
M          = 8;        % number of mics
mic_radius = 0.042;     % circular array radius (m)
reg_alpha  = 0.1;      % regularization factor inside function (tune: 1e-3..1e-1)
%audioFolder = 'E:\All_mvdr\DNN_DATA\val5\small_office_rt0.15\ex_14d99908_mix.wav';
audioFolder = 'F:\test0_enhanced_out_folder_data_train_alldata_one2_targetone1_audio2_all550okok15data4oka_enhancedok3_diffusion_outlier55_audio_spatial\scene_000005_mix_multich.wav';
[xr, fs] = audioread(audioFolder);
% ---------------- compute superdirective weights ------------------
% returns Wopt_all: M x K (K = Nfft/2+1) and Fvv: K x M x M
[Wopt_all, Fvv] = superdirective(angle_az, angle_el, Nfft, fs, M, mic_radius, c, reg_alpha);

% provide a DS-equivalent copy if needed
Wopt_ds = Wopt_all;
% ---------------- prepare signals (these must exist in workspace) ----------------
% xr, noise_only, speech_clean : expected shape M x Nsamps
% If they come as transposed compared to your previous runs, adjust accordingly.
xr = xr';            % you used xr=xr' originally; preserve that behavior
% % noise_only = noise_only';
% % speech_clean = speech_clean';

% ---------------- parameters for frame-based processing (keep your logic) ---
N = length(xr(1,:));             % total samples per channel
L = Nfft/2;                      % hop (you used this)
G = fix((N - 2*L + L) / L);      % number of blocks (same as original)
y1 = zeros(N, 1);                % BF output (time-domain)
y1_ds = zeros(N, 1);
buffer1 = zeros(1, L);
buffer1_ds = zeros(1, L);
buffer2 = zeros(M, L);
% window: you used a sinusoidal window
win = sin((0.5 : Nfft-0.5) / Nfft * pi);   % row vector

% initialization of overlap history
xo = xr(:, 1:L);


% STFT buffers
K_m = Nfft/2 + 1;
Y = zeros(1, K_m);
Y_ds = zeros(1, K_m);
Y1 = zeros(M, K_m);
% energy arrays
E_Y = zeros(K_m, 1);
E_Yds = zeros(K_m, 1);

% decision-directed state (preserve your params; cleaned)
alpha_ps = 0.8; alpha_xi = 0.9;
floor_xi = 1e-6;
G_prev = ones(K_m, 1);
gamma_prev = ones(K_m, 1);
Ps_smooth = 1e-6 * ones(K_m, 1);
Pn_smooth = 1e-6 * ones(K_m, 1);
% main frame loop - preserved structure (cleaned variable names & indexing)
tic
for i = 1:G
    pos = L * i + 1;
    xk = xr(:, pos:pos+L-1);


    % overlap-add input block (concatenate previous L and current L)
    x = [xo, xk] .* win;            % size: M x Nfft (transposed before fft)


    % update history
    xo = xk; 
    % FFT along time for each mic: compute X as M x K matrix
    X0 = fft(x', Nfft);   % Nfft x M

    X = X0.';   % M x Nfft


% %     % mask & PSD tensors expected to be provided from caller (mask(:,i), PhiS_time(:,:,:,i), PhiN_time(:,:,:,i))
% %     mask1 = mask(:, i);                      % your mask (if used)
% %     PhiS_time1 = PhiS_time(:, :, :, i);      % M x M x K
% %     PhiN_time1 = PhiN_time(:, :, :, i);

    % build frequency-domain outputs Y and Y_ds (bins 1..K_m)
    xxF_Xf = zeros(M, K_m);   % store narrowband vectors across bins
    for n = 1:K_m
        Xf = X(:, n);                 % M x 1 narrowband across mics
        xxF_Xf(:, n) = Xf;

 % beamformer outputs per freq
        Y(n) = (xxF_Xf(:, n).' * (Wopt_all(:, n)));
        Y_ds(n) = (xxF_Xf(:, n).' * (Wopt_ds(:, n)));

        % instantaneous energy/power (for later)
        E_Y(n) = abs(Y(n));
        E_Yds(n) = abs(Y_ds(n)).^2;
    end
     % normalize for visualization if needed
    E_Y_norm = E_Y / max(E_Y + eps);

    % elementwise multiplication used later: produce time-domain M-channel beamformer partials
    Y1 = xxF_Xf .* (Wopt_all);   % M x K_m, because Wopt_all is MxK

%     % now compute per-bin projected powers using PhiN_time1 and PhiS_time1
%     for n = 1:K_m
%         % compute beamformer-projected noise and signal power
%         w = Wopt_all(:, n);
%         P_n = abs( w' * squeeze(PhiN_time1(:, :, n)) * w ) + eps;
%         P_s = abs( w' * squeeze(PhiS_time1(:, :, n)) * w ) + eps;
%         % you can use P_n, P_s for gain computations
%     end
 % reconstruct full-spectrum symmetric bins for IFFT
    Y(Nfft:-1:K_m+1) = conj(Y(2:K_m-1));
    Y_ds(Nfft:-1:K_m+1) = conj(Y_ds(2:K_m-1));
    Y1(:, Nfft:-1:K_m+1) = conj(Y1(:, 2:K_m-1));

    % IFFT and overlap-add: convert back to time-domain blocks
    temp1 = real(ifft(Y)) .* win;
    temp1_ds = real(ifft(Y_ds)) .* win;
    temp2 = real(ifft(Y1.')) .* win.';   % yields Nfft x M, transpose later
    temp2 = temp2.';
     % overlap-add into output buffers (exact same positions as you had)
    y1(pos-L:pos-1) = buffer1 + temp1(1:L);
    y1_ds(pos-L:pos-1) = buffer1_ds + temp1_ds(1:L);
    y2(:, pos-L:pos-1) = buffer2 + temp2(:, 1:L);

    % update buffers for next overlap
    buffer1 = temp1(L+1:2*L);
    buffer1_ds = temp1_ds(L+1:2*L);
    buffer2 = temp2(:, L+1:2*L);
end
toc
% final outputs (same variable names you used)
y = y1;
y01 = y1_ds;
r01 = y2;
 % Prepare output directory
    outDir = ['F:\Gridd_superdirectiveok1'];
%     outDir = ['Real_audio_output_fixedAZ_50deg_FIR5'];
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Save beamformed output named by original elevation
%     outName = sprintf('%d.wav', el_in);
    outName = sprintf('%d.wav', angle_az);
    audiowrite(fullfile(outDir, outName), y, fs);


