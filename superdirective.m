%% FUNCTION: superdirective_CM20
function [H_opt, Fvv] = superdirective(azimuth_angle, elevation_angle, Nfft, fs, M, mic_radius, c, reg_alpha)
% SUPERDIRECTIVE_CM20 computes superdirective weights H_opt and diffuse matrix Fvv.
% Signature required: [H_opt, Fvv] = superdirective_CM20(azimuth_angle, elevation_angle, Nfft)
% Optional args: fs (default 16000), M (default 19), mic_radius (default 0.042), c (default 340),
%                reg_alpha (default 1e-3 .. larger for stronger regularization)
%
% H_opt: M x K complex weights (K = Nfft/2 + 1)
% Fvv:  K x M x M diffuse coherence (unregularized)

% set defaults if missing
if nargin < 4 || isempty(fs), fs = 16000; end
if nargin < 5 || isempty(M), M = 19; end
if nargin < 6 || isempty(mic_radius), mic_radius = 0.042; end
if nargin < 7 || isempty(c), c = 340; end
if nargin < 8 || isempty(reg_alpha), reg_alpha = 1e-3; end

K = Nfft/2 + 1;
F = linspace(0, fs/2, K);    % frequency vector (Hz)
F(1) = max(F(1), (1e-1)*fs/Nfft);  % avoid exact zero in denominator

% circular microphone angles gamma (rad)
gamma = (0:(M-1))' * (2*pi / M);   % Mx1

% radii vector R (M x 1)
R = ones(M,1) * mic_radius;

% precompute pairwise Euclidean distances between mics on the circle
theta = gamma;  % same naming as original
dij = zeros(M, M);
for i = 1:M
    for j = 1:M
        dij(i,j) = sqrt( (R(i)*cos(theta(i)) - R(j)*cos(theta(j)))^2 + (R(i)*sin(theta(i)) - R(j)*sin(theta(j)))^2 );
    end
end
% allocate outputs
H_opt = zeros(M, K);          % M x K
Fvv = zeros(K, M, M);         % K x M x M

% convert angles to radians
fai0 = deg2rad(azimuth_angle);
seta0 = deg2rad(elevation_angle);
% loop over frequencies and compute Fvv(k,:,:), steering w0 and superdirective H0
for k = 1:K
    fk = F(k);
    % steering vector (M x 1)
    % w0 = exp(-1j * 2*pi * R * F(k) / c * sin(seta0) .* cos(fai0 - gamma));
    phase = -1j * 2*pi * fk / c .* ( R * sin(seta0) .* cos(fai0 - gamma) ); % Mx1 (complex)
    w0 = exp(phase);
    w0 = exp(-1j*2*pi * F(k)/c .* ( R .* cos(seta0) .* cos(fai0 - gamma) ) );
     % compute diffuse coherence matrix Fk (M x M)
    Fk = zeros(M, M);
%     for i = 1:M
%         for j = 1:M
%             if i == j
%                 Fk(i,j) = 1;
%             else
%                 x = 2 * pi * fk * dij(i,j) / c;
%                 if abs(x) < 1e-8
%                     Fk(i,j) = 1.0;
%                 else
%                     Fk(i,j) = sin(x) / x;   % your original formulation (sin(x)/x)
%                 end
%             end
%         end
%     end
%% diffuse noise field MSC
    f0 = F(k);
    % f0(1) = 1e-8;
    Fvv = zeros(Nfft/2+1,M,M);
    for i = 1:M
        for j = 1:M   
            if i == j
                Fk(i,j) = 1;
            else
                % dij = 2*R*sin(pi*abs(i-j)/M);
                % 127麦同心圆阵双麦间距(欧式距离)
                dij=sqrt((R(i)*cos(theta(i))-R(j)*cos(theta(j))).^2+(R(i)*sin(theta(i))-R(j)*sin(theta(j))).^2); 
                Fk(i,j) = sin(2*pi*f0*dij*1.0/c)./(2*pi*f0*dij*1.0/c);
            end
        end
    end

      % store unregularized Fvv
    Fvv(k,:,:) = Fk;

    %Fk = Fvv;

    % regularize Fk for numerical stability (proportional to trace)
    traceFk = real(trace(Fk));
    reg_val = reg_alpha * max(traceFk / M, 1e-8);
    Fk_reg = Fk + reg_alpha * eye(M);

    % compute superdirective weights H0 = inv(Fk_reg) * w0 / (w0' * inv(Fk_reg) * w0)
    invFk_w0 = Fk_reg \ w0;               % solve to avoid explicit inverse
    denom = (w0' * invFk_w0);
    if abs(denom) < 1e-12
        H0 = invFk_w0;   % fallback, unnormalized
    else
        H0 = invFk_w0 / denom;
    end
    H0 = invFk_w0 / denom;

    % store column (M x 1)
    H_opt(:, k) = H0;
    %H_opt(:,k) = w0/M;
end

end  % end function superdirective

    
