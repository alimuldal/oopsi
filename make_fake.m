function [F, V, P] = make_fake(V, P, seed)

% algorithm variables
% ----------------------------
if nargin < 1,              V   = struct;       end
if ~isfield(V,'Ncells'),    V.Ncells = 1;       end     % # of cells in image
if ~isfield(V,'T'),         V.T = 1000;         end     % # of time steps
if ~isfield(V,'dt'),        V.dt = 0.02;        end     % frame duration
if ~isfield(V,'w'),         V.w = 1;            end     % ROI width (px)
if ~isfield(V,'h'),         V.h = 1;            end     % ROI height (px)
if ~isfield(V,'Npixels'),   V.Npixels = V.w*V.h;end     % # of pixels in ROI

if ~isfield(V,'fast_plot'), V.fast_plot = 1;    end     % plotting on
if ~isfield(V,'fast_iter_max'), V.fast_iter_max = 10;    end     % plotting on

% neuron model parameters
% ----------------------------
if nargin < 2,              P = struct;         end

% scaling parameter / spatial mask
if ~isfield(P, 'a'),
    if V.Npixels > 1,
        if ~isfield(P, 'mask_sigma'),   P.mask_sigma = 10;  end
        % 2D Gaussian mask
        [x, y] = meshgrid(-V.h/2:1:V.h/2, -V.w/2:1:V.w/2);
        xs = x^2;
        ys = y^2;
        two_ss = 2*V.mask_sigma^2;
        mask = exp(-1*((xs + ys)/two_ss));
        mask = mask(:);
        mask = mask ./ sum(mask);                % normalize to sum to 1
        P.a = mask;
    else
        P.a = 1;
    end
end

% baseline
if ~isfield(P, 'b'),
    if V.Npixels > 1,
        % standard deviation of background
        if ~isfield(P, 'bg_intensity'), P.bg_intensity = 0.1;   end
        P.b = randn(V.Npixels, 1) .* P.bg_intensity;
    else
        P.b = 0.1;
    end
end

% standard deviation of noise
if ~isfield(P, 'sig')
    if V.Npixels > 1
        P.sig = 0.001;
    else
        P.sig = 0.1;
    end
end

% decay parameter
if ~isfield(P, 'gam')
    tau = 1;        % decay time constant
    P.gam = exp(-V.dt / tau) * ones(V.Ncells, 1);
end

% mean firing rate (Hz)
if ~isfield(P, 'lam'),  P.lam = ones(V.Ncells, 1);	end

% set the seed for the random number generator
if nargin == 3, rng(seed); end


% generate fake data
% ----------------------------

V.n = zeros(V.T, V.Ncells);
V.C = zeros(V.T, V.Ncells);
for ii = 1:V.Ncells
    % poisson spikes
    V.n(:, ii) = poissrnd(P.lam(ii) * V.dt, [V.T, 1]);

    % simulated calcium
    V.C(:, ii) = filter(1, [1, -P.gam(ii)], V.n(:, ii));
end

% scaling
F = P.a * V.C';

% baseline
F = bsxfun(@plus, F, P.b);

% noise
epsilon = randn(V.Npixels, V.T) * P.sig;
F = bsxfun(@plus, F, epsilon);