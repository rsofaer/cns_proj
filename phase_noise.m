% Input/Output parameters:
% num_samp desired number of output samples
% f0 reference frequency (must be in Hz.)
% dbc_per_hz power per hertz relative to carrier at ref. freq.
% num_taps number of filter taps in AR 1/f filter
% (optional; default = 100)
% 
% pn phase-modulated 1/f process
% theta 1/
function [pn, theta] = phase_noise(num_samp, f0, dbc_per_hz, num_taps)
  % Check input.

  if dbc_per_hz >= 0
  error('Power per Hz. must be negative.');
  elseif f0 <= 0
  error('Reference frequency must be positive.');
  end

  if nargin < 4
  num_taps = 100;
  end


  % Generate white noise. Apply gain for desired dBc/Hz. Warn user
  % if gain is too large (gain thresholds have been chosen somewhat
  % arbitrarily -- needs work).

  gain = sqrt(2*pi * f0 * 10^(dbc_per_hz/10));
  wn = gain * randn(1,num_samp);

  fprintf('Gain applied to white noise = %f.\n', gain);
  if gain >= 1
  fprintf('WARNING: Narrowband approximation no longer valid.\n');
  elseif gain >= .5
  fprintf('WARNING: Narrowband approximation on the verge of collapse.\n');
  end


  % Generate 1/f AR filter and apply to white noise to produce 1/f
  % noise.

  a = zeros(1,num_taps);
  a(1) = 1;
  for ii = 2:num_taps
  a(ii) = (ii - 2.5) * a(ii-1) / (ii-1);
  end
  theta = filter(1,a,wn);


  % Phase modulate.

  pn = exp(i*theta);

  return;
end
