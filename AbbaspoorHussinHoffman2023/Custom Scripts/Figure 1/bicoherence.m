function [bic,waxis] = bicoherence(y, varargin)
%BICOHER - Direct (FD) method for estimating bicoherence
%	[bic,waxis] = bicoherence (y,  varagin)
%	y     - data vector or time-series
%	nfft - fft length [default = power of two > segsamp]
%	       actual size used is power of two greater than 'nsamp'
%	wind - specifies the time-domain window to be applied to each
%	       data segment; should be of length 'segsamp' (see below);
%		otherwise, the default Hanning window is used.
%	segsamp - samples per segment [default: such that we have 8 segments]
%	        - if x is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap, allowed range [0,99]. [default = 50];
%	        - if x is a matrix, overlap is set to 0.
%	bic     - estimated bicoherence: an nfft x nfft array, with origin
%	          at the center, and axes pointing down and to the right.
%	waxis   - vector of frequencies associated with the rows and columns
%	          of bic;  sampling frequency is assumed to be 1.
% This function is slightly changed from HOSA biocoher version.
% Saman Abbaspoor, 04/25/2021

%%
p = inputParser;
addParameter(p,'nfft',128,@isnumeric)
addParameter(p,'wind',0,@isnumeric)         
addParameter(p,'nsamp',0,@isnumeric)    
addParameter(p,'overlap',50,@isnumeric)           
parse(p,varargin{:})


% Assign Parameters
nfft       = p.Results.nfft;
wind       = p.Results.wind;
nsamp      = p.Results.nsamp;
overlap    = p.Results.overlap;

% --------------------- parameter checks -----------------------------

[ly, nrecs] = size(y);
if ly == 1
    y = y';
    [ly, nrecs] = size(y);
end

overlap = max(0,min(overlap,99));

if (nrecs == 1 & nsamp == 0)
    nsamp = fix(ly/ (8 - 7 * overlap/100));
end

if (nfft < nsamp)   nfft = 2^nextpow2(nsamp); end
if wind ==0; wind = hanning(nsamp); end

overlap  = fix( nsamp * overlap/100);
nadvance = nsamp - overlap;
nrecs    = fix ( (ly*nrecs - overlap) / nadvance);

%% ---------------- accumulate triple products ----------------------

bic  = zeros(nfft,nfft);
Pyy  = zeros(nfft,1);

mask = hankel([1:nfft],[nfft,1:nfft-1] );   % the hankel mask (faster)
Yf12 = zeros(nfft,nfft);
ind  = [1:nsamp];

for k = 1:nrecs
    ys = y(ind);
    ys = (ys(:) - mean(ys)) .* wind;
    Yf = fft(ys,nfft)  / nsamp;
    CYf     = conj(Yf);
    Pyy     = Pyy + Yf .* CYf;
    Yf12(:) = CYf(mask);
    bic = bic + (Yf * Yf.') .* Yf12;
    ind = ind + nadvance;
end

bic     = bic / nrecs;
Pyy     = Pyy / nrecs;
mask(:) = Pyy(mask);
bic = abs(bic).^2 ./ (Pyy * Pyy.' .* mask);
bic = fftshift(bic);

if (rem(nfft,2) == 0)
    waxis = [-nfft/2:(nfft/2-1)]'/nfft;
else
    waxis = [-(nfft-1)/2:(nfft-1)/2]'/nfft;
end

%% ------------ contour plot of magnitude bispectum --------------------
% 
% hold off, clf
% %   contour(bic,4,waxis,waxis),grid
% contour(waxis,waxis,bic,4), grid on
% title('Bicoherence estimated via the direct (FFT) method')
% xlabel('f1'), ylabel('f2')
% set(gcf,'Name','Hosa BICOHER')
% 
% [colmax,row] = max(bic)  ;
% [maxval,col] = max(colmax);
% row = row(col);
% disp(['Max: bic(',num2str(waxis(row)),',',num2str(waxis(col)),') = ', ...
%     num2str(maxval) ])

return
