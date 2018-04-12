%**************************************************************************
%   Hybrid denoising based on the HOS analysis, block thresholding and Wienar 
%   Filtering. 
%**************************************************************************
%
%   [USAGE] 
%   [dn org opt] = wlBlock(data,opt);
% 
%   [INPUTS] 
%   data:    a structure including the information of input data
%   opt:     a structure including parameters needed for the CWT and
%            denoising
% 
%   [OUTPUTS]
%   dn:      a structure including denoised signal and thresholded 
%            coefficients and block threshold_blocksize map. 
%   org:     a structure including wavelet coefficients and scale info 
%            before thresholding. 
%   opt:     parameter structure updated by a rough estimation of the 
%            arrival time. 
% 
%   [References] 
%   [1]  Mousavi S. M., and C. A. Langston (2016). Hybrid Seismic Denoising 
%   Using Wavelet Block Thresholding and Higher Order Statistics, 
%   Bulletin of Seismological Society of America, 106 (4), 1380-1393, 
%   DOI:10.1785/0120150345. 
%
%   [2]  Mousavi S. M., C. A. Langston, and S. P. Horton (2016). Automatic
%   microseismic denoising and onset detection using the synchrosqueezed 
%   continuous wavelet transform, Geophysics, 81 (4), V341-V355, 
%   DOI:10.1785/0120150345. 

%   [3]  Cai, T. (1999). Adaptive wavelet estimation: A block thresholding and 
%   oracle inequality approach, Ann. Stat. 27, 898?924.
%
%   [4]  Cai, T., and H. Zhou (2009). A data-driven block thresholding approach
%   to wave- let estimation, Ann. Stat. 37, no. 2, 569?595, doi: 10.1214/07-AOS538.
%
%-------------------------------------------------------------------------- 
%   S. Mostafa Mousavi 
%   smousavi@memphis.edu
%   Last time modified: June, 2, 2017
%---------------------------------------------------------------------------

%%%   Copyright (c) 2015 S. Mostafa Mousavi
%%%   All rights reserved.
%%%   This software is provided for non-commercial research purposes only. 
%%%   No warranty is implied by this distribution. Permission is hereby granted
%%%   without written agreement and without license or royalty fees, to use, 
%%%   copy, modify, and distribute this code and its documentation, provided 
%%%   that the copyright notice in its entirety appear in all copies of this 
%%%   code, and the original source of this code, 


function [dn org opt] = wlBlock(data,opt);

%%  wavelet Transforming
h = waitbar(0.1,'Wavelet Transforming...');
[wlCoef,as] = cwt_fw(data.noise,opt.type,opt.nv,data.dt);
waitbar(0.8,h,'Wavelet Transforming...'); close(h)
org.c = wlCoef; org.a = as;
[na n] = size(org.c);

if opt.at == 0
%%  a rough arriaval time estimation
% this is based on a simplified form of the characteristic function introduced
% in [2].
ee  = zeros(na,n);
for i = 1:na
v= real(wlCoef(i,:));
a = (v).^2 ;
b = (hilbert(v)).^2; 
ee(i,:) = sqrt(a+b);
end
row = sum(abs(ee));
xr = smooth(row,0.01,'loess');
f = xr - 0.6*max(xr);
f(f < 0 )= 0;
n = length(f);
out = [];

for i = 2:n-1
   if (f(i-1) == 0 & f(i+1)> 0 & f(i) == 0);
    out = [out i+1];
   end
end
opt.at=data.t(out(1))
end

%%  HOS pre-processing
kk = zeros(na,1);
wk = org.c;

for i = 1:na
    w = real(org.c(i,:));
    Vk = 24/length(org.c(i,:)); % eq(12) in [1]
    lk = (sqrt(Vk)/sqrt(1-0.9));
    kr = (sum((w-mean(w)).^4)./length(w)) ./ (var(w,1).^2)-3; % calculating Kurtosis eq (11) in [1]
    if abs(kr) <= lk*opt.gc;
    wk(i,:) = 0;   
    end    
end

wk(isnan(wk))=0;

%% Block Thresholding the Coefficients Scale-by-Scale
% Partitioning;
for i = 1:na
    H{i} = wk(i,:);
end

% Thresholding 
h = waitbar(0.30,'Block Denoising...')

for i = 1:na
    waitbar(i./(na))
    % level dependent noise level estimation from pr-signal noise 
    CH = H{i};
    sigma =  mad(abs(CH(:,1:round(opt.at))))./0.6745; % noise level 
    % estimation as discribed in page 1384 of [1]. 
    [LB{i}, thres{i},T_H{i}] = BtThrsh(H{i},sigma); % block thresholding 
end

% Assembeling 
waitbar(0.30,h,'Preparing the output...')
N = size(wk,2);
t_decom = zeros(1,N);
tlMap = [0,0];

for i = 1:na
    waitbar(i./(na))
    t_decom = [t_decom;T_H{i}];        
    tlMap = [tlMap; LB{i}, thres{i}];  
end
close(h)
clear h

[row col] = size(t_decom);
dnCoef = t_decom(2:row,:);
dnCoef(isnan(dnCoef)) = 0;

%% Recunstructing denoised signal 
x = cwt_iw(dnCoef, opt.type,opt.nv, opt);
dn.c = dnCoef; dn.x = x; dn.m = tlMap;

function [L, lambda,tc] = BtThrsh(subC,sig)
%   This function uses the optimal block size and threshold from Stein's 
%   unbiased risk estimate and threshols coefficients using the block 
%   thresholding improved by Wiener filtering.
%
%   [INPUTS]
%   subC;   wavelet coefficients of the noisy observations
%   sig;    noise level estimation
%
%   [Outputs]
%   L:      the optimal block length 
%   thr:    the optimal threshold (Lamma) 
%   tc:     the thresholded coefficients 
%-------------------------------------------------------------------------- 
%% Estimating the optimal block size and threshold from SURE 
    [L,lambda] = optimalTL(subC,sig); 
    L = round(L);
    
%% partitioning coefficients ino blocks and thresholding coefficients
% in each block by James-Stein's shrinkage rule

    nCol = length(subC);  
    tc = zeros(1,nCol);
    m = 0; Td =0;
    
    for i = 1:nCol            
        m = m + 1; S(m)=(abs(subC(i)).^2)-1;
    end
    s = sum(S(:));
    Td = s/nCol;
    Gammad = (log2(nCol^1.5))./sqrt(nCol);

    if abs(Td) > Gammad
       tc = zeros(1,nCol);
    for j = 1:nCol/L
         
       b2 = (j-1)*L;
       bm = subC(1,b2+1:b2+L);
       % eq (24)in [1] and page 904 of [3]
       a = (1 - lambda*L*(sig)^2./sum(abs(bm(:)).^2)); 
       f = a * (a > 0);
       bt= f*bm; % Block Thresholding
       
       % post-denoising
       % Weiner filtering eq (25)in [1] and page 905 of [3]
       dnWiener = bm.* (abs(bt).^2 ./ (abs(bt).^2 + L*(sig)^2));
       tc(1,b2+1:b2+L) = dnWiener;  
       end
     else
            
       for i = 1:nCol 
       % eq (24)in [1] as on page 904 of [3]
       a = (1 - (2*log10(nCol))./abs(subC(i))^2); 
       f = a * (a > 0);
       bt= f*subC(i);% Block Thresholding
       
       % post-denoising
       % Weiner filtering eq (25)in [1] as on page 905 of [3]
       dnWiener = subC(i).* (abs(bt).^2 ./ (abs(bt).^2 + (sig)^2));
       tc(1,i) = dnWiener; 
       end
    end

      
function [LB,lambda] = optimalTL(sc,sig)
%   This function estimates the optimal threshold and block size using the 
%   Stein's unbiased risk estimation 
%
%   [INPUTS]
%   subb:       noisy subband
%
%   [OUTPUTS]
%   LB:         Optimal block size
%   lambda:     Optimal threshold value
%-------------------------------------------------------------------------- 
% Computting different block sizes based on the restirictions on [4] page 10 
Lmax = ceil(size(sc,2)^0.75); 

p = floor(log2(Lmax));
for i = 1:p
    L_set(i) = 2^i;
end

k = 0;
for L = L_set
    k = k + 1;
    [risk(k),thr(k)] = SURE(sc,L,sig);
end

% Obtain the optimal block size and the corresponding threshold
[guess,ibest] = min(risk);
lambda = thr(ibest); LB = L_set(ibest);


function [guess,lambda] = SURE(sc,L,sig)
% Finding an optimal threshold for a given block size.
%   [INPUTS]
%   subb: inputted noisy subband
%   L: inputted block size
%   sig: noise level estimation
%   
%   [OUTPUTS]
%   guess: outputted minimum risk
%   th: outputted threshold value
%---------------------------------------------------------------------------
sc = real(sc);
% Compute the square sum of all coefficients in a block and generate a threshold series
if L == 1
    sb = sc(:).^2;
    if min(sb) > 0 
       Thres = [sb; 0];
    end
else
    
% the square sum of all atoms in each block
    nCol = length(sc);
    p = mod(nCol, L);
    padSubb = padarray(sc,[0 L-p], 'post');
    blc_n = (length(padSubb)./L);
    sb = zeros(blc_n,1); 
    m = 0;
    
% the square sum of the atoms in each block     
    for i = 1:blc_n            
        bm = padSubb(1,((i-1)*L)+1:i*L);
        m = m + 1; 
        sb(m) = sum(bm(:).^2);               
    end
end

%% Generating the threshold vector
LammaF = floor(2*L*log10(nCol)); % [4], page 10
upth = max(sb):LammaF;
upth = upth';
Thres = [0; abs(sb); upth];
Thres = Thres(find(Thres >= 1)); % Yu et.al (2008)

%% Computing the threshold corresponding to the minimum risk
risk = zeros(length(Thres),1); m = 0;
blc_num = length(sb); Thres = Thres';

for th = Thres
    risk_temp = 0;
    for k = 1:blc_num
    % eq(20) in [1] 
    R = (sig^2)*(L + ((((th^2)*(L.^2))-2*th*L*(L-2))/(sb(k)/(sig^2)))*((sb(k)/(sig^2))>th*L) + ((sb(k)/(sig^2))-2*L)*((sb(k)/(sig^2))<=th*L)); 
    risk_temp = risk_temp+ R; 
    end
    m = m + 1; risk(m) = risk_temp;
end

[guess,ibest] = min(risk);
lambda = Thres(ibest);


function [Wx,as] = cwt_fw(x, type, nv, dt, opt)
% Forward continuous wavelet transform, discretized, as described
% in Mallat, S., Wavelet Tour of Signal Processing 3rd ed.Sec. 4.3.3.
%
% [INPUTS]
%     x: input signal vector.
%  type: wavelet type, string
%    nv: number of voices 
%    dt: sampling period 
%   opt: options structure
%
% [OUTPUTS]
%    Wx: [na x n] size matrix (rows = scales, cols = times)
%    as: na length vector containing the associated scales
%
%---------------------------------------------------------------------------------
%    Modified after a wavelet transform by Eugene Brevdo
%---------------------------------------------------------------------------------
    opt = struct();
    opt.rpadded = 0;

    x = x(:); % Turn into column vector
    n = length(x);

    % Padding the signal 
    N = 2^(1+round(log2(length(x)+eps)));
    n1 = floor((N-n)/2); 
    n2 = n1;if (mod(2*n1+n,2)==1), n2 = n1 + 1; end

    xl = padarray(x(:), n1, 'pre');
    xr = padarray(x(:), n2, 'post');
 
    x = [xl(1:n1); x(:); xr(end-n2+1:end)];

    % Choosing more than this means the wavelet window becomes too short
    noct = log2(N)-1;
    assert(noct > 0 && mod(noct,1) == 0);
    assert(nv>0 && mod(nv,1)==0); 
    assert(dt>0);
    assert(~any(isnan(x)));
    
    na = noct*nv;
    as = 2^(1/nv) .^ (1:1:na);
    
    Wx = zeros(na, N);
    
    x = x(:).';
    xh = fft(x);
    
    % for each octave
    for ai = 1:na
        a = as(ai);
        psih = wfilth(type, N, a, opt);
        xcpsi = ifftshift(ifft(psih .* xh));
        Wx(ai, :) = xcpsi;
    end

    % Shorten W to proper size (remove padding)
    if (~opt.rpadded)
        Wx = Wx(:, n1+1:n1+n);
    end

    % Output a for graphing purposes, scale by dt
    as = as * dt;

function x = cwt_iw(Wx, type, nv, opt)
% The inverse wavelet transform
%
% Implements Eq. (4.67) of Mallat, S., Wavelet Tour of Signal Processing 3rd ed.
%
% Inputs:
%  Wx: wavelet transform of a signal, see help cwt_fw
%  type: wavelet used to take the wavelet transform,
%        see help cwt_fw and help wfiltfn
%  opt: options structure used for forward wavelet transform.
%
% Output:
%  x: the signal, as reconstructed from Wx
%
%---------------------------------------------------------------------------------
%    Modified after a wavelet transform written by Eugene Brevdo
%---------------------------------------------------------------------------------
     opt = struct()
    [na, n] = size(Wx);

    % Padding the signal 
    N = 2^(1+round(log2(n+eps)));
    n1 = floor((N-n)/2); 
    n2 = n1;if (mod(2*n1+n,2)==1), n2 = n1 + 1; end
    Wxp = zeros(na, N);

    Wxp(:, n1+1:n1+n) = Wx;
    Wx = Wxp; clear Wxp;

    noct = log2(N)-1;
    as = 2^(1/nv) .^ (1:1:na);
    
    assert(mod(noct,1) == 0);
    assert(nv>0 && mod(nv,1)==0); 

    % the admissibility coefficient Cpsi
    switch type
      case 'shannon',
        Cpsi = log(2);
      otherwise
        psihfn = wfiltfn(type, opt);
        Cpsi = quadgk(@(x) (conj(psihfn(x)).*psihfn(x))./x, 0, Inf);
    end
    
    % Normalize
    Cpsi = Cpsi / (4*pi);
         
    x = zeros(1, N);
    for ai=1:na
        a = as(ai);
        Wxa = Wx(ai, :);
        psih = wfilth(type, N, a, opt);

        % Convolution theorem 
        Wxah = fft(Wxa);
        xah = Wxah .* psih;
        xa = ifftshift(ifft(xah));
        x = x + xa/a;
    end

     % Take real part and normalize by log_e(a)/Cpsi
     x = log(2^(1/nv))/Cpsi * real(x);

     % Keep the unpadded part
     x = x(n1+1: n1+n);


function [psih] = wfilth(type, N, a, opt)
% Outputs the FFT of the wavelet of family 'type' with parameters
% in 'opt', of length N at scale a: (psi(-t/a))^.
%
% [Inputs]
%   type: wavelet type 
%   N: number of samples to calculate
%   a: wavelet scale parameter 
%   opt: wavelet options 
%   opt.dt: delta t 
%
% [Outputs]
%   psih: wavelet sampling in frequency domain 
%---------------------------------------------------------------------------------
    opt = struct(); 

    k = 0:(N-1);
    xi = zeros(1, N);

    xi(1:N/2+1) = 2*pi/N*[0:N/2];
    xi(N/2+2:end) = 2*pi/N*[-N/2+1:-1];

    psihfn = wfiltfn(type, opt);
    psih = psihfn(a*xi);

    % Normalizing 
    psih = psih * sqrt(a) / sqrt(2*pi);

    % Center around zero in the time domain
    psih = psih .* (-1).^k;


function [psihfn] = wfiltfn(type, opt)
% Wavelet transform function of the wavelet filter in question,
% fourier domain.
%
% [Input]
%   type: string (see below)
%   opt: options structure, e.g. struct('s',1/6,'mu',2)
%
% [Output]
%   psihfn: mother wavelet function ( mexican hat, morlet, shannon, or hermitian)

% Example:
%  psihfn = wfiltfn('bump', struct('mu',1,'s',.5));
%  plot(psihfn(-5:.01:5));
%---------------------------------------------------------------------------------
    switch type
        case 'mhat', % mexican hat
        if ~isfield(opt,'s'), s = 1; else s = opt.s; end
	psihfn = @(w) -sqrt(8)*s^(5/2)*pi^(1/4)/sqrt(3)*w.^2.*exp(-s^2*w.^2/2);
      case 'morlet',
        % can be used with synsq for large enough s (e.g. >5)
        if ~isfield(opt,'mu'), mu = 2*pi; else mu = opt.mu; end
        cs = (1+exp(-mu^2)-2*exp(-3/4*mu^2)).^(-1/2);
        ks = exp(-1/2*mu^2);
        psihfn = @(w)cs*pi^(-1/4)*(exp(-1/2*(mu-w).^2)-ks*exp(-1/2*w.^2));
      case 'shannon',
        psihfn = @(w)exp(-i*w/2).*(abs(w)>=pi & abs(w)<=2*pi);
      case 'hhat', % hermitian hat
        psihfn = @(w)2/sqrt(5)*pi^(-1/4)*w.*(1+w).*exp(-1/2*w.^2);
%       case 'mostafa',
%         load ss
%         if ~isfield(opt,'mu'), mu = 5; else mu = opt.mu; end
%         if ~isfield(opt,'s'), s = 1; else s = opt.s; end
%         psihfnorig = @(w)(0.0720*w.^8)+(0.2746*w.^7)+(0.2225*w.^6)+(-0.2781*w.^5)+(-0.3884*w.^4)+(0.0735*w.^3)+(-0.3354*w.^2)+(-0.0043*w)+(0.3675);
%         psihfn = @(w) psihfnorig((w-mu)/s);
     otherwise
        error('Unknown wavelet type: %s', type);
    end 

    
    function y = mad(x,flag)
%MAD Mean/median absolute deviation. 
%   Y = MAD(X) returns the mean absolute deviation of the values in X.  For
%   vector input, Y is MEAN(ABS(X-MEAN(X)).  For a matrix input, Y is a row
%   vector containing the mean absolute deviation of each column of X.  For
%   N-D arrays, MAD operates along the first non-singleton dimension.
%
%   MAD(X,1) computes Y based on medians, i.e. MEDIAN(ABS(X-MEDIAN(X)).
%   MAD(X,0) is the same as MAD(X), and uses means.
%
%   MAD treats NaNs as missing values, and removes them.
%
%   See also VAR, STD, IQR.

%   References:
%      [1] L. Sachs, "Applied Statistics: A Handbook of Techniques",
%      Springer-Verlag, 1984, page 253.

%   Copyright 1993-2004 The MathWorks, Inc. 
%   $Revision: 2.10.2.2 $  $Date: 2004/01/24 09:34:28 $

% The output size for [] is a special case, handle it here.
if isequal(x,[])
    y = NaN;
    return;
end;

if nargin < 2
    flag = 0;
end

% Figure out which dimension nanmean will work along.
sz = size(x);
dim = find(sz ~= 1, 1);
if isempty(dim)
    dim = 1;
end

% Need to tile the output of nanmean to center X.
tile = ones(1,ndims(x));
tile(dim) = sz(dim);

if flag
    % Compute the median of the absolute deviations from the median.
    y = nanmedian(abs(x - repmat(nanmedian(x), tile)));
else
    % Compute the mean of the absolute deviations from the mean.
    y = nanmean(abs(x - repmat(nanmean(x), tile)));
end


function y = nanmedian(x,dim)
%NANMEDIAN Median value, ignoring NaNs.
%   M = NANMEDIAN(X) returns the sample median of X, treating NaNs as
%   missing values.  For vector input, M is the median value of the non-NaN
%   elements in X.  For matrix input, M is a row vector containing the
%   median value of non-NaN elements in each column.  For N-D arrays,
%   NANMEDIAN operates along the first non-singleton dimension.
%
%   NANMEDIAN(X,DIM) takes the median along the dimension DIM of X.
%
%   See also MEDIAN, NANMEAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.12.2.2 $  $Date: 2004/01/24 09:34:33 $

if nargin == 1
    y = prctile(x, 50);
else
    y = prctile(x, 50,dim);
end


function m = nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along the dimension DIM of X. 
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.13.4.2 $  $Date: 2004/01/24 09:34:32 $

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end



function y = prctile(x,p,dim)
%PRCTILE Percentiles of a sample.
%   Y = PRCTILE(X,P) returns percentiles of the values in X.  P is a scalar
%   or a vector of percent values.  When X is a vector, Y is the same size
%   as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%   the i-th row of Y contains the P(i)-th percentiles of each column of X.
%   For N-D arrays, PRCTILE operates along the first non-singleton
%   dimension.
%
%   Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM.  The
%   DIM'th dimension of Y has length LENGTH(P).
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = prctile(x,50); % the median of x
%      y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
%
%   See also IQR, MEDIAN, NANMEDIAN, QUANTILE.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.12.4.4 $  $Date: 2004/01/24 09:34:55 $

if ~isvector(p) || numel(p) == 0
    error('stats:prctile:BadPercents', ...
          'P must be a scalar or a non-empty vector.');
elseif any(p < 0 | p > 100)
    error('stats:prctile:BadPercents', ...
          'P must take values between 0 and 100');
end

% Figure out which dimension prctile will work along.
sz = size(x);
if nargin < 3 
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1; 
    end
    dimArgGiven = false;
else
    % Permute the array so that the requested dimension is the first dim.
    nDimsX = ndims(x);
    perm = [dim:max(nDimsX,dim) 1:dim-1];
    x = permute(x,perm);
    % Pad with ones if dim > ndims.
    if dim > nDimsX
        sz = [sz ones(1,dim-nDimsX)];
    end
    sz = sz(perm);
    dim = 1;
    dimArgGiven = true;
end

% If X is empty, return all NaNs.
if isempty(x)
    if isequal(x,[]) && ~dimArgGiven
        y = nan(size(p),class(x));
    else
        szout = sz; szout(dim) = numel(p);
        y = nan(szout,class(x));
    end

else
    % Drop X's leading singleton dims, and combine its trailing dims.  This
    % leaves a matrix, and we can work along columns.
    nrows = sz(dim);
    ncols = prod(sz) ./ nrows;
    x = reshape(x, nrows, ncols);

    x = sort(x,1);
    nonnans = ~isnan(x);

    % If there are no NaNs, do all cols at once.
    if all(nonnans(:))
        n = sz(dim);
        if isequal(p,50) % make the median fast
            if rem(n,2) % n is odd
                y = x((n+1)/2,:);
            else        % n is even
                y = (x(n/2,:) + x(n/2+1,:))/2;
            end
        else
            q = [0 100*(0.5:(n-0.5))./n 100]';
            xx = [x(1,:); x(1:n,:); x(n,:)];
            y = zeros(numel(p), ncols, class(x));
            y(:,:) = interp1q(q,xx,p(:));
        end

    % If there are NaNs, work on each column separately.
    else
        % Get percentiles of the non-NaN values in each column.
        y = nan(numel(p), ncols, class(x));
        for j = 1:ncols
            nj = find(nonnans(:,j),1,'last');
            if nj > 0
                if isequal(p,50) % make the median fast
                    if rem(nj,2) % nj is odd
                        y(:,j) = x((nj+1)/2,j);
                    else         % nj is even
                        y(:,j) = (x(nj/2,j) + x(nj/2+1,j))/2;
                    end
                else
                    q = [0 100*(0.5:(nj-0.5))./nj 100]';
                    xx = [x(1,j); x(1:nj,j); x(nj,j)];
                    y(:,j) = interp1q(q,xx,p(:));
                end
            end
        end
    end

    % Reshape Y to conform to X's original shape and size.
    szout = sz; szout(dim) = numel(p);
    y = reshape(y,szout);
end
% undo the DIM permutation
if dimArgGiven
     y = ipermute(y,perm);  
end

% If X is a vector, the shape of Y should follow that of P, unless an
% explicit DIM arg was given.
if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p)); 
end
    

function [c,ww] = smooth(varargin)
%SMOOTH  Smooth data.
%   Z = SMOOTH(Y) smooths data Y using a 5-point moving average.
%
%   Z = SMOOTH(Y,SPAN) smooths data Y using SPAN as the number of points used
%   to compute each element of Z.
%
%   Z = SMOOTH(Y,SPAN,METHOD) smooths data Y with specified METHOD. The
%   available methods are:
%
%           'moving'   - Moving average (default)
%           'lowess'   - Lowess (linear fit)
%           'loess'    - Loess (quadratic fit)
%           'sgolay'   - Savitzky-Golay
%           'rlowess'  - Robust Lowess (linear fit)
%           'rloess'   - Robust Loess (quadratic fit)
%
%   Z = SMOOTH(Y,METHOD) uses the default SPAN 5.
%
%   Z = SMOOTH(Y,SPAN,'sgolay',DEGREE) and Z = SMOOTH(Y,'sgolay',DEGREE)
%   additionally specify the degree of the polynomial to be used in the
%   Savitzky-Golay method. The default DEGREE is 2. DEGREE must be smaller
%   than SPAN.
%
%   Z = SMOOTH(X,Y,...) additionally specifies the X coordinates.  If X is
%   not provided, methods that require X coordinates assume X = 1:N, where
%   N is the length of Y.
%
%   Notes:
%   1. When X is given and X is not uniformly distributed, the default method
%   is 'lowess'.  The 'moving' method is not recommended.
%
%   2. For the 'moving' and 'sgolay' methods, SPAN must be odd.
%   If an even SPAN is specified, it is reduced by 1.
%
%   3. If SPAN is greater than the length of Y, it is reduced to the
%   length of Y.
%
%   4. In the case of (robust) lowess and (robust) loess, it is also
%   possible to specify the SPAN as a percentage of the total number
%   of data points. When SPAN is less than or equal to 1, it is
%   treated as a percentage.
%
%   For example:
%
%   Z = SMOOTH(Y) uses the moving average method with span 5 and
%   X=1:length(Y).
%
%   Z = SMOOTH(Y,7) uses the moving average method with span 7 and
%   X=1:length(Y).
%
%   Z = SMOOTH(Y,'sgolay') uses the Savitzky-Golay method with DEGREE=2,
%   SPAN = 5, X = 1:length(Y).
%
%   Z = SMOOTH(X,Y,'lowess') uses the lowess method with SPAN=5.
%
%   Z = SMOOTH(X,Y,SPAN,'rloess') uses the robust loess method.
%
%   Z = SMOOTH(X,Y) where X is unevenly distributed uses the
%   'lowess' method with span 5.
%
%   Z = SMOOTH(X,Y,8,'sgolay') uses the Savitzky-Golay method with
%   span 7 (8 is reduced by 1 to make it odd).
%
%   Z = SMOOTH(X,Y,0.3,'loess') uses the loess method where span is
%   30% of the data, i.e. span = ceil(0.3*length(Y)).
%
%   See also SPLINE.

%   Copyright 2001-2012 The MathWorks, Inc.

if nargin < 1
    error(message('curvefit:smooth:needMoreArgs'));
end

if nargout > 1 % Called from the GUI cftool
    ws = warning('off', 'all'); % turn warning off and record the previous warning state.
    [lw,lwid] = lastwarn;
    lastwarn('');
else
    ws = warning('query','all'); % Leave warning state alone but save it so resets are no-ops.
end

% is x given as the first argument?
if nargin==1 || ( nargin > 1 && (length(varargin{2})==1 || ischar(varargin{2})) )
    % smooth(Y) | smooth(Y,span,...) | smooth(Y,method,...)
    is_x = 0; % x is not given
    y = varargin{1};
    y = y(:);
    x = (1:length(y))';
else % smooth(X,Y,...)
    is_x = 1;
    y = varargin{2};
    x = varargin{1};
    y = y(:);
    x = x(:);
end

% is span given?
span = [];
if nargin == 1+is_x || ischar(varargin{2+is_x})
    % smooth(Y), smooth(X,Y) || smooth(X,Y,method,..), smooth(Y,method)
    is_span = 0;
else
    % smooth(...,SPAN,...)
    is_span = 1;
    span = varargin{2+is_x};
end

% is method given?
method = [];
if nargin >= 2+is_x+is_span
    % smooth(...,Y,method,...) | smooth(...,Y,span,method,...)
    method = varargin{2+is_x+is_span};
end

t = length(y);
if t == 0
    c = y;
    ww = '';
    if nargout > 1
        ww = lastwarn;
        lastwarn(lw,lwid);
        warning(ws);  % turn warning back to the previous state.
    end
    return
elseif length(x) ~= t
    warning(ws); % reset warn state before erroring
    error(message('curvefit:smooth:XYmustBeSameLength'));
end

if isempty(method)
    diffx = diff(x);
    if uniformx(diffx,x,y)
        method = 'moving'; % uniformly distributed X.
    else
        method = 'lowess';
    end
end

% realize span
if span <= 0
    warning(ws); % reset warn state before erroring
    error(message('curvefit:smooth:spanMustBePositive'));
end
if span < 1, span = ceil(span*t); end % percent convention
if isempty(span), span = 5; end % smooth(Y,[],method)

idx = 1:t;

sortx = any(diff(isnan(x))<0);   % if NaNs not all at end
if sortx || any(diff(x)<0) % sort x
    [x,idx] = sort(x);
    y = y(idx);
end

c = NaN(size(y));
ok = ~isnan(x);
switch method
    case 'moving'
        c(ok) = moving(x(ok),y(ok),span);
    case {'lowess','loess','rlowess','rloess'}
        robust = 0;
        iter = 5;
        if method(1)=='r'
            robust = 1;
            method = method(2:end);
        end
        c(ok) = lowess(x(ok),y(ok),span, method,robust,iter);
    case 'sgolay'
        if nargin >= 3+is_x+is_span
            degree = varargin{3+is_x+is_span};
        else
            degree = 2;
        end
        if degree < 0 || degree ~= floor(degree) || degree >= span
            warning(ws); % reset warn state before erroring
            error(message('curvefit:smooth:invalidDegree'));
        end
        c(ok) = sgolay(x(ok),y(ok),span,degree);
    otherwise
        warning(ws); % reset warn state before erroring
        error(message('curvefit:smooth:unrecognizedMethod'));
end

c(idx) = c;

if nargout > 1
    ww = lastwarn;
    lastwarn(lw,lwid);
    warning(ws);  % turn warning back to the previous state.
end
%--------------------------------------------------------------------
function c = moving(x,y, span)
% moving average of the data.

ynan = isnan(y);
span = floor(span);
n = length(y);
span = min(span,n);
width = span-1+mod(span,2); % force it to be odd
xreps = any(diff(x)==0);
if width==1 && ~xreps && ~any(ynan), c = y; return; end
if ~xreps && ~any(ynan)
    % simplest method for most common case
    c = filter(ones(width,1)/width,1,y);
    cbegin = cumsum(y(1:width-2));
    cbegin = cbegin(1:2:end)./(1:2:(width-2))';
    cend = cumsum(y(n:-1:n-width+3));
    cend = cend(end:-2:1)./(width-2:-2:1)';
    c = [cbegin;c(width:end);cend];
elseif ~xreps
    % with no x repeats, can take ratio of two smoothed sequences
    yy = y;
    yy(ynan) = 0;
    nn = double(~ynan);
    ynum = moving(x,yy,span);
    yden = moving(x,nn,span);
    c = ynum ./ yden;
else
    % with some x repeats, loop
    notnan = ~ynan;
    yy = y;
    yy(ynan) = 0;
    c = zeros(n,1);
    for i=1:n
        if i>1 && x(i)==x(i-1)
            c(i) = c(i-1);
            continue;
        end
        R = i;                                 % find rightmost value with same x
        while(R<n && x(R+1)==x(R))
            R = R+1;
        end
        hf = ceil(max(0,(span - (R-i+1))/2));  % need this many more on each side
        hf = min(min(hf,(i-1)), (n-R));
        L = i-hf;                              % find leftmost point needed
        while(L>1 && x(L)==x(L-1))
            L = L-1;
        end
        R = R+hf;                              % find rightmost point needed
        while(R<n && x(R)==x(R+1))
            R = R+1;
        end
        c(i) = sum(yy(L:R)) / sum(notnan(L:R));
    end
end
%--------------------------------------------------------------------
function c = lowess(x,y, span, method, robust, iter)
% LOWESS  Smooth data using Lowess or Loess method.
%
% The difference between LOWESS and LOESS is that LOWESS uses a
% linear model to do the local fitting whereas LOESS uses a
% quadratic model to do the local fitting. Some other software
% may not have LOWESS, instead, they use LOESS with order 1 or 2 to
% represent these two smoothing methods.
%
% Reference: 
% [C79] W.S.Cleveland, "Robust Locally Weighted Regression and Smoothing
%    Scatterplots", _J. of the American Statistical Ass._, Vol 74, No. 368 
%    (Dec.,1979), pp. 829-836.
%    http://www.math.tau.ac.il/~yekutiel/MA%20seminar/Cleveland%201979.pdf


n = length(y);
span = floor(span);
span = min(span,n);
c = y;
if span == 1
    return;
end

useLoess = false;
if isequal(method,'loess')
    useLoess = true;
end

diffx = diff(x);

% For problems where x is uniform, there's a faster way
isuniform = uniformx(diffx,x,y);
if isuniform
    % For uniform data, an even span actually covers an odd number of
    % points.  For example, the four closest points to 5 in the
    % sequence 1:10 are {3,4,5,6}, but 7 is as close as 3.
    % Therefore force an odd span.
    span = 2*floor(span/2) + 1;

    c = unifloess(y,span,useLoess);
    if ~robust || span<=2
        return;
    end
end

% Turn off warnings when called from command line (already off if called from
% cftool).
ws = warning( 'off', 'MATLAB:rankDeficientMatrix' );
cleanup = onCleanup( @() warning( ws ) );


ynan = isnan(y);
anyNans = any(ynan(:));
seps = sqrt(eps);
theDiffs = [1; diffx; 1];

if isuniform
    % We've already computed the non-robust smooth, so in preparation for
    % the robust smooth, compute the following arrays directly
    halfw = floor(span/2);
    
    % Each local interval is from |halfw| below the current index to |halfw|
    % above
    lbound = (1:n)-halfw;
    rbound = (1:n)+halfw;
    % However, there always has to be at least |span| points to the right of the
    % left bound
    lbound = min( n+1-span, lbound );
    % ... and at least |span| points to the left of the right bound
    rbound = max( span, rbound );
    % Furthermore, because these bounds index into vectors of length n, they
    % must contain valid indices
    lbound = max( 1, lbound );
    rbound = min( n, rbound );
    
    % Since the input is uniform we can use natural numbers for the input when
    % we need them.
    x = (1:numel(x))';
else
    if robust
        % pre-allocate space for lower and upper indices for each fit,
        % to avoid re-computing this information in robust iterations
        lbound = zeros(n,1);
        rbound = zeros(n,1);
    end

    % Compute the non-robust smooth for non-uniform x
    for i=1:n
        % if x(i) and x(i-1) are equal we just use the old value.
        if theDiffs(i) == 0
            c(i) = c(i-1);
            if robust
                lbound(i) = lbound(i-1);
                rbound(i) = rbound(i-1);
            end
            continue;
        end
        
        % Find nearest neighbours
        idx = iKNearestNeighbours( span, i, x, ~ynan );
        if robust
            % Need to store neighborhoods for robust loop
            lbound(i) = min(idx);
            rbound(i) = max(idx);
        end
        
        if isempty(idx)
            c(i) = NaN;
            continue
        end

        x1 = x(idx)-x(i); % center around current point to improve conditioning
        d1 = abs(x1);
        y1 = y(idx);

        weight = iTricubeWeights( d1 );
        if all(weight<seps)
            weight(:) = 1;    % if all weights are 0, just skip weighting
        end

        v = [ones(size(x1)) x1];
        if useLoess
            v = [v x1.*x1]; %#ok<AGROW> There is no significant growth here
        end
        
        v = weight(:,ones(1,size(v,2))).*v;
        y1 = weight.*y1;
        if size(v,1)==size(v,2)
            % Square v may give infs in the \ solution, so force least squares
            b = [v;zeros(1,size(v,2))]\[y1;0];
        else
            b = v\y1;
        end
        c(i) = b(1);
    end
end

% now that we have a non-robust fit, we can compute the residual and do
% the robust fit if required
maxabsyXeps = max(abs(y))*eps;
if robust
    for k = 1:iter
        r = y-c;
        
        % Compute robust weights
        rweight = iBisquareWeights( r, maxabsyXeps ); 
        
        % Find new value for each point.
        for i=1:n
            if i>1 && x(i)==x(i-1)
                c(i) = c(i-1);
                continue;
            end
            if isnan(c(i)), 
                continue; 
            end
            
            idx = lbound(i):rbound(i);
            if anyNans
                idx = idx(~ynan(idx));
            end
            % check robust weights for removed points
            if any( rweight(idx) <= 0 )
                idx = iKNearestNeighbours( span, i, x, (rweight > 0) );
            end
            
            x1 = x(idx) - x(i);
            d1 = abs(x1);
            y1 = y(idx);

            weight = iTricubeWeights( d1 );
            if all(weight<seps)
                weight(:) = 1;    % if all weights 0, just skip weighting
            end

            v = [ones(size(x1)) x1];
            if useLoess
                v = [v x1.*x1]; %#ok<AGROW> There is no significant growth here
            end
            
            % Modify the weights based on x values by multiplying them by
            % robust weights.
            weight = weight.*rweight(idx);
            
            v = weight(:,ones(1,size(v,2))).*v;
            y1 = weight.*y1;
            if size(v,1)==size(v,2)
                % Square v may give infs in the \ solution, so force least squares
                b = [v;zeros(1,size(v,2))]\[y1;0];
            else
                b = v\y1;
            end
            c(i) = b(1);
        end
    end
end
%--------------------------------------------------------------------
function c=sgolay(x,y,f,k)
% savitziki-golay smooth
% (x,y) are given data. f is the frame length to be taken, should
% be an odd number. k is the degree of polynomial filter. It should
% be less than f.

% Reference: Orfanidis, S.J., Introduction to Signal Processing,
% Prentice-Hall, Englewood Cliffs, NJ, 1996.

n = length(x);
f = floor(f);
f = min(f,n);
f = f-mod(f-1,2); % will subtract 1 if frame is even.
diffx = diff(x);
notnan = ~isnan(y);
nomissing = all(notnan);
if f <= k && all(diffx>0) && nomissing, c = y; return; end
hf = (f-1)/2; % half frame length

idx = 1:n;
if any(diffx<0) % make sure x is monotonically increasing
    [x,idx]=sort(x);
    y = y(idx);
    notnan = notnan(idx);
    diffx = diff(x);
end
% note that x is sorted so max(abs(x)) must be abs(x(1)) or abs(x(end));
% already calculated diffx for monotonic case, so use it again. Only
% recalculate if we sort x.
if nomissing && uniformx(diffx,x,y)
    v = ones(f,k+1);
    t=(-hf:hf)';
    for i=1:k
        v(:,i+1)=t.^i;
    end
    [q,~]=qr(v,0);
    ymid = filter(q*q(hf+1,:)',1,y);
    ybegin = q(1:hf,:)*q'*y(1:f);
    yend = q((hf+2):end,:)*q'*y(n-f+1:n);
    c = [ybegin;ymid(f:end);yend];
    return;
end

% non-uniformly distributed data
c = y;

% Turn off warnings when called from command line (already off if called from
% cftool).
ws = warning('off', 'all');
[lastwarnmsg,lastwarnid]=lastwarn;

for i = 1:n
    if i>1 && x(i)==x(i-1)
        c(i) = c(i-1);
        continue
    end
    L = i; R = i;                          % find leftmost and rightmost values
    while(R<n && x(R+1)==x(i))
        R = R+1;
    end
    while(L>1 && x(L-1)==x(i))
        L = L-1;
    end
    HF = ceil(max(0,(f - (R-L+1))/2));     % need this many more on each side

    L = min(n-f+1,max(1,L-HF));            % find leftmost point needed
    while(L>1 && x(L)==x(L-1))
        L = L-1;
    end
    R = min(n,max(R+HF,L+f-1));            % find rightmost point needed
    while(R<n && x(R)==x(R+1))
        R = R+1;
    end

    tidx = L:R;
    tidx = tidx(notnan(tidx));
    if isempty(tidx)
        c(i) = NaN;
        continue;
    end
    q = x(tidx) - x(i);   % center to improve conditioning
    vrank = 1 + sum(diff(q)>0);
    ncols = min(k+1, vrank);
    v = ones(length(q),ncols);
    for j = 1:ncols-1
        v(:,j+1) = q.^j;
    end
    if size(v,1)==size(v,2)
        % Square v may give infs in the \ solution, so force least squares
        d = [v;zeros(1,size(v,2))]\[y(tidx);0];
    else
        d = v\y(tidx);
    end
    c(i) = d(1);
end
c(idx) = c;

lastwarn(lastwarnmsg,lastwarnid);
warning(ws);
%--------------------------------------------------------------------
function ys = unifloess(y,span,useLoess)
%UNIFLOESS Apply loess on uniformly spaced X values

y = y(:);

% Omit points at the extremes, which have zero weight
halfw = (span-1)/2;              % halfwidth of entire span
d = abs((1-halfw:halfw-1));      % distances to pts with nonzero weight
dmax = halfw;                    % max distance for tri-cubic weight

% Set up weighted Vandermonde matrix using equally spaced X values
x1 = (2:span-1)-(halfw+1);
weight = (1 - (d/dmax).^3).^1.5; % tri-cubic weight
v = [ones(length(x1),1) x1(:)];
if useLoess
    v = [v x1(:).^2];
end
V = v .* repmat(weight',1,size(v,2));

% Do QR decomposition
[Q,~] = qr(V,0);

% The projection matrix is Q*Q'.  We want to project onto the middle
% point, so we can take just one row of the first factor.
alpha = Q(halfw,:)*Q';

% This alpha defines the linear combination of the weighted y values that
% yields the desired smooth values.  Incorporate the weights into the
% coefficients of the linear combination, then apply filter.
alpha = alpha .* weight;
ys = filter(alpha,1,y);

% We need to slide the values into the center of the array.
ys(halfw+1:end-halfw) = ys(span-1:end-1);

% Now we have taken care of everything except the end effects.  Loop over
% the points where we don't have a complete span.  Now the Vandermonde
% matrix has span-1 points, because only 1 has zero weight.
x1 = 1:span-1;
v = [ones(length(x1),1) x1(:)];
if useLoess
    v = [v x1(:).^2];
end
for j=1:halfw
    % Compute weights based on deviations from the jth point,
    % then compute weights and apply them as above.
    d = abs((1:span-1) - j);
    weight = (1 - (d/(span-j)).^3).^1.5;
    V = v .* repmat(weight(:),1,size(v,2));
    [Q,~] = qr(V,0);
    alpha = Q(j,:)*Q';
    alpha = alpha .* weight;
    ys(j) = alpha * y(1:span-1);

    % These coefficients can be applied to the other end as well
    ys(end+1-j) = alpha * y(end:-1:end-span+2);
end
%--------------------------------------------------------------------
function isuniform = uniformx(diffx,x,y)
%ISUNIFORM True if x is of the form a:b:c

if any(isnan(y)) || any(isnan(x))
    isuniform = false;
else
    isuniform = all(abs(diff(diffx)) <= eps*max(abs([x(1),x(end)])));
end
%--------------------------------------------------------------------
function idx = iKNearestNeighbours( k, i, x, in )
% Find the k points from x(in) closest to x(i)

if nnz( in ) <= k
    % If we have k points or fewer, then return them all
    idx = find( in );
else
    % Find the distance to the k closest point
    d = abs( x - x(i) );
    ds = sort( d(in) );
    dk = ds(k);
    
    % Find all points that are as close as or closer than the k closest point
    close = (d <= dk);
    
    % The required indices are those points that are both close and "in"
    idx = find( close & in );
end
%--------------------------------------------------------------------
% Bi-square (robust) weight function
function delta = iBisquareWeights( r, myeps )
% Convert residuals to weights using the bi-square weight function.
% NOTE that this function returns the square root of the weights

% Only use non-NaN residuals to compute median
idx = ~isnan( r );
% And bound the median away from zero
s = max( 1e8 * myeps, median( abs( r(idx) ) ) );
% Covert the residuals to weights
delta = iBisquare( r/(6*s) );
% Everything with NaN residual should have zero weight
delta(~idx) = 0;
function b = iBisquare( x )
% This is this bi-square function defined at the top of the left hand
% column of page 831 in [C79]
% NOTE that this function returns the square root of the weights
b = zeros( size( x ) );
idx = abs( x ) < 1;
b(idx) = abs( 1 - x(idx).^2 );
%--------------------------------------------------------------------
% Tri-cubic weight function
function w = iTricubeWeights( d )
% Convert distances into weights using tri-cubic weight function.
% NOTE that this function returns the square-root of the weights.
%
% Protect against divide-by-zero. This can happen if more points than the span
% are coincident.
maxD = max( d );
if maxD > 0
    d = d/max( d );
end
w = (1 - d.^3).^1.5;
