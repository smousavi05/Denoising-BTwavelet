%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The following is a demo to run the Hybrid Seismic Denoising code. This
%  is a sample code and depending on your data might need some
%  modifications. 
%
%  Referece:
%  Mousavi S. M., and C. A. Langston (2016). Hybrid Seismic denoising Using 
%  Higher Order Statistics and Improved Wavelet Block Thresholding, Bulletin 
%  of the Seismological Society of America,106 (4), 1380-1393, 
%  doi:10.1785/0120150345. 
%  
%  Mostafa Mousavi
%  Center for Earthquake Research and Information (CERI), 
%  University of Memphis, Memphis, TN.
%  smousavi@memphis.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all 
tic

% Field seismic data: 'real'
% Synthetic test: 'synth'
test='synth';

% Type od the mother wavelet used for CWT calculation: 'morlet' 'shannon' 'mhat' 'hhat'
opt.type = 'morlet'; 

% Number of voices. You can sellect 8,16,32,64. Higher values result in finer
% CWT but makes the denoising slower.
opt.nv = 16;                 

% Guassian correction factor. This corrects the uncertainties for the 
% estimation of the guassianity and amplifies the pre-processing effects.
% It should be selected highh enough to remove the strong noises outside
% of the signal's frequency band but not too high to remove signal's energy. 
% value of 1.0 means no correction. 
% selecting 0.0 means turning off the pre-processing step.
opt.gc=1; 

% Signal's arrival time estimation
% if you know the onset arrival time inset it in the following variable, 
% otherwise set it to 0.0 and a rough estimation of onset time will be made
% automatically
opt.at = 0;    % Onset arrival time (s), 0 means automatic. 

switch test
    case 'real'
data.nm ='ex_real';  % Microearthquake

disp('Field Seismic Data')   

% Read the noisy data 
[data.t,data.noise,data.hdr] = read_sac(data.nm);
data.dt = data.hdr.times.delta        % delta

%% Denoising vi Block Thresholding
% block thresholding 
[dn org opt] = wlBlock(data,opt);

Xnoisy = data.noise/max(data.noise);
denoised_norm = dn.x/max(dn.x);

%% Plotting outputs
% Original Record Signal
figure;
subplot 221
hold on; 
plot(data.t,Xnoisy)
xlim([min(data.t) max(data.t)]);
grid on
grid minor
title('Original Record ','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',12); 
ylabel('Amplitude (count)','FontSize',12)
ax = gca;
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
hold off
clear title xlabel ylabel

subplot 223
imagesc(data.t, org.a, abs(org.c));
xlim([min(data.t) max(data.t)]);
clim=get(gca,'clim');
clim=[-0.03 0.5]*clim(2);
set(gca,'clim',clim)
title({'CWT Scalogram'},'Rotation',0,'FontSize',14); 
xlabel({'Time  (s)'},'FontSize',12)
ylabel('Scale (1/Hz)','FontSize',12)
ax = gca;
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
hold off
clear title xlabel ylabel

subplot 222
hold on; 
plot(data.t,denoised_norm)
grid on
grid minor
title('Denoised Signal','Rotation',0,'FontSize',14);
xlabel({'Time (s)'},'FontSize',14); 
xlim([min(data.t) max(data.t)]);
ylabel('Amplitude (count)','FontSize',12)
ax = gca;
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
hold off
clear title xlabel ylabel
subplot 224
imagesc(data.t, org.a, abs(dn.c)); 
title('CWT')
clim=get(gca,'clim');
clim=[-0.03 0.5]*clim(2);
set(gca,'clim',clim)
ax = gca;
title({'CWT Scalogram'},'Rotation',0,'FontSize',14); 
xlabel({'Time (s)'},'FontSize',12)
ylabel('Scale (1/Hz)','FontSize',12)
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
hold off
clear title xlabel ylabel

% cheking the SNR  micro
disp('SNR of the original signal')
SNRo = sqrt(mean((Xnoisy(round(opt.at):round(opt.at)*2)).^2))/sqrt(mean((Xnoisy(1:round(opt.at))).^2))
disp('SNR after the denoising')
SNRdn = sqrt(mean((denoised_norm(round(opt.at):round(opt.at)*2)).^2))/sqrt(mean((denoised_norm(1:round(opt.at))).^2))


    case 'synth'      
        
disp('Synthetic Seismic Data') 
data.nm ='ex_synth.mat';   % waveform name

% Read the original 
h = waitbar(0,'Loading...');
load(data.nm);
data.x=synt_80_z;
data.t = linspace(0,(100),length(data.x));
data.dt = 0.1;
close(h)
clear h

% % Random Noise + Real Noise 
load('ex_synth+Noise.mat');
data.noise =syntNoisy3_z;

tic
%%  Denoising vi Block Thresholding
[dn org opt] = wlBlock(data,opt);
toc

data.x=data.x/max(data.x);
Xnoisy = data.noise/max(data.noise);
denoised_norm = dn.x/max(dn.x);

figure;
subplot 311
data.x=data.x/max(data.x);
plot(data.t,data.x)
grid on
grid minor
title('Original Record Signal','Rotation',0,'FontSize',14);
ylabel('Normalized Amplitude','FontSize',12)

ax = gca;
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
hold off
clear title xlabel ylabel

subplot 312
hold on; grid on;
plot(data.t(1:length(data.noise)),data.noise); axis tight
grid on
grid minor
title('Synthetic Signal + Real Seismic Noise','Rotation',0,'FontSize',14);
xlim([min(data.t) max(data.t)]);
ax = gca;
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
hold off
clear title xlabel ylabel

subplot 313
plot(data.t(1:end-1),denoised_norm); axis tight
grid on
grid minor
title('Denoised Signal via Block-Thresholding','Rotation',0,'FontSize',14);
ax = gca;
ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';
hold off
clear title xlabel ylabel

RMSEB = sqrt(mean((data.x(1:end-1) - denoised_norm(:)).^2)) 
snN = sqrt(mean((data.noise(250:450)).^2))/sqrt(mean((data.noise(1:200)).^2))
snrB = sqrt(mean((denoised_norm(250:450)).^2))/sqrt(mean((denoised_norm(1:200)).^2))
[cB,lgs] = xcorr(denoised_norm,data.x(1:end-1),'coeff');
max(cB)

end


function [t,x,hdr] = read_sac(filename)
%   [USAGE] 
%   [t,x,hdr] = read_sac('filename');

%   [INPUTS] 
%   filename:  filename ... the filename/path of the SAC file. 
% 
%   [OUTPUTS]
%   t:   sac function time axis  
%   x:   sac function data
%   hdr: a structure containing header information
%        hdr.data
%        hdr.times
%        hdr.station
%        hdr.event
%        hdr.eventstation
%        hdr.user
%        hdr.info
%        hdr.response
%
%-------------------------------------------------------------------------- 
%   Last time modified: Sep, 25, 2015
%-------------------------------------------------------------------------- 
% check that the file exists
if nargin <1, error('ERROR: filename was not given.'); end

% check that the filename/path is a string
assert(ischar(filename), ...
    'ERROR: read_sac only accepts string filenames.');
%--------------------------------------------------------------------------
% load the sac file
fid=fopen(filename, 'rb');

if (fid==-1)
  disp('can not open input data file format, press CTRL-C to exit \n');
  pause
end

head1=fread(fid, [5, 14], 'float32');
head2=fread(fid, [5, 8], 'int32');
head3=fread(fid, [24, 8], 'char');
head1=head1'; head2=head2'; head3=head3';
npts=head2(2, 5);
x=fread(fid, npts, 'float32');
fclose(fid);
%--------------------------------------------------------------------------
% Get the headers   
% hdr.data
      hdr.data.npts = head2(2,5);
      hdr.data.scale = head1(1,4);
      
% hdr.times
      hdr.times.delta   = head1(1,1); % time increment
      hdr.times.b   = head1(2,1); % begin time
      hdr.times.e   = head1(2,2); % end time
      hdr.times.o   = head1(2,3); % event origin marker
      hdr.times.a   = head1(2,4); % first arrival (P) marker
      hdr.times.t0  = head1(3,1); % time pick 0 (S) marker
      hdr.times.t1  = head1(3,2); % user-defined time pick 1
      hdr.times.t2  = head1(3,3); % user-defined time pick 2
      hdr.times.t3  = head1(3,4); % user-defined time pick 3
      hdr.times.t4  = head1(3,5); % user-defined time pick 4
      hdr.times.t5  = head1(4,1); % user-defined time pick 5
      hdr.times.t6  = head1(4,2); % user-defined time pick 6
      hdr.times.t7  = head1(4,3); % user-defined time pick 7
      hdr.times.t8  = head1(4,4); % user-defined time pick 8
      hdr.times.t9  = head1(4,5); % user-defined time pick 9
      hdr.times.k0  = char(head3(2,9:16)); % event origin time string
      hdr.times.ka  = char(head3(2,17:24)); % first arrival time string
      hdr.times.kt0 = char(head3(3,1:8)); % user-defined pick string 0
      hdr.times.kt1 = char(head3(3,9:16)); % user-defined pick string 1
      hdr.times.kt2 = char(head3(3,17:24)); % user-defined pick string 2
      hdr.times.kt3 = char(head3(4,1:8)); % user-defined pick string 3
      hdr.times.kt4 = char(head3(4,9:16)); % user-defined pick string 4
      hdr.times.kt5 = char(head3(4,17:24)); % user-defined pick string 5
      hdr.times.kt6 = char(head3(5,1:8)); % user-defined pick string 6
      hdr.times.kt7 = char(head3(5,9:16)); % user-defined pick string 7
      hdr.times.kt8 = char(head3(5,17:24)); % user-defined pick string 8
      hdr.times.kt9 = char(head3(6,1:8)); % user-defined pick string 9
      hdr.times.kf = char(head3(6,9:16)); % end of event time string

% hdr.event
      hdr.event.evla = head1(8,1); % event latitude
      hdr.event.evlo = head1(8,2); % event longitude
      hdr.event.evel = head1(8,3); % event elevation 
      hdr.event.evdp = head1(8,4); % event depth
      hdr.event.nzyear = head2(1,1); % GMT time year
      hdr.event.nzjday = head2(1,2); % event time year(Julian)
      hdr.event.nzhour = head2(1,3); % event time hour
      hdr.event.nzmin = head2(1,4); % event time minute
      hdr.event.nzsec = head2(1,5); % event time second
      hdr.event.nzmsec = head2(2,1); % event time millisecond
      hdr.event.kevnm = char(head3(1,9:24)); % event name
      hdr.event.mag = head1(8,5); % event magnitude
      hdr.event.imagtyp = []; % magnitude type
      hdr.event.imagsrc = []; % source of magnitude information

% hdr.station
      hdr.station.stla = head1(7,2); % station latitude
      hdr.station.stlo = head1(7,3); % station longitude
      hdr.station.stel = head1(7,4); % station elevation
      hdr.station.stdp = head1(7,5); % station depth 
      hdr.station.cmpaz = head1(12,3); % component azimuth relative to north
      hdr.station.cmpinc = head1(12,4); % component "incidence angle" reletive to the vertical
      hdr.station.kstnm = char(head3(1,1:8)); % station name
      hdr.stations.kcmpnm = char(head3(7,17:24)); % channel name
      hdr.stations.knetwk = char(head3(8,1:8)); % network name

% hdr.evntstation
      hdr.evsta.dist = head1(11,1); % source receiver distance (km)
      hdr.evsta.az = head1(11,2); % event-station azimuth
      hdr.evsta.baz = head1(11,3); % event-station back azimuth
      hdr.evsta.gcarc = head1(11,4); % great circle distance (deg)
    
% hdr.user
      hdr.user.data = [head1(9,1:5),head1(10,1:5)]; % user-defined variable
      hdr.user.label = [char(head3(6,17:24)),char(head3(7,1:8)),char(head3(7,9:16))];

% hdr.info
      hdr.info.iftype = head2(4,1); % type of file
      hdr.info.idep = head2(4,2); % type of independent variable
      hdr.info.iztype = head2(4,3); % reference time equivalence
      hdr.info.iinst = head2(4,5); % type of recording instrument
      hdr.info.istreg = head2(5,1); % station geographic region
      hdr.info.ievreg = head2(5,2); % event geographic region 
      hdr.info.ievtyp = head2(5,3); % type of event
      hdr.info.iqual = head2(5,4); % quality of data 
      hdr.info.isynth = head2(5,5); % synthetic data flag 

% hdr.response
      hdr.response = [head1(5,2:5),head1(6,1:5),head1(7,1)]; % intrument response parameters

t = [hdr.times.b:hdr.times.delta:(hdr.data.npts-1)*hdr.times.delta+hdr.times.b]';

end