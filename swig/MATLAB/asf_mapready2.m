clear;
close all;
clc;
%
speedOfLight=299792458.0;
% Generate the list of jar files
jarDirectory = '/media/arwillis/My Passport/SAR/ASF/ASF_MapReady/MATLAB';
dirResult = dir([jarDirectory, filesep, '*.jar']);
jars = cellfun( @(name)fullfile(jarDirectory, name), ...
    {dirResult.name}, 'UniformOutput', false );
% Add jars to the dynamic classpath
javaaddpath( jars );
% Check classpath
%javaclasspath('-dynamic')
import swig.asf.*;
%methodsview('swig.asf.asf_mapready')
%import swig.ASFMainProgram;
%methodsview('swig/ASFMainProgram')
rObj = javaObject('swig.ASFMapReadyJava');
%rObj = javaObject('swig.asf.dataset_sum_rec');
%javaMethod('main',swig.ASFMainProgram,'')
%javaMethod('loadNativeLibraryFromJar',swig.ASFMainProgram);
%inBaseName = '/media/arwillis/My Passport/SAR/data/ASF/E1_19698_STD_L0_F303/E1_19698_STD_L0_F303.000';
%inBaseName = '/media/arwillis/My Passport/SAR/data/ASF/E2_10001_STD_L0_F279/E2_10001_STD_L0_F279.000';
inBaseName = '/media/arwillis/My Passport/SAR/data/ASF/E1_04693_STD_L0_F284/E1_04693_STD_L0_F284.000';
%inBaseName = '/media/arwillis/My Passport/SAR/data/ASF/J1_12380_STD_L0_F361/J1_12380_STD_L0_F361.000';
%inBaseName = '/media/arwillis/My Passport/SAR/data/ASF/R1_60910_ST1_L0_F374/R1_60910_ST1_L0_F374.000';
leaderFilename = strcat(inBaseName,'.ldr');
rawCEOSFilename = strcat(inBaseName,'.raw');
CEOS_Data=javaMethod('parseCEOS','swig.ASFMapReadyJava',leaderFilename);
BinState_Data = javaMethod('parseBinState','swig.ASFMapReadyJava',rawCEOSFilename);
asfMetaData = javaMethod('asfImportCEOS','swig.ASFMapReadyJava',inBaseName);
%
% range parameters
%
asfSARData = asfMetaData.getSar();
rng_samp_rate = CEOS_Data.getRng_samp_rate()*1e+06; % Hz
rng_samp_rate = asfMetaData.getSar().getRange_sampling_rate();
%pulse_dur = CEOS_Data.getRng_length() *1e-06; % 
pulse_dur = BinState_Data.getPulsedur();
pulse_dur = asfMetaData.getSar().getPulse_duration();
%chirp_slope = 15.55*1e6/(CEOS_Data.getRng_length()*1e-6);
chirp_slope = BinState_Data.getSlope();
chirp_slope = 15.55*1e6/asfMetaData.getSar().getPulse_duration();
%rng_samp_rate = 1.896e+07;
%pulse_dur = 3.71e-05;
%chirp_slope = 4.1779e+11;
%
% azimuth parameters
%
PRF = CEOS_Data.getPrf(); % Hz
PRF = asfMetaData.getSar().getPrf();
radar_wavelength = CEOS_Data.getWave_length(); % m
radar_wavelength = asfMetaData.getSar().getWavelength();
SC_vel = BinState_Data.getVel();

%PRF=1679.902394

%radar_wavelength=0.0566666;
%SC_vel=7125.;
%
% compute the range to the radar reflectors
%
% slant range to first pixel
near_range = BinState_Data.getRange_gate()*speedOfLight/2.0;
near_range = asfMetaData.getSar().getSlant_range_first_pixel();
dr = speedOfLight/(2.*rng_samp_rate);
range = near_range+2700*dr;
%
% use the doppler centroid estimated from the data and the
% doppler rate from the spacecraft velocity and range
%
fdcarray = CEOS_Data.getCrt_dopcen();
fdc = javaMethod('double_array1d_getitem','swig.asf.asf_mapready',fdcarray,0);
fdcarray = asfMetaData.getSar().getRange_doppler_coefficients()
fdc = javaMethod('double_array1d_getitem','swig.asf.asf_mapready',fdcarray,3);
fr = 2*SC_vel*SC_vel/(range*radar_wavelength);
%
% Data content specifics
%
numLines = CEOS_Data.getSc_lin()*2;           % number of pulses received 
numLines = asfMetaData.getGeneral().getLine_count();
numSamplesPerLine = BinState_Data.getNSamp(); % samples per pulse (5616 for ERS-1/2)
numSamplesPerLine = asfMetaData.getGeneral().getSample_count();
%
% DC bias for complex values
%
%  The complex signal is digitized at 5 bits per pixel so the numbers range from 0 to 31;
%  This has a built-in bias, i.e., DC offset, of 15.5 for both I/Q
%  channels.
%
I_bias = BinState_Data.getI_BIAS();
Q_bias = BinState_Data.getQ_BIAS();
%
% get some sar data
%
%[cdata,nrow,ncol] = read_rawsar(rawCEOSFilename, I_bias, Q_bias,1,24000);
[cdata,nrow,ncol] = read_rawsar(rawCEOSFilename, I_bias, Q_bias,5000,numLines);
cdata = fliplr(cdata);
%
% generate the range reference function
%
[cref,fcref]=rng_ref(ncol,rng_samp_rate,pulse_dur,chirp_slope);
%
% take the fft of the SAR data
%
fcdata=fft(cdata);
%
% multiply by the range reference function
%
cout=0.*fcdata;
for k=1:nrow
    ctmp=fcdata(:,k);
    ctmp=fcref.*ctmp;
    cout(:,k)=ctmp;
end
clear cdata
%
% now take the inverse fft
%
odata=ifft(cout);
clear cout
%
% generate the azimuth reference function
%
[cazi,fcazi]=azi_ref(nrow,PRF,fdc,fr);
%
% take the column-wise fft of the range-compressed data
%
fcdata=fft(odata');
%
% multiply by the azimuth reference function
%
cout=0.*fcdata;
for k=1:ncol
    ctmp=fcdata(:,k);
    ctmp=fcazi.*ctmp;
    cout(:,k)=ctmp;
end
%
% now take the inverse fft and plot the data
%
odata=ifft(cout);
clear cout fcdata;
%figure(2)
% map=ones(21,3);
% for k=1:21
%     level=0.15*(k+2);
%     level=min(level,1);
%     map(k,:)=map(k,:).*level;
% end
% colormap(map);

% map=ones(21,3);
% for k=1:21;
%     level=0.05*(k+8);
%     level=min(level,1);
%     map(k,:)=map(k,:).*level;
% end
% colormap(map);
% imagesc(abs(odata));

colormap('gray');
amplitude = abs(odata);
%imshow(20*log(abs(odata)),[]);
[ampHist,binCenters] = hist(amplitude(:),1255); 
epdfVals = ampHist/sum(ampHist);
%plot(binCenters,epdfVals);
%hold on;
ecdfVals=cumsum(ampHist)/numel(amplitude);
%plot(binCenters,ecdfVals);
alpha_trim = 0.05;
[maxVal,maxIdx]=min(abs(ecdfVals-(1-alpha_trim)));
[minVal,minIdx]=min(abs(ecdfVals-(alpha_trim)));
amplitude(amplitude>binCenters(maxIdx)) = binCenters(maxIdx);
amplitude(amplitude<binCenters(minIdx)) = binCenters(minIdx);
amplitude=uint8(255*(amplitude-binCenters(minIdx))/binCenters(maxIdx));
x = linspace(0, 0.048 * 550, 551);
y = linspace(0, 0.048 * 550, 551);
imagesc(amplitude);
xlabel('range')
ylabel('azimuth')
title('range and azimuth compressed')
%imagesc(x,y,amplitude);
%imagesc(amplitude, 'XData', [0 25], 'YData', [0 25]);
%imshow(amplitude, [binCenters(minIdx), binCenters(maxIdx)]);
%subplot(2,2,2),imagesc(abs(odata));
% hold on
% plot(x0,y0,'o')
xlabel('range')
ylabel('azimuth')
title('range and azimuth compressed')
%axis([2600,2900,1000,1200])

%****************************************************
function [cref,fcref]=rng_ref(nfft,fs,pulsedur,slope)
%
% routine to compute ERS chirp and its fourier transform
%
% input
% fs - sampling frequency, ts=1./fs
% pulsedur - pulse duration
% slope - chirp slope
%
% set the constants and make npts be odd
%
npts=floor(fs*pulsedur);
ts=1./fs;
if(mod(npts,2.0) == 0.0)
    npts=npts+1;
end
%
% compute the reference function
%
npt2=floor(npts/2.);
t=ts*(-npt2:npt2);
phase=pi*slope*t.*t;
cref1=exp(1i*phase);
%
% pad the reference function to nfft
%
cref=[cref1,zeros(1,nfft-npts)]';
%
% compute the fourier transform
%
fcref=fft(cref)/nfft;
end

%****************************************************
function [cazi,fcazi]=azi_ref(nazi,PRF,fdc,fr)
%
% routine to compute ERS azimuthal chirp and its fourier transform
%
% input
% nazi - number of points in azimuth
% PRF - pulse repitition frequency, ts=1./fs
% fdc - doppler centriod frequency
% fr - doppler frequency rate
%
% set the constants and make npts be odd
%
% For the ERS radar, the synthetic aperture is 1296 points long so we need to read at least
% that many rows into the memory of the computer
npts=min(nazi-1,1296);
ts=1./PRF;
if(mod(npts,2.0) == 0.0)
    npts=npts+1;
end
%
% compute the azimuth chirp function
%
npt2=floor(npts/2.);
t=ts*(-npt2:npt2);
phase=-2.*pi*fdc*t+pi*fr*t.*t;
cazi1=exp(1i*phase);
%
% pad the function to nfft
%
cazi=[cazi1(npt2:npts),zeros(1,nazi-npts), ...
    cazi1(1:npt2-1)]';
%
% compute the fourier transform
%
fcazi=fft(cazi)/nazi;
end

%****************************************************
function [csar,nrow,ncol]=read_rawsar(sar_file, I_bias, Q_bias, offset, nlines)
%
% routine to read and unpack ERS SAR data in DPAF format
%
% each line is 11644 bytes which decomposes as follows:
% 412 bytes of timing information followed by 11232 bytes
% of complex data consisting of 5616 pairs of 1-byte I/Q biased complex 
% valued 8-bit unsigned data
%
fid=fopen(sar_file,'r');
recordSize = 11644;
if (exist('offset','var')==false)
    offset = 0;
end
if (exist('nlines','var')==false)
    fseek(fid, 0, 'eof');    
    filesize = ftell(fid);
    fseek(fid, 0, 'bof');
    nlines = filesize/recordSize;
end
foffset = offset*recordSize;
fseek(fid, foffset, 0);
%
% read the bytes
%
[sar,nsar]=fread(fid, nlines*recordSize, 'uchar');
sar=reshape(sar,recordSize,nsar/recordSize);
nrow=nsar/recordSize;
ncol=5616;
st=fclose(fid);
%
% extract the real and imaginary parts
% and remove the mean value
%
csar=complex(sar(413:2:11643,:)-I_bias,sar(414:2:11644,:)-Q_bias);
end