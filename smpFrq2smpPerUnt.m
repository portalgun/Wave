function smpPerUnt = smpFrq2smpPerUnt(smpCpu)

% function smpFrq2smpPerUnt(smpCpu)
%
%   example call: smpFrq2smpPerUnt(smpFrq(10,10))
%
% converts frequency samples to samples per unit
%
% smpCpu:     frq samples in cycles per unit (e.g. cyc/deg, cyc/meter, etc)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smpPerUnt:  samples per unit                (e.g. deg, meter, etc)
%
%       *** see smpPos2smpPerUnt.m ***

% TOLERANCE
tol     = 1e-12;
% ENSURE ROUNDING DOESN'T MESS UP UNIQUE'S FUNCTION
smpCpu    = roundDec(smpCpu,tol);
% UNIQUE FREQUENCY SAMPLES
cpuUnq = unique(smpCpu(:));
% CONVERT TO POSITION SAMPLES
untUnq = smpFrq2smpPos(cpuUnq);
% DIFFERENCE BETWEEN POSITION SAMPLES IN UNITS
untPerSmp = diff(untUnq(1:2));

% SAMPLES PER UNIT
smpPerUnt = roundDec( 1./untPerSmp, 1e-7);
