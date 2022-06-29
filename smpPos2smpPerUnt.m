function smpPerUnt = smpPos2smpPerUnt(smpUnt)

% function smpPos2smpPerUnt(smpUnt)
%
%   example call: smpPos2smpPerUnt(smpPos(10,10))
%
% converts position samples to smp per unit
%
% smpUnt:      position samples in units of deg, img, meters, etc.
% %%%%%%%%%%
% smpPerUnt:   samples per unit (e.g. pixPerDeg)
%
%       *** see smpFrq2smpPerUnt.m ***   

% TOLERANCE
tol     = 1e-12;
% ENSURE ROUNDING DOESN'T MESS UP UNIQUE'S FUNCTION
smpUnt    = roundDec(smpUnt,tol);
% UNIQUE POSITION SAMPLES
untUnq = unique(smpUnt(:));
% DIFFERENCE BETWEEN POSITION SAMPLES IN UNITS
untPerSmp = diff(untUnq(1:2));
% SAMPLES PER UNIT
smpPerUnt = roundDec( 1./untPerSmp, 1e-7);
