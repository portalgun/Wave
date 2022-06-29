function [g,G] = grating1Dcompound(x,x0,frqCpd,A,phsDeg,dskDmDeg,rmpDmDeg,bNormalize,bPLOT)

% function [g,G] = grating1Dcompound(x,y,x0,y0,frqCpd,A,ortDeg,phsDeg,dskDmDeg,rmpDmDeg,dskDmDegY,rmpDmDegY,bNormalize,bPLOT)    
%
%   example call: [g G] = grating1Dcompound(Wave.smpPos(128,128),0,[3 9],[1 1/3],[0],0.5,0.5,1,1); 
%
% implementation based on Geisler_GaborEquations.pdf in /VisionNotes
%
% x:          x position  in visual degrees matrix
%             []    -> 1 deg patch size, 128 smpPerDeg
%             [n]   -> 1 deg patch size,   n smpPerDeg
%             [1xn] -> x(end)-x(1) patch size, length(x) samples
% x0:         x position  of gabor center                      [  scalar  ]
% frqCpd:     frequency in cycles per deg                      [ 1 x nCmp ]
% A     :     amplitude 
%             [scalar]   -> assigns same amplitude to all components
%             [1 x nCmp] -> unique amplitude for each component
% phsDeg:     phase in deg                                
%             [scalar]   -> assigns same phase to all components
%             [1 x nCmp] -> unique phase for each component
% dskDmDeg:   window disk size in degrees for x dimension       [ scalar ]
% rmpDmDeg:   window ramp size in degrees for y dimension       [ scalar ]
% bNormalize: 1 -> normalize to vector magnitude of 1
%             0 -> don't
% bPLOT:      1 -> plot
%             0 -> not
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g:          gabor
% G:          gabor in frequency domain 
%             where G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));

if ~exist('bPLOT','var') || isempty(bPLOT) bPLOT = 0; end
if isempty(x)
   if isempty(x)
      x = Wave.smpPos(128,128);    
   elseif isscalar(x)
      x = Wave.smpPos(x,x);
   elseif Vec.is(x)
      x = x;    
   end
end

% DEFAULT PARAMETER VALUES
if isempty(x0); x0 = 0; end
if isempty(frqCpd); error('gabor2DcompoundBW: SPECIFY frqCpd!'); end
if isempty(A); A = ones([1 length(frqCpd)]); end
if isempty(phsDeg); phsDeg = zeros([1 length(frqCpd)]); end
if isempty(dskDmDeg); dskDmDeg = (max(x(1,:))-min(x(1,:)))./2; end
if isempty(rmpDmDeg); rmpDmDeg = (max(x(1,:))-min(x(1,:)))./2; end

% IF phsDeg IS SCALAR
if length(phsDeg)==1
   phsDeg = phsDeg.*ones([1 length(frqCpd)]); 
end

% IF A IS SCALAR
if length(A)==1
   A = A.*ones([1 length(frqCpd)]); 
end

% MAKE SURE PARAMETER VECTORS ARE SAME LENGTH
if ~(length(frqCpd)==length(phsDeg) && length(phsDeg)==length(A))
   error(['gabor2DcompoundBW: NUMBER OF ELEMENTS MUST BE SAME IN frqCpd,' ...
          ', phsDeg, and A!']); 
end

% CONVERT DISK AND RMP SIZES FROM DEGREES TO PIXELS
dskDmPix = floor(dskDmDeg*size(x,2));
rmpDmPix = floor(rmpDmDeg*size(x,2));

% CREATE FLATTOP WINDOW
cosWindow1D = cosWindowFlattop([1 length(x)],dskDmPix,rmpDmPix,1,0);

for i = 1:length(frqCpd)
    % GABOR
    gCmp(:,i) = A(i).*cosWindow1D.* ...
                        cos( (2.*pi.*frqCpd(i).*x) + phsDeg(i).*pi./180);
end

% ADD COMPONENTS
g = sum(gCmp,2);

% NORMALIZE TO MAGNITUDE OF 1
if bNormalize
    g = g./sqrt(sum(g(:).^2));
end

if nargout > 1 
    G = fftshift(fft2(ifftshift(g)))./sqrt(numel(g));
end

if bPLOT
   %%%%%%%%%%%%%%%%%%%%%%%
   % PLOT GABOR IN SPACE %
   %%%%%%%%%%%%%%%%%%%%%%%
   figure('position',[411  523  1014  531]);
   s1=subplot(1,2,1); hold on
   plot(x,g);
   % surf(X,Y,g,'edgecolor','none','facealpha',.8);
   Fig.format('X (deg)','Y (deg)','Space',0,0,18,14);
   axis square;
   % view([-39    16]);
   view([0 90]);
%  rotate3d
   axis xy
   xlim(minmax(x));
   
   % COMPUTE FOURIER TRANSFORM
   G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % PLOT GABOR IN FREQUENCY %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   subplot(1,2,2);
   
   % FIND NATIVE SAMPLING FREQUENCY IN EITHER DIRECTION
    Wave.smpPosPerUnitX = Wave.smpPos2smpPerUnit(x(1,:));
   % NUMBER OF SAMPLES IN EITHER DIRECTION
    nx = size(x,2);
    % FREQUENCIES
    fxcpd = Wave.smpFrq(Wave.smpPosPerUnitX,nx);
    %[U,V]=meshgrid(fxcpd,fycpd); 
    
    plot(fxcpd,abs(G));
    Fig.format('U (cpd)','V (cpd)','Frequency',0,0,18,14);
    axis square;
    axis tight;
    axis xy
    view([0 90]);
   
    figure(gcf); suptitle(['\Phi=' '[' num2str(phsDeg,3) ']'],22);
    subplot(s1);
end

end