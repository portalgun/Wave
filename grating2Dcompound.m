function [g,G] = grating2Dcompound(x,y,x0,y0,frqCpd,A,ortDeg,phsDeg,dskDmDegX,rmpDmDegX,dskDmDegY,rmpDmDegY,bNormalize,bPLOT)

% function [g,G] = grating2Dcompound(x,y,x0,y0,frqCpd,A,ortDeg,phsDeg,dskDmDegX,rmpDmDegX,dskDmDegY,rmpDmDegY,bNormalize,bPLOT)    
%
%   example call: [g G] = grating2Dcompound(Wave.smpPos(128,128),[],0,0,[3 9],[1 1/3],[0],[0],0.5,0.5,0.5,0.5,1,1);    
%
% implementation based on Geisler_GaborEquations.pdf in /VisionNotes
%
% x:          x position  in visual degrees matrix
%             []    -> 1 deg patch size, 128 smpPerDeg
%             [n]   -> 1 deg patch size,   n smpPerDeg
%             [1xn] -> x(end)-x(1) patch size, length(x) samples
% y:          y position  in visual degrees matrix       
% x0:         x position  of gabor center                      [  scalar  ]
% y0:         y position  of gabor center                      [  scalar  ]
% frqCpd:     frequency   in cycles per deg                    [1 x nComp ]
% A     :     amplitude 
%             [scalar] -> assigns same amplitude to all components
%             [1 x nComp] -> unique amplitude for each component
% ortDeg:     orientation in deg                               
%             [scalar] -> assigns same orientation to all components
%             [1 x nComp] -> unique orientation for each component
% phsDeg:     phase       in deg                                
%             [scalar] -> assigns same phase to all components
%             [1 x nComp] -> unique phase for each component
% dskDmDegX:  window disk size in degrees for x dimension       [ scalar ]
% rmpDmDegX:  window ramp size in degrees for y dimension       [ scalar ]
% dskDmDegY:  window disk size in degrees for x dimension       [ scalar ]
% rmpDmDegY:  window ramp size in degrees for y dimension       [ scalar ]
% bNormalize: 1 -> normalize to vector magnitude of 1
%             0 -> don't
% bPLOT:      1 -> plot
%             0 -> not
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g:          gabor
% G:          gabor in frequency domain 
%             where G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));

if ~exist('bPLOT','var') || isempty(bPLOT) bPLOT = 0; end
if isempty(x) || isempty(y)
   if isempty(x)
   [x, y] = meshgrid(Wave.smpPos(128,128));    
   elseif isscalar(x)
   [x, y] = meshgrid(Wave.smpPos(x,x));
   elseif Vec.is(x)
   [x, y] = meshgrid(x);    
   elseif min(size(x)) > 1
   y = x';    
   end
end
if Vec.is(x) && Vec.is(y) 
   [x, y]=meshgrid(x,y); 
end
% DEFAULT PARAMETER VALUES
if isempty(x0); x0 = 0; end
if isempty(y0); y0 = 0; end
if isempty(frqCpd); error('gabor2DcompoundBW: SPECIFY frqCpd!'); end
if isempty(A); A = ones([1 length(frqCpd)]); end
if isempty(ortDeg); ortDeg = zeros([1 length(frqCpd)]); end
if isempty(phsDeg); phsDeg = zeros([1 length(frqCpd)]); end
if isempty(dskDmDegX) && isempty(rmpDmDegX) && isempty(dskDmDegY) && isempty(rmpDmDegY)
    bWINDOW = 0;
else
    bWINDOW = 1;
end
if isempty(dskDmDegX); dskDmDegX = (max(x(1,:))-min(x(1,:)))./2; end
if isempty(rmpDmDegX); rmpDmDegX = (max(x(1,:))-min(x(1,:)))./2; end
if isempty(dskDmDegY); dskDmDegY = (max(y(:,1))-min(y(:,1)))./2; end
if isempty(rmpDmDegY); rmpDmDegY = (max(y(:,1))-min(y(:,1)))./2; end

% IF ortDeg IS SCALAR
if length(ortDeg)==1
   ortDeg = ortDeg.*ones([1 length(frqCpd)]); 
end

% IF phsDeg IS SCALAR
if length(phsDeg)==1
   phsDeg = phsDeg.*ones([1 length(frqCpd)]); 
end

% IF A IS SCALAR
if length(A)==1
   A = A.*ones([1 length(frqCpd)]); 
end

% MAKE SURE PARAMETER VECTORS ARE SAME LENGTH
if ~(length(frqCpd)==length(ortDeg) && length(ortDeg)==length(phsDeg) && length(phsDeg)==length(A))
   error(['gabor2DcompoundBW: NUMBER OF ELEMENTS MUST BE SAME IN frqCpd, ortDeg,' ...
          ', phsDeg, and A!']); 
end

% CONVERT DISK AND RMP SIZES FROM DEGREES TO PIXELS
dskDmPixX = floor(dskDmDegX*size(x,2)./(max(x(1,:))-min(x(1,:))));
rmpDmPixX = floor(rmpDmDegX*size(x,2)./(max(x(1,:))-min(x(1,:))));
dskDmPixY = floor(dskDmDegY*size(y,1)./(max(y(:,1))-min(y(:,1))));
rmpDmPixY = floor(rmpDmDegY*size(y,1)./(max(y(:,1))-min(y(:,1))));

% CREATE FLATTOP WINDOW
cosWindow1DX = cosWindowFlattop([1 size(x,2)],dskDmPixX,rmpDmPixX,1,0);
cosWindow1DY = cosWindowFlattop([1 size(y,1)],dskDmPixY,rmpDmPixY,1,0);
cosWindow2D = cosWindow1DY'*cosWindow1DX;

for i = 1:length(frqCpd)
    if bWINDOW==1
        % GABOR
        gCmp(:,:,i) = A(i).*cosWindow2D.* ...
                            cos( (2.*pi.*frqCpd(i).*x) + phsDeg(i).*pi./180);
    elseif bWINDOW==0
        % GABOR
        gCmp(:,:,i) = A(i).*cos( (2.*pi.*frqCpd(i).*x) + phsDeg(i).*pi./180);
    else
        error(['grating2Dcompound: invalid bWINDOW parameter ' num2str(bWINDOW) '!']);
    end
end

% ADD COMPONENTS
g = sum(gCmp,3);

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
   imagesc(x(1,:),y(:,1)',g);
   % surf(X,Y,g,'edgecolor','none','facealpha',.8);
   Fig.format('X (deg)','Y (deg)','Space',0,0,18,14);
   axis square;
   grid on
   % view([-39    16]);
   view([0 90]);
%  rotate3d
   caxis(max(abs(minmax(g)))*[-1 1])
   axis xy
   xlim(minmax(x));
   ylim(minmax(y));
   
   % COMPUTE FOURIER TRANSFORM
   G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % PLOT GABOR IN FREQUENCY %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   subplot(1,2,2);
   
   % FIND NATIVE SAMPLING FREQUENCY IN EITHER DIRECTION
    Wave.smpPosPerUnitX = Wave.smpPos2smpPerUnit(x(1,:));
    Wave.smpPosPerUnitY = Wave.smpPos2smpPerUnit(y(:,1));
   % NUMBER OF SAMPLES IN EITHER DIRECTION
    nx = size(x,2);
    ny = size(y,1);
    % FREQUENCIES
    fxcpd = Wave.smpFrq(Wave.smpPosPerUnitX,nx);
    fycpd = Wave.smpFrq(Wave.smpPosPerUnitY,ny);
    %[U,V]=meshgrid(fxcpd,fycpd); 
    
    imagesc(fxcpd,fycpd,abs(G));
    Fig.format('U (cpd)','V (cpd)','Frequency',0,0,18,14);
    axis square;
    axis tight;
    axis xy
    view([0 90]);
    cb = colorbar;
   
    figure(gcf); suptitle(['\Theta=' '[' num2str(ortDeg,3) ']' ...
                           '\Phi=' '[' num2str(phsDeg,3) ']'],22);
    set(cb,'position',[ 0.9159    0.1714    0.0247    0.6405]);
    subplot(s1);
    colormap(cmapBWR(256));
end

end