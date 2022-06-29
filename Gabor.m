classdef Gabor < handle
% References :
% http://www.peterkovesi.com/matlabfns/PhaseCongruency/Docs/convexpl.html
% Batista, Jorge, and Helder Araújo. "Stereoscopic depth perception using a model based on the primary visual cortex." PloS one 8.12 (2013): e80745.
% Shift operator : https://stackoverflow.com/questions/25827916/matlab-shifting-an-image-using-fft
properties
    % TYPE
    bLog=0
    bRadial=0

    % X - size or grid
    %   example call:  [x y]= meshgrid(Wave.smpPos(128,128));
    %                  G = gabor2D(x,y,0,0,3,0,0,.15,.15,0,1);
    szDeg
    szPix
        %OR
    x % grid

    ctr            % x position  of gabor center                        [ 1 x N  ]
                   % x position  in visual degrees matrix
                   % [n]   -> 1 deg patch size,   n sampPerDeg
                   % [1xn] -> x(end)-x(1) patch size, length(x) samples


    % FREQUENCY
    fCpd           % frequency                                          [ 1 x N  ]

    % ORIENTATION
    ort
    ortDeg

    % PHS
    phs
    phsDeg

    %
    sigmaDeg
    sigmaCpd    % ~bLog - standard deviation of spatial gaussian envelope in directions:
                %     (1) bandpass
                %     (2) lowpass
                %  bLog
                %     (1) sigma parameter of log Gabor frequency tuning curve
                %         radial log-gaussian: exp((-(log(f/frqCpd)).^2) ./ (2 * log(sigmaCpd)^2));
                %     (2) standard deviation of Gaussian-shaped orientation tuning curve
    BWOct       % octave bandwidth
    BWOrt       % orientation bandwidth
    BWOrtDeg


    bNormalize=0

% COMPUTED
    nDim
    grid
    %fGrid

    fsCpd

% OUT
    g

end
properties(Hidden)
    bFlip=false
    bFast=false
end
properties(Constant,Hidden)
    SREQ={'BWOct','BWOrt','BWOrtDeg','sigmaDeg','sigmaCpd'};
end
methods(Static)
    function obj=fast(varargin)
        obj=Gabor();
        obj.bFast=true;
        obj.parse(varargin{:});
        obj.get();
    end
end
methods
    function obj=Gabor(varargin)
        if nargin < 1
            return
        end
        obj.parse(varargin{:});
        obj.get();
    end
end
methods(Access=private)
    function get(obj)
        if obj.bLog
            obj.get_log_gabor();
        elseif obj.bRadial
            obj.get_radial_gabor();
        else
            obj.get_gabor();
        end
        if obj.bNormalize
            obj.normalize();
        end
    end
    function parse(obj,varargin)
        P={...
            'fCpd',      [],'Num.is';
            'ctr',       [],'Num.is_e';
            ...
            'x',       [],'Num.is_e';
            'szPix',      [],'Num.isInt';
            'szDeg',      [],'Num.isInt';
            ...
            'BWOct',    [],'Num.is_e';
            'BWOrt',    [],'Num.is_e';
            'BWOrtDeg', [],'Num.is_e';
            'sigmaDeg', [],'Num.is_e';
            'sigmaCpd', [],'Num.is_e';
            ...
            'ort',      [],'Num.is_e';
            'ortDeg',   [],'Num.is_e';
            'phs',      [],'Num.is_e';
            'phsDeg',   [],'Num.is_e';
            ...
            'bLog',      0,'Num.is_e';
            'bRadial',   0,'Num.is_e';
            'bNormalize',0,'Num.is_e';
        };
        if obj.bFast
            if length(varargin)==1 && isstruct(varargin{1})
                flds=fieldnames(varargin{1});
                S=varargin{1};
                for i = 1:length(flds)
                    obj.(flds{i})=S.(flds{i});
                end
            else
                for i = 1:2:nargin-1
                    fld=varargin{i};
                    obj.(fld)=varargin{i+1};
                end
            end
        else
            if length(varargin)==1 && isstruct(varargin)
                Args.parse(obj,P,varargin{1});
            else
                Args.parse(obj,P,varargin{:});
            end
            if all(cellfun(@(x) isempty(obj.(x)), Gabor.SREQ))
                error('No defined envelop width: %s',strjoin(Gabor.SREQ,', '));
            end
        end
        obj.nDim=sum(obj.szPix~=1);
        n=length(obj.szPix);

        flds1D={'bLog','bNormalize','ortDeg'};
        for i = 1:length(flds1D)
            fld=flds1D{i};
            if numel(obj.(fld))>1
                ind=any(obj.(fld));
                if ~any(ind)
                    ind=1;
                end
                obj.(fld)=obj.(fld)(ind);
            end
        end


        % pos
        % XXX 2 double for nyquist?

        %szPix=num2cell(obj.szPix);
        if isempty(obj.szDeg)
            obj.szDeg=ones(1,obj.nDim);
        end
        if numel(obj.fCpd) < obj.nDim
            obj.fCpd=[obj.fCpd zeros(1,obj.nDim-numel(obj.fCpd))];
        end

        % fsCpd
        obj.fsCpd=cell(n,1);
        for i = 1:n
            p=obj.szPix(i);
            d=obj.szDeg(i);
            smps{i}=Wave.smpPos(p,p)*d; % XXX SLOW 2
            obj.fsCpd{i} = Wave.smpFrq(smps{i}); % XXX SLOW 1
        end

        % Grid
        obj.x=cell(obj.nDim,1);
        [obj.x{:}]=ndgrid(smps{:}); % XXX SLOW 3

        % CTR
        if isempty(obj.ctr)
            obj.ctr=zeros(1,obj.nDim);
        end
        if numel(obj.ctr) < obj.nDim
            n=obj.nDim-numel(obj.ctr);
            obj.ctr=[Vec.row(obj.ctr) zeros(1,n)];
        end
        %obj.ctr=num2cell(obj.ctr);

        % BWOrt
        if ~isempty(obj.BWOrtDeg)
            obj.BWOrt=deg2rad(obj.BWOrtDeg);
        end

        % sigmaDeg
        if isempty(obj.sigmaDeg)
            if ~isempty(obj.BWOrt)
                [obj.sigmaDeg,obj.sigmaCpd] = gabor.BWOrt2sigma(obj.fCpd,obj.BWOrt,obj.bLog);
            elseif ~isempty(obj.BWOct)
                [obj.sigmaDeg,obj.sigmaCpd] = gabor.BWOct2sigma(obj.fCpd, obj.BWOct,obj.bLog);
            end
        elseif numel(obj.sigmaDeg)==1 && obj.nDim > 1
            obj.sigmaDeg=repmat(obj.sigmaDeg,1,obj.nDim);

        end

        % BWOrt
        if isempty(obj.BWOrt)
            % XXX thetaDeg
            for i = 1:obj.nDim
                obj.BWOrt(i)=Gabor.sigma2BWOrt(obj.fCpd(i),obj.sigmaDeg(i),obj.ortDeg,obj.bLog);
            end
        end
        % BWOct
        if isempty(obj.BWOct)
            for i = 1:obj.nDim
                obj.BWOct(i)=Gabor.sigma2BWOct(obj.fCpd(i),obj.sigmaDeg(i),[],obj.bLog);
            end
        end
        % sigmaCpd
        if isempty(obj.sigmaCpd)
            obj.sigmaCpd=deg2cpd(obj.sigmaDeg);
        end

        % ort
        if isempty(obj.ort)
            if ~isempty(obj.ortDeg)
                obj.ort = deg2rad(obj.ortDeg);
            else
                obj.ort = 0;
            end
        end
        % phs
        if isempty(obj.phs)
            if ~isempty(obj.phsDeg)
                obj.phs=deg2rad(obj.phsDeg);
            elseif isempty(obj.phsDeg)
                obj.phs=0;
            end
        end

        %% GRID
        %frqs=cell(obj.nDim,1);
        %for i = 1:obj.nDim
        %    frqs{i}=Wave.smpPos(obj.x{i}(i,:));
        %    obj.fsCpd = Wave.smpFrq(frqs{i});
        %end
        %[obj.grid{1:obj.nDim}]=meshgrid(frqs{:});

        %if obj.bLog
        %    [obj.fGrid{1:obj.nDim}]=meshgrid(obj.fsCpd{:});
        %end
    end
    function get_log_gabor(obj)
        % PERFORM IFFT TO GET SPATIAL DOMAIN & APPLY SPATIAL SHIFT OPERATOR IN FREQUENCY DOMAIN
        raw = fftshift(ifft2(ifftshift(obj.get_mag .* obj.get_ang .* exp(-1i*2*pi.*(obj.ctr(1).*obj.fGrid{1} + obj.ctr(2).*obj.fGrid{2})))));

        % APPLY REQUIRED PHASE SHIFT TO SPATIAL DOMAIN FILTER
        obj.g = cos(obj.phs).* real(raw) + sin(obj.phs) .* imag(raw);

        % RETRIEVE FREQUENCY DOMAIN REPRESENTATION TO CHECK
        % LGtest = ifftshift(fft2(fftshift(g)));
    end
    function get_radial_gabor(obj)

        r = sqrt((obj.grid{1}-obj.ctr(1)).^2 + (obj.grid{2}-obj.ctr(2)).^2);

        % GABOR = envelope * carier
        obj.g = exp( -0.5.*( (r-obj.ctr(1))./(obj.sigmaDeg) ).^2 ) .* ...
                cos( 2.*pi.*obj.fCpd.*(r-obj.x0)  + obj.phs );
    end
    function get_gabor(obj)
        % ROTATED POSITIONS: apply rotation matrix
        yp     =  (obj.x{1}-obj.ctr(1)).*cos(obj.ort) + (obj.x{2}-obj.ctr(2)).*sin(obj.ort);
        xp     = -(obj.x{1}-obj.ctr(1)).*sin(obj.ort) + (obj.x{2}-obj.ctr(2)).*cos(obj.ort);

        obj.g = exp(-0.5.*((yp./obj.sigmaDeg(1)).^2 + ...
                           (xp./obj.sigmaDeg(2)).^2)) ...
                .* cos( (2.*pi.*obj.fCpd(1).*yp) + obj.phs(1) ) ...
                .* cos( (2.*pi.*obj.fCpd(2).*xp) + obj.phs(2) );

        % 1D
        %g = exp(          -0.5.*( (x-x0)./  sigmaDeg).^2       ).* ...
        %    cos( 2.*pi.*frqCpd.*(x-x0) + phaseDeg.*pi./180  );
    end
%% HELPERS
    function shft=get_shift_opterator(obj)
        shft=exp(-1i*2*pi.*(obj.x0.*obj.U + obj.y0.*obj.V));
    end
    function ang=get_ang(obj)
        angDeg                        = atan2d(obj.fGrid{1},obj.fGrid{2});               % ANGLE (SIGN CONVENTION: ANTICLOCKWISE POSITIVE)
        ds                            = sind(angDeg) * cosd(obj.ortDeg) - cosd(angDeg) * sind(obj.ortDeg);
        dc                            = cosd(angDeg) * cosd(obj.ortDeg) + sind(angDeg) * sind(obj.ortDeg);
        dtheta                        = abs(atan2d(ds,dc));

        ang                         =  exp((-dtheta.^2) ./ (2 * obj.sigmaDeg^2));
    end
    function mag=get_mag(obj)
        % RADIAL
        PctrRC = fliplr(floor(size(obj.x)/2+1));                    %CENTER COORDINATES IN ROW-COLUMN
        r                             = sqrt(obj.fGrid{1}^2 + obj.fGrid{2}^2);  % RADIUS
        r(PctrRC(1),PctrRC(2))        = 1;                          % ENSURE NONZERO CENTER TO AVOID ISSUES WHILE TAKING LOG

        mag                         = exp((-(log(r/obj.fCpd)).^2) ./ (2 * log(obj.sigmaCpd)^2));
        mag(PctrRC(1),PctrRC(2))    = 0;
    end
    function G=get_fft(obj)
        if obj.bRadial
            [ycpd, xcpd] = ndgrid( Wave.smpFrq(1./diff(obj.x(1,1:2)) ,length(obj.x)) );
            rcpd = sqrt(xcpd.^2 + ycpd.^2);

            G = exp(-0.5.* ((rcpd-obj.fCpd)./obj.sigmaCpd(1)).^2);
        else
            G = fftshift(fft2(fftshift(obj.g)))./sqrt(numel(obj.g));

        end
    end
%% UTIL
    function gnew=copy_fun(obj)
        flds={ ...
            'bLog', ...
            'bRadial', ...
            'szDeg', ...
            'szPix', ...
            'x', ...
            'ctr', ...
            'fCpd', ...
            'ort', ...
            'ortDeg', ...
            'phs', ...
            'phsDeg', ...
            'sigmaDeg', ...
            'sigmaCpd', ...
            'BWOct', ...
            'BWOrt', ...
            'BWOrtDeg', ...
            'bNormalize', ...
            'nDim', ...
            'grid', ...
            'fsCpd', ...
            'g' ...
        };
        gnew=Gabor();
        for i = 1:length(flds)
            gnew.(flds{i})=obj.(flds{i});
        end
    end
    function gnew=combine_fun(obj,g)
        % XXX
        gnew=Gabor();
        if obj.bLog ~= g.bLog
            gnew.bLog=nan;
        else
            gnew.bLog=obj.bLog;
        end
        if ~isequal(obj.ctr,g.ctr)
            gnew.ctr=nan;
        else
            gnew.ctr=obj.ctr;
        end
        % szDeg
        % szPix
        gnew.fsCpd=obj.fsCpd;
        gnew.x=obj.x;
    end
end
methods
    function gnew=plus(obj,g)
        if isnumeric(g)
            gnew=obj.copy_fun();
            gnew.g=gnew.g+g;
        else
            gnew=obj.combine_fun(g);
            gnew.g=obj.g+g.g;
        end
    end
    function R=dot(obj,g)
        R=obj.g(:)'*g.g(:);
    end
    function out=normalize(obj)
        g = obj.g ./ norm(obj.g(:));
        if nargout < 1
            obj.g=g;
        else
            out=g;
        end
    end

    function plot(obj)
        g=gcf;
        p=g.Position(1:2);
        set(g,'position',[p  1014  531]);

        s1=subplot(1,2,1); hold on
        hold off;
        obj.plot_space();

        subplot(1,2,2);
        hold off;
        cb=obj.plot_frq();

        % XXX needed?
        subplot(s1);
        colormap(cmapBWR(256));
        str=obj.get_param_str();
        g=gcf;
        ind=arrayfun(@(x) isa(x,'matlab.graphics.illustration.subplot.Text') && strcmp(x.Type,'subplottext'),g.Children);
        if ~any(ind)
            sgtitle(str);
        else
            c=g.Children(find(ind,1,'first'));
            c.String=str;
        end
       %set(cb,'position',[ 0.9159    0.1714    0.0247    0.6405]);
       %set(cb,'position',[ 0.9159    0.145    0.0247    0.6405]);
    end
    function plot1D(obj)
       % XXX
       % COMPUTE FOURIER TRANSFORM
       Xcpd = Wave.smpPos2Wave.smpFrq(x);
       G = fftshift(fft(ifftshift(g)))./sqrt(numel(g));

       % PLOT GABOR IN SPACE
       figure('position',[411  523  1014  531]);
       subplot(1,2,1);
       plot(x,g,'k','linewidth',2);
       Fig.format('Position (deg)','Weight',[],0,0,20,16);
       axis square;
       ylim(max(abs(ylim)).*[-1 1]);
       Axis.writeText(.02,.3,{['Filter1']},'ratio',15,'left');
       Axis.writeText(.02,.25,{['\mu= '      num2str(x0.*60,2) ' arcmin']},'ratio',15,'left');
       Axis.writeText(.02,.2, {['\sigma  = ' num2str(sigmaDeg.*60,2) ' arcmin']},'ratio',15,'left');
       Axis.writeText(.02,.15,{['f  = '      num2str(frqCpd,3) ' cpd']},'ratio',15,'left');
       Axis.writeText(.02,.1,{['\phi  = '    num2str(phaseDeg,3) ' deg']},'ratio',15,'left');
       Axis.writeText(.02,.05,{['BWoct='        num2str(BWoct,2) ' oct']},'ratio',15,'left');

       % PLOT GABOR IN FREQUENCY
       subplot(1,2,2);
       plot(Xcpd,abs(G),'ko-','linewidth',2);
       Fig.format('Frequency (cpd)','Amplitude',[],0,0,20,16);
       axis square;
       xlim([0 max(abs(Xcpd))]);
       if bNormalize
        ylim([0 .5]);
       else
        ylim([0 Num.round(max(abs(G))+.5,.5)]);
       end

       figure(gcf);
       sgtitle(['Gabor Filter'],20);

    end
    function plot_space(obj)
       % PLOT GABOR IN SPACE %
       Gabor.plotSpace(obj.g,obj.x{1}(:),obj.x{2}(:));
    end
    function CB=plot_frq(obj)
        CB=Gabor.plotFrq(obj.get_fft(),obj.fsCpd{1},obj.fsCpd{2});
    end
    function out=get_param_str(obj)
       bwoct=tostrfun(obj.BWOct,2);
       bwort=tostrfun(obj.BWOrt.*180./pi,2);
       ortdeg=tostrfun(obj.ortDeg,3);
       phsdeg=tostrfun(obj.phsDeg,3);
       frqcpd=tostrfun(obj.fCpd,3);
       sigm=tostrfun(obj.sigmaDeg,3);
       x=tostrfun(obj.ctr,0);
       s='    ';
       out=[ 'BW_{oct}=' bwoct '    '...
             'BW_{\theta}=' bwort '\circ'  newline ...
              ...
             '\Theta=' ortdeg '\circ' s...
             '\Phi=' phsdeg '\circ' s...
             '\omega= ' frqcpd ' cpd' newline...
             'x= ' x s...
             '\sigma= ' sigm '\circ' ...
             ];
       function out=tostrfun(val,sig)
           if sig > 0
               val=round(abs(val),sig);
           else
               val=abs(val);
           end
           out=strrep(Num.toStr(val),',',', ');
           if numel(val) > 1
               out=[ '[' out ']' ];
            end
       end
    end
    function plot_steps(obj)
        figure('position',[486         751        1732         594]);

        titl1=sprintf(' f_{0} = %.2f cpd , \sigma_{f}  = %.2f cpd',obj.fCpd, obj.sigmaCpd);
        titl2=sprintf(' \theta_{0} = %.2f^{\circ}, \sigma_{\theta} = %.2f^{\circ}',obj.ortDeg,obj.sigmaDeg);
        mag=obj.get_ang;
        ang=obj.get_mag;
        shf=mag .* ang .* obj.get_shift_opterator();
        plot_fun(1,abs(mag), titl1);
        plot_fun(2,abs(ang), titl2);
        plot_fun(3,abs(mag.*ang), 'Log Gabor amp. spec.');
        plot_fun(4,real(shf),'Shifted spec.: Real part');
        plot_fun(5,imag(shf),'Shifted spec.: Img part');

        function plot_fun(i,C,titl)
            h=subplot(1,5,i);
            imagesc(obj.U(1,:),obj.V(:,1)',imag(LGshf));
            Fig.format('U (cpd)','V (cpd)',titl,0,0,20,20);
            imagesc(obj.U(1,:),obj.V(:,1)',C);
            axis square; axis tight; axis xy; view([0 90]);
            caxis(caxisLim);colormap(h,jet(256)); colorbar;
            % colormap(hL,jet(256)); colorbar;
        end
    end
end
methods(Static)
    function sigmaCpd=BW2sigma(BWCpd,frqCpd,bLog)
        % function sigmaCpd = bandwidth2sigmaLogGabor(BWfrqCpd,frqCpd)
        %
        % example call : sigmaCpd = bandwidth2sigmaLogGabor(3,2);
        %
        % CALCULATES SIGMA-FREQUENCY PARAMETER FOR A LOG GABOR FROM FREQUENCY BANDWIDTH
        %
        % BWfrqCpd        :  frequency bandwidth, cpd
        % frqCpd          :  preferred("center") frequency of log Gabor frequency tuning curve
        % ********************************************
        % sigmaCpd     : sigma parameter of log Gabor frequency tuning curve i.e. exp((-(log(f/frqCPD)).^2) ./ (2 * log(sigmaCpd/frqCPD)^2));

        if nargin < 3
            bLog=true; % XXX
        end
        if bLog
            % CALCULATE HALF-HEIGHT FREQUENCIES
            fHi            =   (  BWCpd + sqrt(BWfrqCpd.^2 + 4.*(frqCpd.^2)) )./2;       % UPPER HALF WIDTH FREQUENCY IN LINEAR SPACE
            fLo            =   fHi - BWCpd;                                              % LOWER HALF WIDTH REQUENCY IN LINEAR SPACE

            % FWHH OF FREQUENCY TUNING CURVE IN LOG SPACE
            FWHH           =  log(fHi/fLo);

            % STANDARD DEVIATION OF GAUSSIAN IN LOG-SPACE WITH GIVEN FHWW
            sigmaGaussLog  = FWHH./(2.*sqrt(log(4)));

            % CONVERT GAUSSIAN STD INTO SIGMA-FREQUENCY PARAMENTER OF LOG GABOR
            sigmaCpd    = exp(sigmaGaussLog);
        end
    end
    function [sigmaDeg,sigmaCpd]=BWOrt2sigma(f0cpd,BWOrt,bLog)
        % function [sigmaDeg sigmaCpd] = bandwidthOrt2sigma(f0cpd, BWOrt)
        %
        %   example call: sigmaDeg = bandwidthOrt2sigma(4, 42.*pi./180)
        %
        % see derivation in Proof_GaussianBandwidth*.doc in ../VisionNotes/
        %
        % returns standard deviation in low pass direction of 2D gabor
        % given an orientation bandwidth in radians
        %
        % f0cpd:    carrier frequency in cycles/deg
        % BWOrt:    orientation bandwidth in radians
        %           the arc subtended by the width of the gaussian envelope
        %           in the low pass direction of frequency space
        % %%%%%%%%%%%%%%%%%%%%%%
        % sigmaDeg: gaussian standard deviation in the SPACE     domain (deg)
        % sigmaCpd: gaussian standard deviation in the FREQUENCY domain (cpd)
        %
        % ***                   see sigma2bandwidthOct.m                       ***

        if bLog
            % XXX no carrier?
            sigmaDeg = rad2deg(BWort)/(2.*sqrt(log(4)));
        else
            sigmaDeg = sqrt(log(4))./(2.*pi.*f0cpd.*tan(BWOrt./2));
        end
        % XXX correct for log?
        sigmaCpd = 1./(2.*pi.*sigmaDeg);
    end
    function [sigmaDeg, sigmaCpd] = BWOct2sigma(f0cpd, BWOct, bLog)

        % function [sigmaDeg, sigmaCpd] = bandwidthOct2sigma(f0cpd, BWOct)
        %
        %   example call: sigmaDeg = bandwidthOct2sigma(4, 1.5)
        %
        % see derivation in Proof_GaussianBandwidth_*.doc in ../VisionNotes/
        %
        % returns the gaussian standard deviation in the space domain OR
        % in the frequency domain, given the carrier frequency and the
        % octave bandwidth of a gabor function
        %
        % f0cpd:    carrier frequency in cycles/deg
        % BWOct:    octave bandwidth is given by log2( fLo/fHi )
        %           of the gabor's amplitude spectrum
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sigmaDeg: gaussian standard deviation in SPACE     domain (deg)
        % sigmaCpd: gaussian standard deviation in FREQUENCY domain (cpd)
        %
        %             *** see sigma2bandwidthOct.m ***
        if nargin < 3
            bLog=false;
        end
        if bLog
            if BWoct <= .05
                error(['bandwidthOct2sigmaLogGabor: WARNING! BWoct=' num2str(BWoct,'%.2f') ', but it may not be less than 1.0']);
            end

            % OCTAVE BANDWIDTH
            % XXX no carrier?
            sigmaCpd = 2.^(0.5.*BWoct./sqrt(log(4)));
        else
            % STANDARD DEVIATION IN FREQUENCY DOMAIN
            sigmaCpd = f0cpd.*(2.^BWOct-1)./( sqrt(log(4)).*(2.^BWOct + 1) );
        end

        % STANDARD DEVIATION IN SPACE DOMAIN
        % XXX correct for log?
        sigmaDeg=cpd2deg(sigmaCpd);
    end
    function BWfrqCpd=sigma2BW(sigmaCpd,fCpd,bLog)

        if nargin < 3 || isempty(bLog)
            bLog=true; % XXX
        end
        if sigmaCpd <= 1
            error(['sigma2bandwidthOctLogGabor: WARNING! sigmaCpd=' num2str(sigmaCpd,'%.2f') ', but it may not be less than 1.0']);
        end
        if bLog


            mu          = log(fCpd);                                                % MEAN OF *GAUSSIAN* FREQUENCY TUNING CURVE IN *LOG SPACE*
            % sigma       = log(sigmaCpd./fCpd);  % STANDARD DEVIATION OF *GAUSSIAN* FREQUENCY TUNING CURVE IN *LOG SPACE*
            sigma       = log(sigmaCpd);

            [~,logfLoCpd,logfHiCpd]  = Gabor.gaussFWHH(mu,sigma);                                        % FWHH OF FREQUENCY TUNING CURVE IN *LOG SPACE*

            % IF log(sigmaCpd/fCpd) < 0, SWAP logfLoCpd & logfHiCpd
            if log(sigmaCpd) < 0
                logfTmpCpd = logfLoCpd;
                logfLoCpd  = logfHiCpd;
                logfHiCpd  = logfTmpCpd;
            end

            fHi         = exp(logfHiCpd);  % HIGHER HH FREQUENCY IN *LINEAR SPACE*
            fLo         = exp(logfLoCpd);  % LOWER HH FREQUENCY IN *LINEAR SPACE*

            BWfrqCpd    =  fHi - fLo;      % FREQUENCY BANDWIDTH IN *LINEAR SPACE*
        end

    end
    function BWOct=sigma2BWOct(frqCpd,sigmaDeg,sigmaCpd,bLog)
    % function BWoct = sigma2bandwidthOct(frqCpd, sigmaDeg)
    %
    %   example call: BWoct = sigma2bandwidthOct(2, 1)
    %
    % octave bandwidth of a gabor function given the carrier frequency
    % and the standard deviation of the gaussian envelope in space
    %
    % see derivation in Proof_GaussianBandwidth_*.doc in ../VisionNotes/
    %
    % frqCpd:    gabor carrier frequency                          (cyc/deg)
    % sigmaDeg: standard deviation of spatial gaussian envelope  (deg)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BWOct:    octave bandwidth is log2( fHi / fLo )
    %           of the gabor's amplitude spectrum
    %
    % ***                   see bandwidthOct2sigma.m                       ***

        if nargin < 4 || isempty(bLog)
            bLog=false;
        end

        if bLog
            % INPUT CHECKING
            if sigmaCpd <= 1
                error(['sigma2bandwidthOctLogGabor: WARNING! sigmaCpd=' num2str(sigmaCpd,'%.2f') ', but it may not be less than 1.0']);
            end

            % OCTAVE BANDWIDTH
            BWOct = 2.*sqrt(log(4)).*log2( sigmaCpd );


        else
            sigmaCpd=deg2cpd(sigmaDeg);

            % HI AND LO FREQUENCYS AT FWHM
            fHiCpd = frqCpd + sigmaCpd.*sqrt(log(4));
            fLoCpd = frqCpd - sigmaCpd.*sqrt(log(4));

            % fHiCpd = frqCpd + sqrt(log(4))./(2.*pi.*sigmaDeg);
            % fLoCpd = frqCpd - sqrt(log(4))./(2.*pi.*sigmaDeg);

            % OCTAVE BANDWIDTH
            BWOct = log2( fHiCpd ./ fLoCpd );
            % BWoct = log2( ( frqCpd + sqrt(log(4))./(2.*pi.*sigmaDeg)) ./ ( frqCpd - sqrt(log(4))./(2.*pi.*sigmaDeg)) );
        end

    end
    function BWort=sigma2BWOrt(fCpd,sigmaDeg,thetaDeg,bLog)
    % function [sigmaDeg sigmaCpd] = sigma2bandwidthOrt(frqCpd, sigmaDeg)
    %
    %   example call: sigmaDeg = sigma2bandwidthOrt(4, .2)
    %
    % see derivation in Proof_GaussianBandwidth*.doc in ../VisionNotes/
    %
    % returns standard deviation in low pass direction of 2D gabor
    % given an orientation bandwidth in radians
    %
    % frqCpd:   carrier frequency in cycles/deg
    % sigmaDeg: standard deviation of Gaussian envelope
    %           in low pass direction of space domain (deg)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BWort:    orientation (linear) bandwidth is the angle (radians) subtended
    %           by the width of the gaussian envelope in the low pass direction
    %           of frequency space
    %
    % ***                   see sigma2bandwidthOrt.m                       ***

        if nargin < 4 || isempty(bLog)
            bLog=false;
        end
        if ~exist('thetaDeg','var') || isempty(thetaDeg)
            thetaDeg = 0;
        end
        if bLog
            [BWOrtDeg,thetaLdeg,thetaUdeg] = Gabor.gaussFWHH(thetaDeg,sigmaDeg);
            BWort=BWOrtDeg*pi/180;
            % XXX to radians
        else
            BWort = 2.*atan(sqrt(log(4))./(2.*pi.*sigmaDeg.*fCpd));
        end
    end
    function [FWHH xL xU] = gaussFWHH(mu,sigma)
    % function [FWHH xL xU] = gaussFWHH(mu,sigma)
    %
    %   example call: gaussFWHH(0,1)
    %
    % Full Width at Half-Height (FWHH) of gaussian distribution
    %
    % mu:    mean of gaussian
    % sigma: standard deviation of gaussian
    %%%%%%%%%%%%%%%
    % FWHH:  full width at half height

        FWHH = 2.*sigma.*sqrt(log(4));

        if nargout > 1
            xL = mu - sigma.*sqrt(log(4));
        end
        if nargout > 2
            xU = mu + sigma.*sqrt(log(4));
        end
    end
    function [gn,x,y,Gn,fxCpd,fyCpd]=cumu(s,D,prop)
        if ~iscell(prop)
            prop={prop};
        end
        n=size(D,1);

        gn=zeros([s.szPix n]);
        if nargout > 1
            Gn=zeros([s.szPix n]);
        end
        for i = 1:n
            for j = 1:size(D,1)
                s.(prop{1})=D(i,:);
            end
            g=Gabor(s);
            gn(:,:,i)=g.normalize;
            x=g.x{1}(:);
            y=g.x{2}(:);
            if nargout > 2
                Gn(:,:,i)=Gabor.getFFT(gn(:,:,i));
                if i == 1 && nargout > 3
                    fxCpd=g.fsCpd{1}(:);
                    fyCpd=g.fsCpd{2}(:);
                end
            end

        end
    end
    function [gn,x,y,Gn,fxCpd,fyCpd]=cumuAdd(s,D,prop)
        if ~iscell(prop)
            prop={prop};
        end
        G=Gabor(s);
        n=size(D,1);

        gn=zeros([s.szPix n]);
        if nargout > 1
            Gn=zeros([s.szPix n]);
        end
        gg=0;
        for i = 1:n
            for j = 1:size(D,1)
                s.(prop{1})=D(i,:);
            end
            g=Gabor(s);
            gg=g+gg;
            gn(:,:,i)=gg.normalize();
            x=g.x{1}(:);
            y=g.x{2}(:);
            if nargout > 2
                Gn(:,:,i)=Gabor.getFFT(gn(:,:,i));
                if i == 1 && nargout > 3
                    fxCpd=g.fsCpd{1}(:);
                    fyCpd=g.fsCpd{2}(:);
                end
            end
        end

    end
    function CB=plotSpace(g,x,y, bFormat)
        if nargin < 4
            bFormat=true;
        end
        imagesc(x,y,g);
        if bFormat
            Fig.format('X (deg)','Y (deg)','Space',0,0,18,14);
            Fig.format('','');
            axis square;
            axis tight;
            axis xy;
        end
        grid on;
        %view([0 270]);

        m=max(abs(g(:))).*[-1 1];
        caxis(m);
        colormap(cmapBWR(256));

        %zlim(m);
        %cb = colorbar;
        cb=0;
        if nargout > 0
            CB=cb;
        end
    end
    function CB=plotFrq(G,fxCpd,fyCpd, bFormat)
        if nargin < 4
            bFormat=true;
        end
        imagesc(fxCpd,fyCpd,abs(G));
        if bFormat
            Fig.format('U (cpd)','V (cpd)','Frequency',0,0,18,14);
            axis square;
            axis tight;
            grid on;
            axis xy;
        end

        cb=0;
        %cb = colorbar;
        %colormap(cmapBWR(256));
        if nargout > 0
            CB=cb;
        end

        %view([0 90]);
    end
    function G=getFFT(g)
        G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));
    end
    function test()

        Opts=struct();
        Opts.bLog=false;
        Opts.bRadial=false;

        Opts.pos
        Opts.ctr

        Opts.fCpd=1;

        %Opts.ort
        Opts.ortDeg=90;

        %Opts.phs
        Opts.phsDeg=0;

        Opts.BWOct=8.1;
        %Opts.BWOrt
        %Opts.BWOrtDeg
        Opts.sigmaDeg=1;
        %Opts.sigmaCpd

        Opts.bNormalize

        obj.plot();
    end

end
end
