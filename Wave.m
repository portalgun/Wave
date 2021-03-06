classdef Wave < handle
methods(Static)
    function smpPrUnit=poss2smpPerUnit(posSmps)
        % XXX smpFrq2smpPerUnt
        % function smpPerUnit = smpFrq2smpPerUnit(frqCpu)
        %
        % Determine sampling rate from frequency samples
        %
        % Example call: bNonZeroOnly = 0; smpPerUnit = smpFrq2smpPerUnit(smpFrq(13,5,bNonZeroOnly));
        %%%%%%%%%%%
        % frqCpu:       sampled frequencies in cycles per unit
        % smpPerUnit:   sampling rate
        tol = 4;
        if length(unique(round(diff(posSmps),tol))) > 1 && length(unique(round(diff(posSmps),tol)))
            error('pos2smpPerUnit.m : ERROR! Code requires uniformly spaced samples as input');
        end

        %IF POSITIONS START AT ZERO, THEN CENTER
        if posSmps(1) == 0
            ctrInd  = floor(length(posSmps)/2 + 1);
            posSmps = posSmps - posSmps(ctrInd);
        end

        % FOR CENTERED POSITIONS...
                     % neg for even pos for odd
        smpPrUnit = (-1 * ~mod(length(posSmps),2)) * 0.5*max(length(posSmps))./min(posSmps);
    end
    function smpPrUnit=frqs2smpPerUnit(frqCpu)
        % function [smpPerUnit] = smpPos2smpPerUnit(posSmps)
        %
        % Determine sampling rate from positional samples
        %
        % Example call : bStartAtZero = 1; smpPerUnit =  smpPos2smpPerUnit(smpPos(12,27,bStartAtZero));
        %
        % posSmps    : position samples assumed to be generated by smpPos.m
        % smpPrUnit : number of samples per unit
        %


        % THROW ERROR IF SAMPLES ARE NONUNIFORMLY SPACED
        tol = 4;
        if length(unique(round(diff(frqCpu),tol))) > 1 && length(unique(round(diff(frqCpu),tol)))
            error('smpPos2smpPerUnit.m : ERROR! Code requires uniformly spaced samples as input');
        end

        if min(frqCpu) > 0
            negFrqs = length(frqCpu)+ mod(length(frqCpu,2)):-1:1;
            frqCpu  = [ -negFrqs 0 frqCpu];
        end

        % FOR FREQUENCY SAMPLES WITH NEGATIVE SAMPLES ALSO PRESENT
        modu=mod(length(frqCpu),2);
        if modu==0
            smpPrUnit = 2.*max(frqCpu).*((length(frqCpu)) / (length(frqCpu) - 2));
        elseif modu==1
            smpPrUnit = 2.*max(frqCpu);
        end

    end
    function frqCpu = smpFrq(frqSmpsORsmpPerUnit,nSmp,bNonZeroOnly)

        % function frqCpu = Wave.smpFrq(smpPerUnit,numSmp,bNonZeroOnly)
        %
        %   example call: cpd = Wave.smpFrq(128,128,1)
        %
        % sampled frequencies given a sampling rate and a patch
        % NOTE: frequencies should always be sampled at zero
        %
        % smpPerUnit:   sampling rate
        % numSmp:       number of samples
        % bNonZeroOnly: 1 -> returns only the non-zero frequencies
        %               0 -> returns all frequencies (default)
        %
        %               see Wave.smpPos.m
        %
        %%%%%%%%%%%
        % frqCpu:       sampled frequencies in cycles per unit

        if numel(frqSmpsORsmpPerUnit) > 1
            frqSmps=frqSmpsORsmpPerUnit;
            smpPerUnit=Wave.frqs2smpPerUnit(frqSmps);
            if nargin < 2 || isempty(nSmp)
                nSmp=length(frqSmps);
            end
        elseif numel(frqSmpsORsmpPerUnit) ==1
            smpPerUnit=frqSmpsORsmpPerUnit;
        end


        if ~exist('bNonZeroOnly','var') || isempty(bNonZeroOnly)
            bNonZeroOnly=0;
        end

        % IF NUM PIX IS EVEN
        if mod(nSmp,2) == 0
            minPos = -smpPerUnit/2;
            maxPos =  smpPerUnit/2 - smpPerUnit/max(nSmp);
        % IF NUM PIX IS ODD
        elseif mod(nSmp,2) == 1
            minPos = -smpPerUnit/2;
            maxPos =  smpPerUnit/2;
        else
            error(['Wave.smpFrq: WARNING! num pix must be integer valued. numPix = ' num2str(nSmp)]);
        end
        frqCpu   = linspace(minPos,maxPos,max(nSmp));

        if bNonZeroOnly
            frqCpu = frqCpu(frqCpu>0);
        end
    end
    function posUnit = smpPos(posSmpsORsmpPerUnit,nSmp,bStartAtZero)
        %smpPerUnt=[128,68,53];
        %nSmp=[128,72,53];
        %%   example call: posUnt = Wave.smpPos([128,72,53],[128 72 53])
        %%   example call: posUnt = Wave.smpPos([128,72,53],50)
        %%   example call: posUnt = Wave.smpPos(128,128)
        % function posUnt = Wave.smpPos(smpPerUnt,nSmp,bStartAtZero)
        %
        %   example call: posUnt = Wave.smpPos(128,128)
        %
        % spatial positions of samples given a sampling rate and a patch
        % spatial positions are always sampled at zero
        %
        % smpPerUnt:    sampling rate (e.g. smp/deg, smp/sec, etc)
        % nSmp:         number of samples
        % bStartAtZero: boolean determining how samples are positioned
        %               0 -> best for space sampling: e.g. [-3 -2 -1  0  1  2] (default)
        %               1 -> best for time  sampling: e.g. [ 0  1  2  3  4  5]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % posUnt:       spatial positions of sample in units
        %
        %                      *** see Wave.smpFrq.m ***
        %

        if numel(posSmpsORsmpPerUnit) > 1
            posSmps=posSmpsORsmpPerUnit;
            smpPerUnit=Wave.poss2smpPerUnit(posSmps);
            if nargin < 2 || isempty(nSmp)
                nSmp=length(posSmps);
            end
        elseif numel(posSmpsORsmpPerUnit) ==1
            smpPerUnit=posSmpsORsmpPerUnit;
        end
        smpPerUnit=smpPerUnit(:);
        nSmp=nSmp(:);

        assert(Num.isInt(nSmp,'all'),'nSmp must be integer valued');

        if ~exist('bStartAtZero','var') || isempty(bStartAtZero)
            bStartAtZero = 0;
        end

        n=double(mod(nSmp,2)==1);

        posMinUnit = -0.5*(nSmp-n)./smpPerUnit;
        posMaxUnit =  0.5*(nSmp-n)./smpPerUnit - (1./smpPerUnit).*~n;

        if numel(nSmp)==1 & numel(smpPerUnit)==1
            posUnit=linspace(posMinUnit,posMaxUnit,nSmp)-bStartAtZero.*posMinUnit;
            return
        elseif numel(nSmp)==1 & numel(smpPerUnit)>1
            nSmp=repmat(nSmp,numel(smpPerUnit),1);
        elseif numel(nSmp)>1 & numel(smpPerUnit)==1
            smpPerUnit=repmat(nSmp,numel(smpPerUnit),1);
        end
        posUnit    = arrayfun(@(x,y,z) linspace(x,y,z)-bStartAtZero*x ,posMinUnit,posMaxUnit,nSmp,'UniformOutput',false);
    end

    function smpPosUnt = smpFrq2smpPos(smpFrqCpu)
        % function smpPosUnt = smpFrq2smpPos(smpFrqCpu)
        %
        %   example call:
        %
        % convert sample positions to sample frequencies
        %
        % smpFrqCpu: sample frequencies in cyc per units (e.g. cpd, cpp, etc.)
        %%%%%%%%%%
        % smpPosUnt: corresponding samples in position units (e.g. deg, pixels, etc.)

        if ~isvector(smpFrqCpu)
        error(['smpFrq2smpPos: WARNING! smpFrqCpu must be a vector. size(smpFrqCpu)=[' num2str(size(smpFrqCpu,1)) ' ' num2str(size(smpFrqCpu,2)) ']']);
        end

        % UNIQUE SAMPLE FREQUENCIES
        smpFrqUnq = unique(smpFrqCpu(:));
        % SAMPLE SEPARATION
        smpSepUnt = diff(smpFrqCpu(1:2));
        % SAMPLE FREQUENCY
        smpPosUnt = linspace(-.5./smpSepUnt,+.5./smpSepUnt-(1./smpSepUnt)./length(smpFrqCpu),length(smpFrqCpu));
    end

    function smpFrqCpu = smpPos2smpFrq(smpPosUnt)

        % function smpFrqCpu = smpPos2smpFrq(smpPosUnt)
        %
        %   example call: smpPos2smpFrq(smpPos(10,10))
        %
        % convert sample positions to sample frequencies
        %
        % smpPosUnt:  sample positions in units (e.g. deg, pixels, etc.)
        %%%%%%%%%%
        % smpFrqCpu: corresponding samples in frequency units (e.g. cpd, cpp, etc.)



        if numel(smpPosUnt)>1
            % SMP SEPARATION
            smpSepUnt = diff(smpPosUnt(1:2));
            % SMP FREQUENCY
            smpFrqCpu = linspace(-.5./smpSepUnt,+.5./smpSepUnt-(1./smpSepUnt)./length(smpPosUnt),length(smpPosUnt));
        else
            smpFrqCpu = 0;
        end
    end

    function smpPosUntDwn = smpPosDownsample(smpPosUnt,dnK)

        % function Wave.smpPosUntDwn = Wave.smpPosDownsample(Wave.smpPosUnt,dnK)
        %
        %   example call: Wave.smpPosDownsample(Wave.smpPos(104,104),8), Wave.smpPos(13,13)
        %
        % Wave.smpPosUnt:    sample positions in an arbtirary unit
        % dnK:          downsampling factor
        % %%%%%%%%%%%%%%%%%%%%
        % Wave.smpPosUntDwn: downsampled positions


        if dnK < 1
            error(['smpPosDownsample: WARNING! invalid dnK value. dnK=' num2str(dnK)])
        end

        % APPLY DOWNSAMPLING & AVERAGING TO SENSOR LOCATIONS
        modu=mod(length(smpPosUnt)./dnK,2);
        if  modu==0
            indSmp = 1:dnK:numel(smpPosUnt);
        elseif modu == 1
            indSmp = floor(dnK/2+1):dnK:numel(smpPosUnt);
        else
            error(['smpPosDownsample: WARNING! invalid dnK value. dnK=' num2str(dnK)]);
        end

        % DOWNSAPMLE
        smpPosUntDwn = smpPosUnt(indSmp);
    end

end
end
