classdef SISOChan < matlab.System
    % SISOChan Multi-path channel simualtor

    properties
        % Multi-path channel parameters
        dlyPath;   % (npath,1) vector of path delays in seconds
        gainPath;  % (npath,1) vector of path gains in dB
        phasePath;  % (npath,1) vector of path phases in radians

        fsamp;  % Sample rate in Hz

        % Fractional delay class
        fracDly;

    end

    methods
        function obj = SISOChan(opt)
            % Constructor
            arguments
                opt.fsamp (1,1) {mustBeNumeric} = 1.0;
                opt.dlyPath (:,1) {mustBeNumeric,mustBeNonempty} = [];
                opt.gainPath (:,1) {mustBeNumeric,mustBeNonempty} = [];
            end

            % Sets the parameters from the arguments
            obj.dlyPath = opt.dlyPath;
            obj.gainPath = opt.gainPath;
            obj.fsamp = opt.fsamp;
        end
    end

    methods (Access = protected)
        function setupImpl(obj)
              % setup:  This is called before the first step.
              % For the SISO channel class, we will use this point to
              % construct the fractional delay object.  
              
              % Creates a fractional delay object 
              obj.fracDly = dsp.VariableFractionalDelay(...
                'InterpolationMethod', 'Farrow','FilterLength',8,...
                'FarrowSmallDelayAction','Use off-centered kernel',...
                'MaximumDelay', 1024);                           
        end
        
        function resetImpl(obj)
            % reset:  Called on the first step after reset or release.
            
            % Reset the fracDly object
            obj.fracDly.reset();
            
            % Initialize phases to a row vector of 
            % dimension equal to the number of paths with uniform values 
            % from 0 to 2pi
            npath = length(obj.gainPath);
            obj.phasePath = 2*pi*rand(npath,1);
        end
        
        function releaseImpl(obj)
            % release:  Called after the release method
            
            % Release the fracDly object
            obj.fracDly.release();
        end
        
        function y = stepImpl(obj, x)
            % step:  Run a vector of samples through the channel
            
                        
            % TODO:  Compute the delay in samples
            %    dlySamp = ...
            dlySamp=obj.fsamp*obj.dlyPath;
        
            % TODO:  Convert the path gains in obj.gainPath, which are
            % in dB to gain in linear scale 
            %    gainLin = ...
            gainLin=db2mag(obj.gainPath); % 전력으로해야하나 진폭으로해야하나?
                                          % 진폭이 근본이니깐.

            % TODO:  Multiply each path gain by the random phase
            % in obj.phase
            %   gainLin = ...
            % 진폭이주어졌을때 랜덤 위상곱하는법.
            %시간지연은 주파수영역에서 X .* exp(-1j*2*pi*f*tau);하기
            %위상지연은 시간영역에서 x_phase = real(x * exp(1j*theta));하기
            gainLin=gainLin.*exp(1j*obj.phasePath); % 이걸원래신호에곱하면 랜덤 위상지연이..
            
            % TODO:  Use the fracDly object to compute delayed versions of
            % the input x.
            %     xdly = ...
            % The resulting xdly should be nsamp x npath.
            npath=length(obj.gainPath);
            xdly=zeros(length(x),npath);
            for i=1:npath
                xdly(:,i)=obj.fracDly(x,dlySamp(i));
            end
            % TODO:  Multiply each path by the linear gain
            %   xdly = xdly.*...
   
            for i=1:npath
                xdly(:,i)=xdly(:,i).*gainLin(i);
            end
            
            
            

            % TODO:  Sum across the paths
            %   y = sum(...);
            y=sum(xdly,2); %1: 열끼리합 , 2:행끼리합 
            %마지막은 sigma 각경로가 모두 적용된 신호..ㅇㅅㅇ
        end
    end
end