classdef  WLANRx <  matlab.System
    properties
        % Packet properties
        psduLen = 2048;  % PSDU length in Bytes
        mcs = 3;         % MCS selection
        cfg;             % WLAN Toolbox configuration
        fsamp;        % Sample rate in MHz

        % Preamble components
        lstf;   % Legacy STF (short training field)
        lltf;   % Legacy LTF (long training field)

        % Packet length including the preamble in samples
        pktLen = 0; 

        % Detection results
        tSTF = 0.25;    % STF threshold
        rhoSTF;     % STF correlation
        pktFound;   % Flag indicating if packet was found
        istf;       % Index of the sample on which STF was found
        rhoLTF;     % LTF correlation
        iltf;       % Index of the sample on which the LTF was found
        snrDetected;  % SNR of detected packet in dB

    end

    methods
        function obj = WLANRx(opt)
            % Constructor
            arguments
                opt.psduLen (1,1) {mustBeInteger} = 512;
                opt.mcs (1,1) {mustBeInteger} = 3;
            end

            % Set properties
            obj.psduLen = opt.psduLen;
            obj.mcs = opt.mcs;

            % Create the packet configuration
            obj.cfg = wlanNonHTConfig('ChannelBandwidth', 'CBW20', ...
                'PSDULength', obj.psduLen,'MCS',obj.mcs);

            % TODO:  Get the sample rate from the
            % wlanSampRate function using obj.cfg.
            %      obj.fsamp = wlanSampleRate(obj.cfg)            

            % Compute the preamble components
            obj.compPreamble();

        end

        function compPreamble(obj)
            % -------wlanLSTF % wlanLLTF 는 해당 설정에맞는 프리엠블을 만들어줍니다,----------
            %   obj.lstf = wlanLSTF(...);
            %   obj.lltf = wlanLLTF(...);   
            obj.lstf=wlanLSTF(obj.cfg);
            obj.lltf=wlanLLTF(obj.cfg);
        end


        function detectSTF(obj, r)
            % Detects the STF using cyclic correlation.  

            % Since we only want to perform detection on the samples 
            % at the beginning of the packet, we remove samples from the 
            % end.  This will ensure that after the pre-amble is detected,
            % we can decode the remainder of the packet
            sampsToRemove = max(0, obj.pktLen - length(obj.lstf) - length(obj.lltf)); % 이거왜 제거하는지?..
            r = r(1:length(r)-sampsToRemove);

            % Get parameters
            n = length(r);
            period = 16; % 16샘플 반복
            len = length(obj.lstf); % STF 160 샘플
            tol = 1e-8;

            % TODO:  Compute rhoSTF as described in wlanPreamble.mlx
            %   obj.rhoSTF = ...
            % 고정지연 계산 알고리즘:  특정(160)구간내에 일정샘플반복이있는가?
            % 
             
              b = ones(len-period,1); % 144샘플
      
              rprod = conj(r(1:n-period)).*r(1+period:n); % 16샘플지연과 상관값구하기
               rprodAvg = conv(rprod, b, 'valid');  % 상관값의 총누적
            
               %  rhoSTF(i) = abs(rprodAvg(i))^2 / rsqAvg(i) / rsqAvg(i+period);
               % rsQ랑 rsqAVG가 필요한거있지 그래서 rpod보고 참고함
              rsq=abs(r(1:n-period)).^2;
              rsq2=abs(r(1+period:end)).^2;
              rsqAvg=conv(rsq,b,'valid');
              rsq2Avg=conv(rsq2,b,'valid');
              obj.rhoSTF=abs(rprodAvg).^2 ./ (rsqAvg+tol) ./ (rsq2Avg + tol); % 고정지연 상관값 이론식 왜?.. 몰라 (배운적없음)
            % TODO: 최대값찾기?
            %    rhoMax = max(obj.rhoSTF)   
            %    obj.istf = index i where obj.rhoSTF(i) is maximum
                [rhoMax,i] = max(obj.rhoSTF);
                 obj.istf=i;
            % TODO:  obj.pktFound to 1 if rhoMax > obj.tSTF.  Otherwise 
            % 0.25 보다큰가?
            obj.pktFound= rhoMax>obj.tSTF;
             
        end

        function pktDetect(obj, r)
            % Packet detection and initial timing synchronization
            % Run the STF correlation and detection
            obj.detectSTF(r);

            % Exit if no packet is found.
            if ~obj.pktFound
                return;
            end

            % The STF is not exact.  So, we need to search for the true
            % initial location over a range of symbols.  The code below
            % extracts data from a range of r around the detected STF.


            n = length(r);
            i1 = max(obj.istf - length(obj.lstf),1); % 신호 첫시작위치=감지된 stf 위치 - lstf
            i2 = min(obj.istf + 4*length(obj.lltf), n); % i2 는 lstf+lltf 샘플 인덱스/ 4를곱해서 충분히 여유를둠
            r1 = r(i1:i2);

            % TODO:  Compute the un-normalized MF of the partial RX signal
            % r1 with the target, obj.lltf;
            %    z = conv(...)    -
            z=conv(r1,conj(flipud(obj.lltf)),'valid'); % ->상관관계식 ㅇㅇ\
                   
            % TODO:  Compute the correlation squared
           b=ones(length(obj.lltf),1);
               tol = 1e-8;
               Er =  conv(abs(r1).^2,b,'valid'); % 개복잡하노
               Ex =  sum(abs(obj.lltf).^2);
               obj.rhoLTF = abs(z).^2 ./ (Er+tol)./(Ex+tol); % 위에 detect STF랑 좀다른데. 두개만 슬라이딩 하나는 상수


            % TODO:  Find the location of the maximum in obj.rhoLTF.  Add
            % the estimated to location to obj.istf.  Store in obj.iltf.
            % This is the estimated start of the packet.
            %
            %    obj.iltf = ...            

           [~,kmax]=max(obj.rhoLTF);
     
            % Add the initial point

            obj.iltf = i1 + kmax-1; % 맞음 -1해줘야되네 구간을 나눠가지고 원본구간에서의 인덱스는-1 해줘야된다

        end


    end
end