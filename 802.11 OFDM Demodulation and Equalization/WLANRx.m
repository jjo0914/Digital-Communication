classdef  WLANRx <  matlab.System
    properties
        % Packet properties
        cfg;        % WLAN Toolbox configuration
        fsamp;      % Sample rate in MHz
        info;       % OFDM info
        nsampSym;   % number of samples per OFDM symbol                
        ofdmFrameOff;    % OFDM frame offset in samples
        nsymData;   % number of data symbols 
        nsymPre = 3;  % number of pre-amble symbols from start of LTF

        % TX and RX legacy LTF field in frequency domain
        symLTFTx;
        symLTFRx;
      
        % Channel estimate outputs
        chanEstRaw; 
        chanEst;
        noiseEst;
        snrEst;

        % Equalized symbols
        eqMethod;   % equalization method 
        symEq;      % payload
        symDataEq;  % data

    end

    methods
        function obj = WLANRx(opt)
            % Constructor
            arguments
                opt.ofdmFrameOff {mustBeInteger} = 4;
                opt.eqMethod = 'inversion';
            end

            % Set properties
            obj.ofdmFrameOff = opt.ofdmFrameOff;
            obj.eqMethod = opt.eqMethod;

            % Create the packet configuration
            obj.cfg = wlanNonHTConfig('ChannelBandwidth', 'CBW20');

            % Get the OFDM info
            obj.info = wlanNonHTOFDMInfo('NonHT-Data');

            % Get the number of samples per sym
            obj.nsampSym = obj.info.FFTLength + obj.info.CPLength;

            % Get the sample rate from the
            obj.fsamp = wlanSampleRate(obj.cfg);

        end

        function setPktData(obj, opt)  
            % Sets the packet data for demodulation
            arguments
                obj;
                opt.symLTFTx (52,2);
                opt.nsymData {mustBeInteger};
            end

            obj.symLTFTx = opt.symLTFTx;
            obj.nsymData = opt.nsymData;

        end

        function compChanEstRaw(obj, r, indltf)
            % Compute the raw channel estimate

            % Set the initial timing of the LTF.  We adjust this a little
            % before the strongest detected path to ensure that all the
            % paths are in the CP
            i0 = max(indltf - obj.ofdmFrameOff,1); % 그냥인덱스를 frmaeoff만큼조정?
                                                   % cp를포함하기위헤 4샘플?
                                                   % 상관위치가 좀뒤에서 나와서 좀앞으로땅긴대

            % TODO:  Get the section of the RX signal r
            % corresponding to the LTF.  This should start at i0
            % and will be same length as the LTF field
            %   rlltf = r(...);
            rlltf=r(i0:i0+160-1);

            % TODO:  Demodulate the lltf using the ofdmdemod() function.
            % You can get the FFT and CP length from obj.info.  Store
            % the values in obj.symLTFRx
            %    obj.symLTFRx = ofdmdemod(...)
            obj.symLTFRx=ofdmdemod(rlltf,obj.info.FFTLength,obj.info.CPLength);
            % 160샘플이니깐 64*2 

            
            % TODO:  Extract the active indices which are actively used
            % based on obj.info.ActiveFFTIndices.  Store these again in 
            %       obj.symLTFRx = obj.symLTFRx(...,:)
                obj.symLTFRx=obj.symLTFRx(obj.info.ActiveFFTIndices,:);

            % TODO:  Compute the raw channel estimate by dividing the RX
            % signal with the TX signal on the LTF.  The resulting matrix,
            % obj.chanEstRaw, should be 52 x 2.
            %    obj.chanEstRaw

            % 채널= 단순나누기, 채널이 행렬? 52*2 ? 2가지경우로 추정하는건지?-> 맞음
            obj.chanEstRaw=obj.symLTFRx./obj.symLTFTx;


        end

        function compChanEst(obj, r, indltf)

            % Compute the raw channel estimate
            obj.compChanEstRaw(r, indltf);
            
            % TODO:  Average chanEstRaw over the two time symbols
            %   chanEstAvg = mean(...)
            chanEstAvg=mean(obj.chanEstRaw,2); % dim:2 -> dim:1 |

            % Smooth the channel estimate over frequency using the
            % smoothdata function with 'sgolay' method with a deg

            % Sgolay: ofdm 간단한 채널추정용? (SNR낮을때 하면좋대) 
            % 원리가뭐임..?ㅜ
            deg = 4;
            obj.chanEst = smoothdata(chanEstAvg, 'sgolay','Degree',deg);
            % -> 진짜 채널이라 가정하네?

            % Compute the residual error 
            d = obj.symLTFRx - obj.chanEst.*obj.symLTFTx;
            % chanEstRAW=R/T 인데 ,chanEstRaw를 스무딩햇으니 잔차가생김 ->근데?
            % 즉 obj.chanEstRaw=추정H , obj.chanEst=ㄹㅇH
            %  R=H^X (추정)   /  R=HX+W (실제) H->스무딩한건데 실제로취급
            % H^X -HX= W (잔차) 노이즈취급 ㅇㅇ 

            % TODO:  Compute the noise estimate by taking the average
            % of abs(d).^2
            %    obj.noiseEst = ...
            % d가 행렬이2개네 ㅇㅅㅇ
            obj.noiseEst = mean( abs(d).^2,'all') % 시간영역에서의 평균?전력?.

            % Compute an estimate of the SNR
            % 그냥 전력끼리 나누네
            obj.snrEst = mean(abs(obj.symLTFRx).^2, 'all') / obj.noiseEst;
            obj.snrEst = pow2db(obj.snrEst);

            
        end

        function symRx = eqSym(obj,r,indltf)

            % Get the starting sample for the for the payload portion
            i0 = indltf - obj.ofdmFrameOff + obj.nsymPre*obj.nsampSym; 
            % i0는 ltf위치에서 조금 앞으로땡기고 +  LTF(2개) + L-SIG(1개) 위치로조정
            % nsampsym 80(1심볼당 8샘플 ㅇㅇ) * nsymPre 3 (LTF(2개) + L-SIG(1개))
            % 즉 io는 payload시작위치

            % Get the number of symbols to decode      
            n = length(r); % r크기 5273..    
            % 전체크기를 80으로 나눠서 r에있는 ofdm심볼갯수를 구함 (원래 48개(data43개) ofdm심볼 보냇음)
            nsym = min( floor((n-i0)/obj.nsampSym), obj.nsymData );

            % Find the final sample
            % io부터 ofdm심볼갯수*80 -1 (소숫점으로나눠지는샘플은 취급x)
            i1 = i0 + obj.nsampSym*nsym-1;
            
            % Get the samples corresponding to the main packet
            r1 = r(i0:i1);

            % TODO:  Per    form the OFDM demod on r1 and extract the active
            % FFT indices and store the raw symbols in symRx
                symRx = ofdmdemod(r1,obj.info.FFTLength,obj.info.CPLength);
                symRx = symRx(obj.info.ActiveFFTIndices,:);
                % 52* nsym / 52*1

            % Equalize the symbols 
            if strcmp(obj.eqMethod, 'inversion')
                % TODO:  Inversion method
                %    obj.symEq = ...
                % 52*nsym / 52*1 = 52*nsym
                obj.symEq=symRx ./obj.chanEst;

            elseif strcmp(obj.eqMethod, 'MMSE')
                % TODO:  MMSE method
                %    obj.symEq = ...
                % 라는데? MMSE:  symEq(i,j) = conj(chanEst(i))*sym(i,j) / (abs(chanEst(i)).^2 + noiseEst)
                % 정해진 해가 있네 유도는몰라
                % 채널= conj(추정채널)*Es /abs(추정).^2 *E+ 노이즈분산(노이즈평균전력)
                % E=1일때 간단히인데.. 송신 심볼의 E는1로 가정하나봐 ->그냥써도된다고?..
                % 요소별 곱 ㅇㅇ
                obj.symEq= (conj(obj.chanEst) ./ (abs(obj.chanEst).^2 + obj.noiseEst)) .* symRx;
            else
                error('Unknown equalization method');
            end

            % TODO:  Extract the equalized data symbols from symEq using the 
            % obj.info.DataIndices
            %    obj.symDataEq = obj.symEq(...)
            obj.symDataEq = obj.symEq(obj.info.DataIndices,:); % 48 * d
        end

    end

    methods (Access = protected)

        function stepImpl(obj, r, indltf)
            % Channel estimation and equalization
            obj.compChanEst(r, indltf);
            obj.eqSym(r,indltf);

        end


    end

end