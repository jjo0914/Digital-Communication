classdef WLANPktAnalyzer <  matlab.System
    % Class for demodulating a transmitted waveform without noise
    % We use this class for extracting the true transmitted symbols
    % and other parameters

    properties
        info;  % OFDM info
        nsampSym;  % number of samples per OFDM symbol

        % Hard-coded constants based on the WLAN spec
        nsymLSTF = 2;  % number symbols in the legacy STF field
        nsymLLTF = 2;  % number symbols in the legacy LTF field
        nsymPre = 5;   % total number of symbols in pre-amble

        % Indices for the pilots and data within the payload
        Idata;
        Ipilot;

        % Extracted fields
        sym;     % OFDM modulation symbols over preamble and payload
        symLTF;  % L-LTF OFDM modulation symbols
        symData; % Data symbols
        symPilot;  % Pilots symbols
        

    end

    methods
        function obj = WLANPktAnalyzer()
            % Constructor            

            % Get the OFDM info
            obj.info = wlanNonHTOFDMInfo('NonHT-Data'); % 'NonHT': NonHT=802.11a.g 'Data ..LTF.. Header'
                                                        
                                                

            % TODO:  Using the FFTLength and CPLength fields in obj.info,
            % get the number of samples per OFDM symbol
            %    obj.nsampSym = ...
            obj.nsampSym = obj.info.FFTLength + obj.info.CPLength;

 
        end

        function extractLTFSym(obj, x)
            % TODO:  Use the ofdmdemod function to extract the OFDM symbols
            % from x.  You can get the FFTLength and CPLength from
            % obj.info
            %    obj.sym = ofdmdemod(...); % 필요한거 심볼fft길이,fft길이,cp길이,, fft윈도우를 cp기준으로 어디인지
            sym=ofdmdemod(x,obj.info.FFTLength,obj.info.CPLength); % 64*48등장
           
            % TODO: Extract the symbols on the active sub-carriers with
            % the obj.info.ActiveFFTIndices field.
            %    obj.sym = obj.sym(...,:);
            obj.sym=sym(obj.info.ActiveFFTIndices,:); % % 4096비트-> mcs1 = 172심볼 ,mcs4=43 ofdm심볼  
                                                  % STF(8us)+ltf 2개(8us) +
                                                  % L_sig1개 (4us) + 43=48개
                                                  % 64* 48 맞음
                                                  % 3 4 가 LTF 맞네 ㅇㅇ
                                                      

            % TODO:  Get the symbols for the LTF field. Recall that the LTF
            % field is on symbols 3 and 4 (the two symbols after the STF)
            %    obj.symLTF = obj.sym(...);
            %? data기준으로 복조했는데?
            obj.symLTF=sym(obj.info.ActiveFFTIndices,[3 4]);

        end

        function extractPayloadSym(obj)

            % TODO:  From the array obj.sym, Get the modulations symbols for the
            % payload.  The payload starts at obj.nsymPre+1
            %    symPayload = obj.sym(...)
            symPayload=obj.sym(:,obj.nsymPre+1:end)
            % TODO:  Get the symbols for pilots and data using the
            % obj.info.PilotIndices and obj.info.DataIndices arrays
            %    obj.symPilot = symPayload(...);
            %    obj.symData = symPayload(...);
            obj.symPilot = symPayload(obj.info.PilotIndices,:);
             obj.symData = symPayload(obj.info.DataIndices,:);
            

        end
    end


    methods (Access = protected)
        function stepImpl(obj,x)
            % Step function:  Extact the OFDM modulation symbols for the
            % LTF and payload fields

            % Extract the LTF symbols
            obj.extractLTFSym(x);

            % Extract the payload symbols
            obj.extractPayloadSym();

        end
    end
end
