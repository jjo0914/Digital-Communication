classdef  WLANTx < matlab.System
    % Wireless LAN transmitter class
    properties
        % Packet properties
        psduLen = 512;   % PSDU length in Bytes
        mcs = 3;         % MCS selection
        cfg;             % WLAN Toolbox configuration
        fsamp;           % Sample rate in Hz
    end
    
    methods
        function obj = WLANTx(opt)
            % Constructor
            arguments
                opt.psduLen (1,1) {mustBeInteger} = 512;
                opt.mcs (1,1) {mustBeInteger} = 3;
            end
           
            % Set properties
            obj.psduLen = opt.psduLen;
            obj.mcs = opt.mcs;
                        
            % ----------802.11패킷설정--------------------
            % 802.11패킷 설정요소
            % 기본:OFDM
            % 채널폭:802.11 표준에따라 설정
            % PSDULength:바이트
            % MCS: 1(BPSK-3/4) 2(QPSK-1/2) 3(QPSK3/4)
    %code rate: 유효 서브캐리어 48=48심볼=96비트 여기서 코딩레이트가 1/2면 48비트가 유효데이터 
           %  1OFDM=N_CP +가드밴드+파일럿서브캐리어+유효서브캐리어+가드밴드
            obj.cfg = wlanNonHTConfig('ChannelBandwidth', 'CBW20', ...
                'PSDULength', obj.psduLen,'MCS',obj.mcs);            

         % ----WLAN 채널 BandWidth('CBW1~120')에 맞는 표준 샘플레이트를 가져옵니다.            
             obj.fsamp = wlanSampleRate(obj.cfg)  ;          
                        
        end
    end
                
    methods (Access = protected)
        function x = stepImpl(obj)
            % ----WLAN 패킷만들기----------
            %  psdeLen길이에 맞춰서 랜덤비트 만들기 (1바이트=8비트)
             bits = randi([0 1],obj.psduLen*8,1); % 4096비트-> mcs1 = 172심볼 ,mcs4=43 ofdm심볼  
                                                  % STF(8us)+ltf 2개(8us) +
                                                  % L_sig1개 (4us) + 43=48개
                                                 
          
            % ---함수! wlanWaveformGenrater 사용하기 using bits and obj.cfg
            x = wlanWaveformGenerator(bits,obj.cfg); %그냥 넘어가자
           
        end              
    end
end