802.11g OFDM수신기(복조+이퀄라이즈) 만들기 



802.11g Waveform -> OFDM 패킷추출 -> 랜덤 다중경로채널 적용 ->채널 추정 -> Symbol Equalization ->  SDR 송수신 


Lab8: 메인 스크립트 \
WLANPktAnalyzer.m    : OFDM패킷 추출 \
WLANRx.m            :  채널추정+이퀄라이제이션 \
RandMPChan.m        : 랜덤 다중 경로 채널 \
WLANTx:		:  802.11 Waveform 생성 
plutoCreateTxR  : SDR송수신 객체 생성
