802.11g Wi-fi 모델 패킷안에있는 프리엠블탐지하기 (WLAN toolbox 설치) \
1. auto-corrlation 으로 STF감지 \
2. matched-filter detector로 LTF감지 \
\
Lab8: 메인 스크립트
WLANTx:  802.11 Waveform 생성
RandDelayChan.m: 1경로지연 및 SNR 채널 적용 
WLANRx.m:  STF 및 LTF감지 
TxFilt.m: 업샘플링  & 대역폭제한 필터
RxFilt.m: 다운샘플링 & 대역폭제한 필터
plutoCreateTxR  : SDR송수신 객체 생성
