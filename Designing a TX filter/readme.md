16QAM심볼생성 -> TxFilter적용 -> 필터적용& 수신 신호의 PSD확인 \
 \
 결과
 ![](lab2result1.jpg) \
 SDR 송수신
 ![](lab2result.gif)
 \
Unit02_lab.m: 메인스크립트 \
TxFilt.m: Waveform 오버샘플링&패스밴드&스탑밴드 필터디자인 \
plutoCreateTxRx.m : SDR 송수샌 객체설정
