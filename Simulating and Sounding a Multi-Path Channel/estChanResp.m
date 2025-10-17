function [hfd,h] = estChanResp(r,xfd,opt)
% estChanResp:  Estimates the channel response using frequency-domain
% correlation
% 
% The model is that the TX repeatedly sends x=ifft(xfd) and the RX 
% receives samples
%     r = h*x + w
% where h is the channel impulse response and w is noise.  
%
% Parameters
% ----------
% r:  RX complex baseband samples of length nfft
% xfd:  TX samples in frequency-domain
% nleft:  Number of samples to the left to place peak at
% normToMean:  Boolean flag indicating that the energy / tap
%     should be normalized to the average energy
%stimate
% Returns
% -------
% hfd:   Frequency domain channel estimate.  hfd(k) is the e
%    of the channel frequency response at f(k) = k/nfft*fsamp
% h:  Time-domain channel impulse response estimate.
arguments
    r (:,1) double;
    xfd (:,1) double;
    
    opt.nleft (1,1) {mustBeInteger} = 8; 
    opt.nright (1,1) {mustBeInteger} = 8;
    opt.normToMean (1,1) = false;

end

% TODO:  Take the FFT of r 
%    rfd = ...
rfd=fft(r);

% TODO:  Estimate the channel frequency response by dividing by the 
% FFT of x, xfd
%    hfd = ...
hfd=rfd./xfd;


% TODO:  Compute the time-domain response
%    h = ...
h=ifft(hfd);

% TODO:  Find the peak location in h
% Then, circularly shift h so that the peak is at position opt.nleft
%    h = circshift(...)
%circshift는 배열을 원형배열로 취급하고 위치를 변경할수있듬
[~,idx]=max(h);
h=circshift(h,opt.nleft-idx); %문제되로 왼쪽8번째위치에위치
                           %(배열,위치,dim) dim 1이면왼쪽이동 2면오른쪽이동
% figure(11);
% plot(abs(h));
% hold on;
% legend('노말라이즈x')


% Normalize to the average energy signal level
%정규화하는 이유는 몰라?유용하대..
%아 이게 r수신할때 스케일링 하냐 안하냐 에 따라 크기가다르거든? 
%근데 true로하면 r의 스케일링이 얼마인지 상관없어 <
if opt.normToMean
    Emean = mean(abs(h).^2);
    h = h / sqrt(Emean);
    % plot(abs(h));
end

end