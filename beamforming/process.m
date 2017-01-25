background=zeros(14,1);
background=mean(abs(ifft(channel_log(1:20,:))));
for ii=30:10000
    plot(max(mean(abs(ifft(channel_log(ii-20:ii,:))))-background,0));
    pause(0.1);
end