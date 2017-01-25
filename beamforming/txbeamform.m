
close all;
clear cam1 cam2
cam1=webcam(1);
cam2=webcam(2);
if exist('nodes','var')
    reinit_nodes=false;
else
    reinit_nodes=true;
end
tic
USE_AGC = false;
USE_EXTERNAL_TRIGGER=false;
NUM_RX=4;
freq=2462;
num_chans=14;
slow_down=false;
ic_enable=false;
gain_boost=false;
chronos=false;
manual_agc=false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the WARPLab experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NUMNODES = 2;band=2.4;
fft_angles_log=[];cal_angles=zeros(1,4);null_factor=0;
ym_log1=[];ym_log2=[];xm_log1=[];xm_log2=[];t_log=[];
fft_rf_csi_log=[];phase_log=[];amp_log=[];fft_log_rx1=[];fft_log_rx2=[];fft_log_rx3=[];fft_log_rx4=[];music_inp_log=[];
fft_tx1=[];fft_tx2=[];fft_tx3=[];fft_tx4=[];fft_rx=[];fft_tx=[];

lambda_log=[];Z_log=[];

%Create a vector of node objects
if(reinit_nodes)
    nodes = wl_initNodes(NUMNODES);
end

%Create a UDP broadcast trigger and tell each node to be ready for it
eth_trig = wl_trigger_eth_udp_broadcast;
wl_triggerManagerCmd(nodes,'add_ethernet_trigger',[eth_trig]);
[T_IN_ETH,T_IN_ENERGY,T_IN_AGCDONE,T_IN_REG,T_IN_D0,T_IN_D1,T_IN_D2,T_IN_D3] =  wl_getTriggerInputIDs(nodes(1));
[T_OUT_BASEBAND, T_OUT_AGC, T_OUT_D0, T_OUT_D1, T_OUT_D2, T_OUT_D3] = wl_getTriggerOutputIDs(nodes(1));

%For the transmit node, we will allow Ethernet to trigger the buffer
%baseband, the AGC, and debug0 (which is mapped to pin 8 on the debug
%header)
nodes(1).wl_triggerManagerCmd('output_config_input_selection',[T_OUT_BASEBAND,T_OUT_D0],[T_IN_ETH,T_IN_REG]);

if(USE_EXTERNAL_TRIGGER)
    %For the receive node, we will allow debug3 (mapped to pin 15 on the
    %debug header) to trigger the buffer baseband, and the AGC
    nodes(2).wl_triggerManagerCmd('output_config_input_selection',[T_OUT_BASEBAND,T_OUT_AGC],[T_IN_D3],[T_IN_ETH,T_IN_REG]);
else
    nodes(2).wl_triggerManagerCmd('output_config_input_selection',[T_OUT_BASEBAND,T_OUT_AGC],[T_IN_ETH,T_IN_REG]);
end

%%

%Get IDs for the interfaces on the boards. Since this example assumes each
%board has the same interface capabilities, we only need to get the IDs
%from one of the boards
[RFA,RFB,RFC,RFD] = wl_getInterfaceIDs(nodes(1));


wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'tx_gains',3,23);
wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'channel',2.4,14);
% wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'tx_lpf_corn_freq',3);
% wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'rx_lpf_corn_freq',3);

if(USE_AGC)
    wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'rx_gain_mode','automatic');
    wl_basebandCmd(nodes,'agc_target',-6);
    wl_basebandCmd(nodes,'agc_trig_delay', 500);
    wl_basebandCmd(nodes,'agc_dco', true);
else
    wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'rx_gain_mode','manual');
    RxGainRF = 3; %Rx RF Gain in [1:3]
    RxGainBB = 5;%Rx Baseband Gain in [0:31]
    
    wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'rx_gains',RxGainRF,RxGainBB+2);
    wl_interfaceCmd(nodes,RFA,'rx_gains',0,18);
end


%We'll use the transmitter's I/Q buffer size to determine how long our
%transmission can be
txLength = nodes(1).baseband.txIQLen;

%Set up the baseband for the experiment
wl_basebandCmd(nodes,'tx_delay',0);
wl_basebandCmd(nodes,'tx_length',txLength);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal processing to generate transmit signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First generate the preamble for AGC. The preamble corresponds to the
% short symbols from the 802.11a PHY standard
shortSymbol_freq = [0 0 0 0 0 0 0 0 1+i 0 0 0 -1+i 0 0 0 -1-i 0 0 0 1-i 0 0 0 -1-i 0 0 0 1-i 0 0 0 0 0 0 0 1-i 0 0 0 -1-i 0 0 0 1-i 0 0 0 -1-i 0 0 0 -1+i 0 0 0 1+i 0 0 0 0 0 0 0].';
shortSymbol_freq = [zeros(32,1);shortSymbol_freq;zeros(32,1)];
shortSymbol_time = ifft(fftshift(shortSymbol_freq));
shortSymbol_time = (shortSymbol_time(1:32).')./max(abs(shortSymbol_time));
shortsyms_rep = repmat(shortSymbol_time,1,32);

preamble = shortsyms_rep;
preamble = preamble(:);

Ts = 1/(wl_basebandCmd(nodes(1),'tx_buff_clk_freq'));
t = [0:Ts:((txLength-length(preamble)-1))*Ts].'; % Create time vector(Sample Frequency is Ts (Hz))

node_tx = nodes(1);
node_rx = nodes(2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit and receive signal using WARPLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rx_gains_log=[];
fft_rf_csi_log=[];

h1_ic=zeros(1,num_chans);h2_ic=zeros(1,num_chans);
cal_obs=50;
cal_ii=cal_obs*num_chans;
null_factor=0;
wl_interfaceCmd(node_tx,RFA+RFB+RFC+RFD,'tx_rx_dis');

wl_interfaceCmd(node_rx,RFA+RFB+RFC+RFD,'rx_en');


wl_basebandCmd(node_rx,RFA+RFB+RFC+RFD,'rx_buff_en');
lts_f=[repmat([1 0 -1 0],1,32) ];
lts_f(:)=0;
lts_f([2 128])=1;
lts_t=ifft(lts_f);
payload_B=lts_t./max([real(lts_t) imag(lts_t)]);
payload_B=repmat(payload_B.',248,1);
txData_B = [preamble;payload_B];
wl_interfaceCmd(node_tx,RFA+RFB+RFC+RFD,'tx_en');
wl_basebandCmd(node_tx,RFA+RFB+RFC+RFD,'tx_buff_en');
txData = [txData_B];
for kk=1:3
    
    for ii=1:num_chans
        wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'channel',band,ii);
        wl_basebandCmd(node_tx,[RFA RFB RFC RFD], 'write_IQ', [txData_B 0*txData_B 0*txData_B 0*txData_B]);
        eth_trig.send();
        rx_IQ = wl_basebandCmd(node_rx,[RFA RFB RFC RFD],'read_IQ', 0, 10000);
        fft_rf=fft(rx_IQ(3500+1+1280+1280+0*128:3500+1280+1280+1*128,:))+fft(rx_IQ(3500+1+1280+1280+1*128:3500+1280+1280+2*128,:))+fft(rx_IQ(3500+1+1280+1280+2*128:3500+1280+1280+3*128,:))+fft(rx_IQ(3500+1+1280+1280+4*128:3500+1280+1280+5*128,:))+fft(rx_IQ(3500+1+1280+1280+5*128:3500+1280+1280+6*128,:));
        fft_rf_csi=fft_rf(:,:).*repmat(lts_f.',1,size(fft_rf,2));
        fft_tx1=[fft_tx1 (fft_rf_csi(2,:)+fft_rf_csi(128,:)).']
        if(manual_agc)
            figure(345)
            plot(abs(rx_IQ))
            pause
        end
        
        wl_basebandCmd(node_tx,[RFA RFB RFC RFD], 'write_IQ', [1*txData_B txData_B 0*txData_B 0*txData_B]);
        eth_trig.send();
        rx_IQ = wl_basebandCmd(node_rx,[RFA RFB RFC RFD],'read_IQ', 0, 10000);
        fft_rf=fft(rx_IQ(3500+1+1280+1280+0*128:3500+1280+1280+1*128,:))+fft(rx_IQ(3500+1+1280+1280+1*128:3500+1280+1280+2*128,:))+fft(rx_IQ(3500+1+1280+1280+2*128:3500+1280+1280+3*128,:))+fft(rx_IQ(3500+1+1280+1280+4*128:3500+1280+1280+5*128,:))+fft(rx_IQ(3500+1+1280+1280+5*128:3500+1280+1280+6*128,:));
        fft_rf_csi=fft_rf(:,:).*repmat(lts_f.',1,size(fft_rf,2));
        fft_tx2=[fft_tx2 (fft_rf_csi(2,:)+fft_rf_csi(128,:)).']
        if(manual_agc)
            figure(345)
            plot(abs(rx_IQ))
            pause
        end
        
        wl_basebandCmd(node_tx,[RFA RFB RFC RFD], 'write_IQ', [1*txData_B 0*txData_B 1*txData_B 0*txData_B]);
        eth_trig.send();
        rx_IQ = wl_basebandCmd(node_rx,[RFA RFB RFC RFD],'read_IQ', 0, 10000);
        fft_rf=fft(rx_IQ(3500+1+1280+1280+0*128:3500+1280+1280+1*128,:))+fft(rx_IQ(3500+1+1280+1280+1*128:3500+1280+1280+2*128,:))+fft(rx_IQ(3500+1+1280+1280+2*128:3500+1280+1280+3*128,:))+fft(rx_IQ(3500+1+1280+1280+4*128:3500+1280+1280+5*128,:))+fft(rx_IQ(3500+1+1280+1280+5*128:3500+1280+1280+6*128,:));
        fft_rf_csi=fft_rf(:,:).*repmat(lts_f.',1,size(fft_rf,2));
        fft_tx3=[fft_tx3 (fft_rf_csi(2,:)+fft_rf_csi(128,:)).']
        if(manual_agc)
            figure(345)
            plot(abs(rx_IQ))
            pause
        end
        
        wl_basebandCmd(node_tx,[RFA RFB RFC RFD], 'write_IQ', [1*txData_B 0*txData_B 0*txData_B 1*txData_B]);
        eth_trig.send();
        rx_IQ = wl_basebandCmd(node_rx,[RFA RFB RFC RFD],'read_IQ', 0, 10000);
        fft_rf=fft(rx_IQ(3500+1+1280+1280+0*128:3500+1280+1280+1*128,:))+fft(rx_IQ(3500+1+1280+1280+1*128:3500+1280+1280+2*128,:))+fft(rx_IQ(3500+1+1280+1280+2*128:3500+1280+1280+3*128,:))+fft(rx_IQ(3500+1+1280+1280+4*128:3500+1280+1280+5*128,:))+fft(rx_IQ(3500+1+1280+1280+5*128:3500+1280+1280+6*128,:));
        fft_rf_csi=fft_rf(:,:).*repmat(lts_f.',1,size(fft_rf,2));
        fft_tx4=[fft_tx4 (fft_rf_csi(2,:)+fft_rf_csi(128,:)).']
        if(manual_agc)
            figure(345)
            plot(abs(rx_IQ))
            pause
        end
    end
    fft_tx=[fft_tx1(2,:) ; fft_tx2(2,:) ; fft_tx3(2,:); fft_tx4(2,:)]; %4 x num_chans
    wl_basebandCmd(node_tx,[RFA RFB RFC RFD], 'write_IQ', [1*txData_B 1*txData_B 0*txData_B 0*txData_B]);
    for ii=1:num_chans
        wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'channel',band,ii);
        eth_trig.send();
        rx_IQ = wl_basebandCmd(node_rx,[RFA RFB RFC RFD],'read_IQ', 0, 10000);
        fft_rf=fft(rx_IQ(3500+1+1280+1280+0*128:3500+1280+1280+1*128,:))+fft(rx_IQ(3500+1+1280+1280+1*128:3500+1280+1280+2*128,:))+fft(rx_IQ(3500+1+1280+1280+2*128:3500+1280+1280+3*128,:))+fft(rx_IQ(3500+1+1280+1280+4*128:3500+1280+1280+5*128,:))+fft(rx_IQ(3500+1+1280+1280+5*128:3500+1280+1280+6*128,:));
        fft_rf_csi=fft_rf(:,:).*repmat(lts_f.',1,size(fft_rf,2));
        fft_rx=[fft_rx (fft_rf_csi(2,:)+fft_rf_csi(128,:)).'];%4 x num_chans
        
    end
    
end
pause
wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'tx_gains',3,30);
wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'rx_gains',RxGainRF,RxGainBB+6);
wl_interfaceCmd(nodes,RFA,'rx_gains',0,18);
peak_angle=0;
for ii=0:num_chans*24*150-1
%     if(mod(ii,num_chans*24*30)==0)
%         'pausing for stimulus'
%         pause
%     end
    if(mod(ii,num_chans)==0)
    peak_angle=peak_angle+15;
    
    d=0.0625;
    wl=3e8/2.4e9;
    phase_diff=repmat(exp(j*2*pi*d*sind(peak_angle)*[0:3]./wl),size(txData_B,1),1);
    capture2d;
    ii
    end
    chan_no=mod(ii,num_chans)+1;
    wl_basebandCmd(node_tx,[RFA RFB RFC RFD], 'write_IQ', [txData_B txData_B./exp(j*angle(fft_tx(2,chan_no))) txData_B./exp(j*angle(fft_tx(3,chan_no))) txData_B./exp(j*angle(fft_tx(4,chan_no)))].*phase_diff);
    wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'channel',band,chan_no);
    eth_trig.send();
    rx_IQ = wl_basebandCmd(node_rx,[RFA RFB RFC RFD],'read_IQ', 0, 10000);
    fft_rf=fft(rx_IQ(3500+1+1280+1280+0*128:3500+1280+1280+1*128,:))+fft(rx_IQ(3500+1+1280+1280+1*128:3500+1280+1280+2*128,:))+fft(rx_IQ(3500+1+1280+1280+2*128:3500+1280+1280+3*128,:))+fft(rx_IQ(3500+1+1280+1280+4*128:3500+1280+1280+5*128,:))+fft(rx_IQ(3500+1+1280+1280+5*128:3500+1280+1280+6*128,:));
    fft_rf_csi=fft_rf(:,:).*repmat(lts_f.',1,size(fft_rf,2));
    fft_rf_csi_dc=fft_rf_csi(128,:)+fft_rf_csi(2,:);
    fft_rf_csi_log=[fft_rf_csi_log fft_rf_csi_dc.'];
    if(slow_down)
        figure(1205)
        plot(abs(rx_IQ))
    end
   
end


wl_basebandCmd(nodes,RFA+RFB+RFC+RFD,'tx_rx_buff_dis');
wl_interfaceCmd(nodes,RFA+RFB+RFC+RFD,'tx_rx_dis');

