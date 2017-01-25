function batch_test()

    time_window=0.0001;
	%Waveform params
	TX_SCALE = 1.0;     %Sale for Tx waveform ([0:1])

    NUMNODES = 2;

    %Create a vector of node objects
    global nodes;
    global phase_offset_log;
    if(isempty(nodes))
        nodes = wl_initNodes(NUMNODES);
    end
    if(isempty(phase_offset_log))
        phase_offset_log=zeros(8,14);
    end
    node_tx = nodes(1);

    %Create a UDP broadcast trigger and tell each node to be ready for it
    eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(nodes,'add_ethernet_trigger',[eth_trig]);
    
    trig_in_ids  = wl_getTriggerInputIDs(node_tx);
    trig_out_ids = wl_getTriggerOutputIDs(node_tx);

    [RFA RFB] = wl_getInterfaceIDs(nodes(1));
	[RFC RFD RFE RFF] = wl_getInterfaceIDs(nodes(2));

    %Set up the interface for the experiment and gain
    wl_interfaceCmd(nodes,'RF_ALL','tx_gains',3,20);
	wl_interfaceCmd(nodes,'RF_ALL','rx_gain_mode','manual');
	wl_interfaceCmd(nodes,'RF_ALL','rx_gains',3,15);

    TX_NUM_SAMPS = 1000;
    node_tx = nodes(1);
    node_rx = nodes(2);
    
    RF_TX=RFA;
    RF_RX =[RFC RFD RFE RFF];
	
	% Precalculation for implementing Wision
	num_ant=4;
	%antenna_num=0:num_ant-1;
	d=0.06; %distances between antennas
	lambda=0.12; %full wave length
	angles=0:pi/12:2*pi-pi/12;
	% Base function used in Wision
	basis=exp(-1j* repmat([0:num_ant-1],24,1).*repmat(cos(angles'),1,num_ant) * 2*pi*d / lambda);
    
    %Set up the baseband for the experiment
    wl_basebandCmd(nodes,'tx_delay',0);
    wl_basebandCmd(nodes,'tx_length',TX_NUM_SAMPS); 

	%LTS for CFO and channel estimation
	lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
	lts_t = ifft(lts_f, 64);

	%Use 30 copies of the 16-sample STS for extra AGC settling margin
	preamble = [lts_t lts_t];
    payload=lts_t./max([real(lts_t) imag(lts_t)]);

	%Construct the full time-domain OFDM waveform
	tx_vec = [preamble payload];

	%Pad with zeros for transmission
	tx_vec_air = [tx_vec zeros(1,TX_NUM_SAMPS-length(tx_vec))];

	%Scale the Tx vector
	tx_vec_air = TX_SCALE .* tx_vec_air ./ max(abs(tx_vec_air));
	
	%% WARPLab Tx/Rx 
	%Write the Tx waveform to the Tx node
	wl_basebandCmd(node_tx,[RF_TX], 'write_IQ', tx_vec_air(:));
    %Enable the Tx and Rx radios
    wl_interfaceCmd(node_tx,RF_TX,'tx_en');
    wl_interfaceCmd(node_rx,RFC+RFD+RFE+RFF,'rx_en');

    %Enable the Tx and Rx buffers
    wl_basebandCmd(node_tx,RF_TX,'tx_buff_en');
    wl_basebandCmd(node_rx,RFC+RFD+RFE+RFF,'rx_buff_en');
    %hold on;

	measure_num=1;
	wl_interfaceCmd(nodes,'RF_ALL','channel',5,14);
	while measure_num<=2000
        %tic
        %Trigger the Tx/Rx cycle at both nodes
        eth_trig.send();

        %Retrieve the received waveform from the Rx node
        rx_mat_air = wl_basebandCmd(node_rx,[RF_RX],'read_IQ', 0, TX_NUM_SAMPS);
		
		rx_vec_air_A = rx_mat_air(:,1).';
		rx_vec_air_B = rx_mat_air(:,2).';
		rx_vec_air_C = rx_mat_air(:,3).';
		rx_vec_air_D = rx_mat_air(:,4).';

        
        lts_corr = abs(conv(conj(fliplr(lts_t)), sign(rx_vec_air_A)));

        %Skip early and late samples
        lts_corr = lts_corr(32:end-32);

        %Find all correlation peaks
        lts_peaks = find(lts_corr > 0.8*max(lts_corr));

        %Select best candidate correlation peak as LTS-payload boundary
        [LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
        [lts_second_peak_index,y] = find(LTS2-LTS1 == length(lts_t));

        %Punt if no valid correlation peak was found
        if(isempty(lts_second_peak_index))
            fprintf('No LTS Correlation Peaks Found!\n');
            continue;
        end
        
        %Re-extract LTS for channel estimate
        rx_lts_A = rx_vec_air_A([1:128]);
		rx_lts_B = rx_vec_air_B([1:128]);
		rx_lts_C = rx_vec_air_C([1:128]);
		rx_lts_D = rx_vec_air_D([1:128]);
        %Calculate channel estimate
        rx_H_est_A = lts_f .* ( fft(rx_lts_A([1:64])) + fft(rx_lts_A([65:128])))/2;
		rx_H_est_B = lts_f .* ( fft(rx_lts_B([1:64])) + fft(rx_lts_B([65:128])))/2;
		rx_H_est_C = lts_f .* ( fft(rx_lts_C([1:64])) + fft(rx_lts_C([65:128])))/2;
		rx_H_est_D = lts_f .* ( fft(rx_lts_D([1:64])) + fft(rx_lts_D([65:128])))/2;
		
		csi_A=(rx_H_est_A(1)+rx_H_est_A(64))/2;
		csi_B=(rx_H_est_B(1)+rx_H_est_B(64))/2;
		csi_C=(rx_H_est_C(1)+rx_H_est_C(64))/2;
		csi_D=(rx_H_est_D(1)+rx_H_est_D(64))/2;
		
		figure(1);
		wision_ans=abs(sum((repmat([csi_A.' csi_B.' csi_C.' csi_D.'],24,1).*basis)'));
        polar(angles,wision_ans);
        
        figure(2);
        subplot(2,1,1);
        plot(real(rx_vec_air_A), 'b');
        axis([0 length(rx_vec_air_A) -TX_SCALE TX_SCALE])
        grid on;
        title('Rx Waveform (I)');

        subplot(2,1,2);
        plot(imag(rx_vec_air_A), 'r');
        axis([0 length(rx_vec_air_A) -TX_SCALE TX_SCALE])
        grid on;
        title('Rx Waveform (Q)');
      
        pause(time_window);
        %toc
        measure_num=measure_num+1;
    end
    %phase_offset(tx_n)=mean(phase_offset_log);
    %Disable the Tx/Rx radios and buffers
    wl_basebandCmd(nodes,'RF_ALL','tx_rx_buff_dis');
    wl_interfaceCmd(nodes,'RF_ALL','tx_rx_dis');
end