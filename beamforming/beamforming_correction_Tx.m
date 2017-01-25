function beamforming_correction_Tx(tx_n)

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
    node_rx = nodes(2);
    %Create a UDP broadcast trigger and tell each node to be ready for it
    eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(nodes,'add_ethernet_trigger',[eth_trig]);
    
    trig_in_ids  = wl_getTriggerInputIDs(node_tx);
    trig_out_ids = wl_getTriggerOutputIDs(node_tx);

    % For the transmit node, we will allow Ethernet to trigger the buffer baseband, the AGC, and debug0 
    % (which is mapped to pin 8 on the debug header)
    node_tx.wl_triggerManagerCmd('output_config_input_selection', [trig_out_ids.BASEBAND, trig_out_ids.EXT_OUT_P0], [trig_in_ids.ETH_A]);

    node_rx.wl_triggerManagerCmd('output_config_input_selection', [trig_out_ids.BASEBAND, trig_out_ids.AGC], [trig_in_ids.EXT_IN_P0, trig_in_ids.EXT_IN_P3]);

    % For the receive node, we enable the debounce circuity on the debug 3 input
    % to deal with the fact that the signal may be noisy.
    node_rx.wl_triggerManagerCmd('input_config_debounce_mode', [trig_in_ids.EXT_IN_P0, trig_in_ids.EXT_IN_P3], true); 

    % Since the debounce circuitry is enabled, there will be a delay at the
    % receiver node for its input trigger. To better align the transmitter and
    % receiver, we can artifically delay the transmitters trigger outputs that
    % drive the buffer baseband and the AGC.
    node_tx.wl_triggerManagerCmd('output_config_delay', [trig_out_ids.EXT_OUT_P0], 0);
    node_tx.wl_triggerManagerCmd('output_config_delay', [trig_out_ids.BASEBAND], 62.5);     % 62.5ns delay

    [RFA,RFB, RFC,RFD] = wl_getInterfaceIDs(nodes(1));
    [RFE,RFF] = wl_getInterfaceIDs(nodes(2));
     %=wl_getInterfaceIDs(nodes(1));

    %Set up the interface for the experiment and gain
    wl_interfaceCmd(nodes,'RF_ALL','tx_gains',3,20);
	wl_interfaceCmd(nodes,'RF_ALL','rx_gain_mode','manual');
	wl_interfaceCmd(nodes,'RF_ALL','rx_gains',3,15);

    TX_NUM_SAMPS = 1000;
    node_tx = nodes(1);
    node_rx = nodes(2);
    switch tx_n
        case 1
            RF_TX=RFA;
        case 2
            RF_TX=RFB;
        case 3
            RF_TX=RFC;
        case 4
            RF_TX=RFD;
    end
    RF_RX = RFE;
    
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
    wl_interfaceCmd(node_rx,RF_RX,'rx_en');

    %Enable the Tx and Rx buffers
    wl_basebandCmd(node_tx,RF_TX,'tx_buff_en');
    wl_basebandCmd(node_rx,RF_RX,'rx_buff_en');
    hold on;
    for channel_num=1:1
        phase_tmp_log=zeros(20,1);
        measure_num=1;
        wl_interfaceCmd(nodes,'RF_ALL','channel',2.4,14);
        while measure_num<=20
        %tic
        %Trigger the Tx/Rx cycle at both nodes
        eth_trig.send();

        %Retrieve the received waveform from the Rx node
        rx_vec_air = wl_basebandCmd(node_rx,[RF_RX],'read_IQ', 0, TX_NUM_SAMPS);
        rx_vec_air = rx_vec_air(:).';

        
        lts_corr = abs(conv(conj(fliplr(lts_t)), sign(rx_vec_air)));

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
        rx_lts = rx_vec_air([1:128]);
        rx_lts1 = rx_lts([1:64]);
        rx_lts2 = rx_lts([65:128]);

        rx_lts1_f = fft(rx_lts1);
        rx_lts2_f = fft(rx_lts2);

        %Calculate channel estimate
        rx_H_est = lts_f .* (rx_lts1_f + rx_lts2_f)/2;
        phase_tmp_log(measure_num)=exp(1i*angle((rx_H_est(1)+rx_H_est(64))/2));
        figure(1);
        plot(1:measure_num,angle(phase_tmp_log(1:measure_num)));
        
        figure(2);
        subplot(2,1,1);
        plot(real(rx_vec_air), 'b');
        axis([0 length(rx_vec_air) -TX_SCALE TX_SCALE])
        grid on;
        title('Rx Waveform (I)');

        subplot(2,1,2);
        plot(imag(rx_vec_air), 'r');
        axis([0 length(rx_vec_air) -TX_SCALE TX_SCALE])
        grid on;
        title('Rx Waveform (Q)');
      
        pause(time_window);
        %toc
        measure_num=measure_num+1;
        end
        phase_offset_log(tx_n,channel_num)=mean(phase_tmp_log);
        figure(3);
        plot(1:channel_num,angle(phase_offset_log(tx_n,1:channel_num)));
    end
    %phase_offset(tx_n)=mean(phase_offset_log);
    %Disable the Tx/Rx radios and buffers
    wl_basebandCmd(nodes,'RF_ALL','tx_rx_buff_dis');
    wl_interfaceCmd(nodes,'RF_ALL','tx_rx_dis');
end