%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wl_example_siso_ofdm_txrx.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%Params:
USE_WARPLAB_TXRX = 1;   %Enable WARPLab-in-the-loop (otherwise sim-only)
WRITE_PNG_FILES = 0;    %Enable writing plots to PNG

%Waveform params
N_OFDM_SYMS = 190;  %Number of OFDM symbols
MOD_ORDER = 16;     %Modulation order (1/4/16 = BSPK/QPSK/16-QAM)
TX_SCALE = 1.0;     %Sale for Tx waveform ([0:1])
INTERP_RATE = 2;        %Interpolation rate (1 or 2)

%OFDM params
SC_IND_PILOTS = [8 22 44 58];   %Pilot subcarrier indices
SC_IND_DATA   = [2:7 9:21 23:27 39:43 45:57 59:64]; %Data subcarrier indices
N_SC = 64;          %Number of subcarriers
CP_LEN = 16;        %Cyclic prefix length
N_DATA_SYMS = N_OFDM_SYMS * length(SC_IND_DATA); %Number of data symbols (one per data-bearing subcarrier per OFDM symbol)

%Rx processing params
FFT_OFFSET = 4;     %Number of CP samples to use in FFT (on average)
LTS_CORR_THRESH = 0.8;  %Normalized threshold for LTS correlation
DO_APPLY_CFO_CORRECTION = 0;    %Enable CFO estimation/correction
USE_PILOT_TONES = 1;    %Enabel phase error correction
DECIMATE_RATE = INTERP_RATE;
tx_n=1;
 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the WARPLab experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    USE_AGC = false;

    NUMNODES = 1;

    %Create a vector of node objects
    global nodes;
    if(isempty(nodes))
        nodes = wl_initNodes(NUMNODES);
    end

    %Create a UDP broadcast trigger and tell each node to be ready for it
    eth_trig = wl_trigger_eth_udp_broadcast;
    wl_triggerManagerCmd(nodes,'add_ethernet_trigger',[eth_trig]);

    %Get IDs for the interfaces on the boards. Since this example assumes each
    %board has the same interface capabilities, we only need to get the IDs
    %from one of the boards
    [RFA,RFB,RFC,RFD] = wl_getInterfaceIDs(nodes(1));

    %Set up the interface for the experiment
    wl_interfaceCmd(nodes,'RF_ALL','tx_gains',3,20);
    wl_interfaceCmd(nodes,'RF_ALL','channel',2.4,11);

	wl_interfaceCmd(nodes,'RF_ALL','rx_gain_mode','manual');
	RxGainRF = 3; %Rx RF Gain in [1:3]
	RxGainBB = 10; %Rx Baseband Gain in [0:31]
	wl_interfaceCmd(nodes,'RF_ALL','rx_gains',RxGainRF,RxGainBB);

    TX_NUM_SAMPS = nodes(1).baseband.txIQLen;
    SAMP_FREQ = wl_basebandCmd(nodes(1),'tx_buff_clk_freq'); 
    node_tx = nodes(1);
    node_rx = nodes(1);
    switch tx_n
        case 1
            RF_TX=RFA;
        case 2
            RF_TX=RFB;
        case 3
            RF_TX=RFC;
    end
    RF_RX = RFB;
    
    %Set up the baseband for the experiment
    wl_basebandCmd(nodes,'tx_delay',0);
    wl_basebandCmd(nodes,'tx_length',TX_NUM_SAMPS); 
    example_mode_string = 'hw';

%Define a halfband 2x interp filter response
interp_filt2 = zeros(1,43);
interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
interp_filt2(22) = 16384;
interp_filt2 = interp_filt2./max(abs(interp_filt2));

%% Define the preamble
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);

%LTS for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

%Use 30 copies of the 16-sample STS for extra AGC settling margin
preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];

%Sanity check inputs
if(INTERP_RATE*((N_OFDM_SYMS * (N_SC + CP_LEN)) + length(preamble)) > TX_NUM_SAMPS)
    fprintf('Too many OFDM symbols for TX_NUM_SAMPS!\n');
    return;
end

%% Generate a payload
tx_data = randi(MOD_ORDER, 1, N_DATA_SYMS) - 1;

%Functions for data -> complex symbol mapping (avoids comm toolbox requirement for qammod)
modvec_bpsk =  (1/sqrt(2))  .* [-1 1];
modvec_16qam = (1/sqrt(10)) .* [-3 -1 +3 +1];

mod_fcn_bpsk = @(x) complex(modvec_bpsk(1+x),0);
mod_fcn_qpsk = @(x) complex(modvec_bpsk(1+bitshift(x, -1)), modvec_bpsk(1+mod(x, 2)));
mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));

%Map the data values on to complex symbols
switch MOD_ORDER
    case 2 %BPSK
        tx_syms = arrayfun(mod_fcn_bpsk, tx_data);
    case 4 %QPSK
        tx_syms = arrayfun(mod_fcn_qpsk, tx_data);
    case 16 %16-QAM
        tx_syms = arrayfun(mod_fcn_16qam, tx_data);      
    otherwise
        fprintf('Invalid MOD_ORDER (%d)!\n', MOD_ORDER);
        return;
end

%Reshape the symbol vector to a matrix with one column per OFDM symbol
tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYMS);

%Define the pilot tones
if(USE_PILOT_TONES)
    pilots = [1 1 -1 1].';
else
    pilots = [0 0 0 0].';
end

%Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);

%% IFFT

%Construct the IFFT input matrix
ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);

%Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(SC_IND_DATA, :) = tx_syms_mat;
ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;

%Perform the IFFT
tx_payload_mat = ifft(ifft_in_mat, N_SC, 1);

%Insert the cyclic prefix
if(CP_LEN > 0)
    tx_cp = tx_payload_mat((end-CP_LEN+1 : end), :);
    tx_payload_mat = [tx_cp; tx_payload_mat];
end

%Reshape to a vector
tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));

%Construct the full time-domain OFDM waveform
tx_vec = [preamble tx_payload_vec];

%Pad with zeros for transmission
tx_vec_padded = [tx_vec zeros(1,(TX_NUM_SAMPS/INTERP_RATE)-length(tx_vec))];

%% Interpolate
if(INTERP_RATE == 1)
    tx_vec_air = tx_vec_padded;
elseif(INTERP_RATE == 2)
    tx_vec_2x = zeros(1, 2*numel(tx_vec_padded));
    tx_vec_2x(1:2:end) = tx_vec_padded;
    tx_vec_air = filter(interp_filt2, 1, tx_vec_2x);
end

%Scale the Tx vector
tx_vec_air = TX_SCALE .* tx_vec_air ./ max(abs(tx_vec_air));
phase_tmp_log=zeros(100000,1);
%% WARPLab Tx/Rx  
    %Write the Tx waveform to the Tx node
for ii=1:100000
    wl_basebandCmd(node_tx,[RF_TX], 'write_IQ', tx_vec_air(:));

    %Enable the Tx and Rx radios
    wl_interfaceCmd(node_tx,RF_TX,'tx_en');
    wl_interfaceCmd(node_rx,RF_RX,'rx_en');

    %Enable the Tx and Rx buffers
    wl_basebandCmd(node_tx,RF_TX,'tx_buff_en');
    wl_basebandCmd(node_rx,RF_RX,'rx_buff_en');

    %Trigger the Tx/Rx cycle at both nodes
    eth_trig.send();

    %Retrieve the received waveform from the Rx node
    rx_vec_air = wl_basebandCmd(node_rx,[RF_RX],'read_IQ', 0, TX_NUM_SAMPS);
    rx_vec_air = rx_vec_air(:).';
   
    %Disable the Tx/Rx radios and buffers
    wl_basebandCmd(nodes,'RF_ALL','tx_rx_buff_dis');
    wl_interfaceCmd(nodes,'RF_ALL','tx_rx_dis');

	%% Decimate
	if(DECIMATE_RATE == 1)
		raw_rx_dec = rx_vec_air;
	elseif(DECIMATE_RATE == 2)  
		raw_rx_dec = filter(interp_filt2, 1, rx_vec_air);
		raw_rx_dec = raw_rx_dec(1:2:end);
	end

	%% Correlate for LTS

	%Complex cross correlation of Rx waveform with time-domain LTS 
	lts_corr = abs(conv(conj(fliplr(lts_t)), sign(raw_rx_dec)));

	%Skip early and late samples
	lts_corr = lts_corr(32:end-32);

	%Find all correlation peaks
	lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));

	%Select best candidate correlation peak as LTS-payload boundary
	[LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
	[lts_second_peak_index,y] = find(LTS2-LTS1 == length(lts_t));

	%Punt if no valid correlation peak was found
	if(isempty(lts_second_peak_index))
		fprintf('No LTS Correlation Peaks Found!\n');
		return;
	end

	%Set the sample indices of the payload symbols and preamble
	payload_ind = lts_peaks(max(lts_second_peak_index))+32;
	lts_ind = payload_ind-160;

		rx_cfo_est_lts = 0;

	%Apply CFO correction to raw Rx waveform
	rx_cfo_corr_t = exp(1i*2*pi*rx_cfo_est_lts*[0:length(raw_rx_dec)-1]);
	rx_dec_cfo_corr = raw_rx_dec .* rx_cfo_corr_t;

	%Re-extract LTS for channel estimate
	rx_lts = rx_dec_cfo_corr(lts_ind : lts_ind+159);
	rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
	rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);

	rx_lts1_f = fft(rx_lts1);
	rx_lts2_f = fft(rx_lts2);

	%Calculate channel estimate
	rx_H_est = lts_f .* (rx_lts1_f + rx_lts2_f)/2;
	
	phase_tmp_log(ii)=angle((rx_H_est(1)+rx_H_est(64))/2);
	
	figure(100);
	plot(1:ii,phase_tmp_log(1:ii));
	
	%% Rx payload processsing

	%Extract the payload samples (integral number of OFDM symbols following preamble)
	payload_vec = rx_dec_cfo_corr(payload_ind : payload_ind+N_OFDM_SYMS*(N_SC+CP_LEN)-1);
	payload_mat = reshape(payload_vec, (N_SC+CP_LEN), N_OFDM_SYMS);

	%Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
	payload_mat_noCP = payload_mat(CP_LEN-FFT_OFFSET+[1:N_SC], :);

	%Take the FFT
	syms_f_mat = fft(payload_mat_noCP, N_SC, 1);

	%Equalize (zero-forcing, just divide by compled chan estimates)
	syms_eq_mat = syms_f_mat ./ repmat(rx_H_est.', 1, N_OFDM_SYMS);

	%Extract the pilots and calculate per-symbol phase error
	pilots_f_mat = syms_eq_mat(SC_IND_PILOTS, :);
	pilot_phase_err = angle(mean(pilots_f_mat.*pilots_mat));
	pilot_phase_corr = repmat(exp(-1i*pilot_phase_err), N_SC, 1);

	%Apply the pilot phase correction per symbol
	syms_eq_pc_mat = syms_eq_mat .* pilot_phase_corr;
	payload_syms_mat = syms_eq_pc_mat(SC_IND_DATA, :);

	%% Demod
	rx_syms = reshape(payload_syms_mat, 1, N_DATA_SYMS);

	demod_fcn_bpsk = @(x) double(real(x)>0);
	demod_fcn_qpsk = @(x) double(2*(real(x)>0) + 1*(imag(x)>0));
	demod_fcn_16qam = @(x) (8*(real(x)>0)) + (4*(abs(real(x))<0.6325)) + (2*(imag(x)>0)) + (1*(abs(imag(x))<0.6325));

	switch(MOD_ORDER)
		case 2 %BPSK
			rx_data = arrayfun(demod_fcn_bpsk, rx_syms);
		case 4 %QPSK
			rx_data = arrayfun(demod_fcn_qpsk, rx_syms);
		case 16 %16-QAM
			rx_data = arrayfun(demod_fcn_16qam, rx_syms);    
	end

	%% Calculate Rx stats

	sym_errs = sum(tx_data ~= rx_data);
	bit_errs = length(find(dec2bin(bitxor(tx_data, rx_data),8) == '1'));
	rx_evm = sqrt(sum((real(rx_syms) - real(tx_syms)).^2 + (imag(rx_syms) - imag(tx_syms)).^2)/(length(SC_IND_DATA) * N_OFDM_SYMS));


	%% Plot Results
	cf = 0;
	%Tx sig
	cf = cf + 1;
	figure(cf); clf;

	subplot(2,1,1);
	plot(real(tx_vec_air), 'b');
	axis([0 length(tx_vec_air) -TX_SCALE TX_SCALE])
	grid on;
	title('Tx Waveform (I)');

	subplot(2,1,2);
	plot(imag(tx_vec_air), 'r');
	axis([0 length(tx_vec_air) -TX_SCALE TX_SCALE])
	grid on;
	title('Tx Waveform (Q)');

	if(WRITE_PNG_FILES)
		print(gcf,sprintf('wl_ofdm_plots_%s_txIQ',example_mode_string),'-dpng','-r96','-painters')
	end

	%Rx sig
	cf = cf + 1;
	figure(cf); clf;
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

	if(WRITE_PNG_FILES)
		print(gcf,sprintf('wl_ofdm_plots_%s_rxIQ',example_mode_string),'-dpng','-r96','-painters')
	end

	%Rx LTS corr
	cf = cf + 1;
	figure(cf); clf;
	lts_to_plot = lts_corr(1:1000);
	plot(lts_to_plot, '.-b', 'LineWidth', 1);
	hold on;
	grid on;
	line([1 length(lts_to_plot)], LTS_CORR_THRESH*max(lts_to_plot)*[1 1], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 2);
	title('LTS Correlation and Threshold')
	xlabel('Sample Index')

	cf = cf + 1;
	figure(cf); clf;

	plot(payload_syms_mat(:),'r.');
	axis square; axis(1.5*[-1 1 -1 1]);
	grid on;
	hold on;

	plot(tx_syms_mat(:),'bo');
	title('Tx and Rx Constellations')
	legend('Rx','Tx')
	
    pause(0.001);
end