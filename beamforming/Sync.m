nodes = wl_initNodes(1);
eth_trig = wl_trigger_eth_udp_broadcast;
nodes.wl_triggerManagerCmd('add_ethernet_trigger', [eth_trig]);

% Read Trigger IDs into workspace
trig_in_ids  = wl_getTriggerInputIDs(nodes(1));
trig_out_ids = wl_getTriggerOutputIDs(nodes(1));

% For the transmit node, we will allow Ethernet to trigger the buffer baseband, the AGC, and debug0 
% (which is mapped to pin 8 on the debug header)
node_tx.wl_triggerManagerCmd('output_config_input_selection', [trig_in_ids.BASEBAND wl_getTriggerInputIDs(nodes(1)).EXT_OUT_P0],wl_getTriggerInputIDs(nodes(1)).ETH_A);
