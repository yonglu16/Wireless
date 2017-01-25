% Implementing Wision
num_ant=2;
antenna_num=0:num_ant-1;
d=0.06; %distances between antennas
lambda=0.12; %full wave length


figure(4);
angles=0:pi/12:2*pi-pi/12;
% Base function used in Wision
basis=exp(-j* repmat([0:num_ant-1],24,1).*repmat(cos(angles'),1,num_ant) * 2*pi*d / lambda);
abs(sum((repmat([csi_processed(1).' csi_processed(2).'],24,1).*basis)'));
plot(rad2deg(angles),katabi);
pause(1);
end
