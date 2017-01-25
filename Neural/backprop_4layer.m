% backprop a per-epoch backpropagation training for a multilayer feedforward
%          neural network.
%   Network = bbackprop(Layers,N,M,SatisfactoryMSE,Input,Desired) 
%   returns Network, a structure having the fields:
%     structure - vector of number of nodes at each layer
%     weights - a cell array specifiying the final weight matrices computed
%     epochs - the epochs required for training
%     mse - the mean squared error at termination
%   weights is a cell array specifying the final weight matrices computed by 
%   minimizing the mean squared error between the Desired output and the 
%   actual output of the network given a set of training samples: Input and 
%   the SatisfactoryMSE (satisfactory mean squared error) 
%
%   Input:
%    Layers - a vector of integers specifying the number of nodes at each
%     layer, i.e for all i, Layers(i) = number of nodes at layer i, there
%     must be at least three layers and the input layer Layers(1) must
%     equal the dimension of each vector in Input, likewise, Layers(end) 
%     must be equal to the dimension of each vector in Desired
%     N - training rate for network learning [0.1 - 0.9]
%     M - momentum for the weight update rule [0.0 - 0.9]
%     SatisfactoryMSE - the mse at which to terminate computation
%     Input - the training samples, a P-by-N matrix, where each Input[p] is
%      a training vector
%     Desired - the desired outputs, a P-by-M matrix where each Desired[p]
%      is the desired output for the corresponding input Input[p]
%
%   This algorithm uses the hyperbolic tangent node function 
%   2/(1+e^(-net)) - 1, suitable for use with bipolar data
%   
%   NOTE: due to its generality this algorithm is not as efficient as a 
%   one designed for a specific problem. If the number of desired layers is 
%   known ahead of time, it is better to a) remove the use of cells for the
%   weight matrices and explicitly create the number of desired weight matrices
%   and 'unfold' the loops inside the presentation loop. That is, explicitly 
%   calculate the input and output of each each layer one by one and subsequently 
%   the modified error and weight matrix modifications b) remove calculation of 
%   momentum if it is not required
%
% Author: Dale Patterson
% $Version: 3.1.2 $ $Date: 2.25.06 $
% 
% function Network = bbackprop(L,n,m,smse,X,D)
%%%%% VERIFICATION PHASE %%%%%
% determine number of input samples, desired output and their dimensions
%each row is a sample
% X=fft(eye(10,10));
% X=(eye(10,10));
% D=[-0.9:0.2:1].';
% load('nlos_txloc.mat')
% ifft_log(1:10,:)=[];
% loc_log(:,1:10)=[];
% for ii=1:size(ifft_log,1)
%    ifft_log(ii,:)=ifft_log(ii,:)./max(abs(ifft_log(ii,:)));
% end
% X_overall=(ifft_log.').';
% D_overall=((loc_log.'/5)-0.5)*2;
% train_loc=[1:70 81:90];test_loc=71:80;
X_overall=abs(useful_ifft)./repmat(max(abs(useful_ifft),[],2),1,size(useful_ifft,2));
% X_overall=[abs(useful_ifft)];
D_overall=(loc_log.'-13.5)/13.5;
train_loc=[1:10 21:30 41:50 61:70 81:90 101:110 111:120 121:130 141:150 161:170 181:190 201:210 221:230 241:250 251:260 261:270];
test_loc=[train_loc+10];test_loc(end-9:end)=[];
X=X_overall(train_loc,:);
X=repmat(X,100,1);
X=X+(randn(size(X)).*[zeros(size(X,1)/2,size(X,2));ones(size(X,1)/2,size(X,2))/105]);
D=repmat(D_overall(train_loc,:),100,1);
L=[14 14 6 1];
n=0.01;
momentum=0.9;
smse=0.00001;
err_log=[];
[P,N] = size(X);
[Pd,M] = size(D);

% make user that each input vector has a corresponding desired output
if P ~= Pd 
    error('backprop:invalidTrainingAndDesired', ...
          'The number of input vectors and desired ouput do not match');
end

% make sure that at least 3 layers have been specified and that the 
% the dimensions of the specified input layer and output layer are
% equivalent to the dimensions of the input vectors and desired output
if length(L) < 3 
    error('backprop:invalidNetworkStructure','The network must have at least 3 layers');
else
    if N ~= L(1) || M ~= L(end)
        e = sprintf('Dimensions of input (%d) does not match input layer (%d)',N,L(1));
        error('backprop:invalidLayerSize', e);
    elseif M ~= L(end)
        e = sprintf('Dimensions of output (%d) does not match output layer (%d)',M,L(end));
        error('backprop:invalidLayerSize', e);    
    end
end

%%%%% INITIALIZATION PHASE %%%%%
nLayers = length(L); % we'll use the number of layers often  

% randomize the weight matrices (uniform random values in [-1 1], there
% is a weight matrix between each layer of nodes. Each layer (exclusive the 
% output layer) has a bias node whose activation is always 1, that is, the 
% node function is C(net) = 1. Furthermore, there is a link from each node
% in layer i to the bias node in layer j (the last row of each matrix)
% because it is less computationally expensive then the alternative. The 
% weights of all links to bias nodes are irrelevant and are defined as 0
w = cell(nLayers-1,1); % a weight matrix between each layer
for i=1:nLayers-2
    if(i==1)
        w{i} = [1 - 2.*rand(L(i+1),L(i)+1) ; zeros(2,L(i)+1)];
    elseif(i==2)
        w{i} = [1 - 2.*rand(L(i+1),L(i)+2) ; zeros(1,L(i)+2)];
        w{i}(:,end)=0;
%         w{i}(:,end-1)=0; %m
    else
        w{i} = [1 - 2.*rand(L(i+1),L(i)+1) ; zeros(1,L(i)+1)];
    end
end
w{end} = 1 - 2.*rand(L(end),L(end-1)+1);

% initialize stopping conditions
mse = Inf;  % assuming the intial weight matrices are bad
epochs = 0;

%%%%% PREALLOCATION PHASE %%%%%
% for faster computation preallocate activation,net,prev_w and sum_w

% Activation: there is an activation matrix a{i} for each layer in the 
% network such that a{1} = the network input and a{end} = network output
% Since we're doing batch mode, each activation matrix a{i} is a 
% P-by-K (P=num of samples,K=nodes at layer i) matrix such that 
% a{i}(j) denotes the activation vector of layer i for the jth input and 
% a{i}(j,k) is the activation(output) of the kth node in layer i for the jth 
% input
a = cell(nLayers,1);  % one activation matrix for each layer
a{1} = [X ones(P,1)]; % a{1} is the input + '1' for the bias node activation
                      % a{1} remains the same throught the computation
for i=2:nLayers-1
    if(i==3)
        a{i} = ones(P,L(i)+2); % third layer includes m and c (P-by-Nodes+1) 
    else
        a{i} = ones(P,L(i)+1); % inner layers include a bias node (P-by-Nodes+1) 
    end
end

a{end} = ones(P,L(end));   % no bias node at output layer

% Net: like activation, there is a net matrix net{i} for each layer
% exclusive the input such that net{i} = sum(w(i,j) * a(j)) for j = i-1
% and each net matrix net{i} is a P-by-K matrix such that net{i}(j) denotes 
% the net vector at layer i for the jth sample and net{i}(j,k) denotes the
% net input at node k of the ith layer for the jth sample
net = cell(nLayers-1,1); % one net matrix for each layer exclusive input
for i=1:nLayers-2;
    if(i==2)
        net{i} = ones(P,L(i+1)+2); % affix m&c nodes
    else
        net{i} = ones(P,L(i+1)+1); % affix bias node 
    end
end
net{end} = ones(P,L(end));

% Since we're using batch mode and momentum, two additional matrices are
% needed: prev_dw is the delta weight matrices at time (t-1) and sum_dw
% is the sum of the delta weights for each presentation of the input
% the notation here is the same as net and activation, that is prev_dw{i} 
% is P-by-K matrix where prev_dw{i} is the delta weight matrix for all samples 
% at time (t-1) and sum_dw{i} is a P-by-K matrix where sum_dw{i} is the
% the sum of the weight matrix at layer i for all samples
prev_dw = cell(nLayers-1,1);
sum_dw = cell(nLayers-1,1);
for i=1:nLayers-1
    prev_dw{i} = zeros(size(w{i})); % prev_dw starts at 0
    sum_dw{i} = zeros(size(w{i}));
end    

% loop until computational bounds are exceeded or the network has converged
% to a satisfactory condition. We allow for 30000 epochs here, it may be
% necessary to increase or decrease this bound depending on the number of 
% training
while mse > smse && epochs < 20000
    % FEEDFORWARD PHASE: calculate input/output off each layer for all samples
    for i=1:nLayers-1
        if(i==2)           
           net{i} = a{i}(:,1:end-2) * (w{i}(:,1:end-2)).';
        else
           net{i} = a{i} * w{i}.'; % compute inputs to current layer
        end
        
        % compute activation(output of current layer, for all layers
        % exclusive the output, the last node is the bias node and
        % its activation is 1
        if i < nLayers-1 % inner layers
            if(i==1)
                %linear activation
                a{i+1} = [net{i}(:,1:end-2) ones(P,2)];
                a{i+1}= [a{i} ones(P,1)];
            elseif (i==2)
                %sqrt activation
                m=real(repmat(w{i}(:,end-1).',P,1));
                k=real(repmat(w{i}(:,end).',P,1));
                im1=2*k./(1+m.*m);
                tof=real(net{i});
                im2=(tof.*tof)-((4.*k.*k.*m.*m)./((1+m.*m).^2));
%                 im2(im2<0)=0;
                a{i+1} = [im1-sqrt(im2)];
%                 a{i+1} = real([net{i}]);
                a{i+1} = [a{i+1}(:,1:end-1) ones(P,1)];
            end
%             a{i+1} = [2./(1+exp(-net{i}(:,1:end-1)))-1 ones(P,1)];
        else             % output layers
            a{i+1} = 2 ./ (1 + exp(-real(net{i}))) - 1;
        end
    end
    
    % calculate sum squared error of all samples
    err = (D-a{end});       % save this for later
    sse = sum(sum(abs(err).^2)); % sum of the error for all samples, and all nodes
    
    % BACKPROPAGATION PHASE: calculate the modified error at the output layer: 
    % S'(Output) * (D-Output) in this case S'(Output) = (1+Output)*(1-Output)
    % then starting at the output layer, calculate the sum of the weight 
    % matrices for all samples: LearningRate * ModifiedError * Activation
    % then backpropagate the error such that the modified error for this
    % layer is: S'(Activation) * ModifiedError * weight matrix
    delta = err .* (1 + a{end}) .* (1 - a{end});
    for i=nLayers-1:-1:1
        sum_dw{i} = n * delta.' * a{i};
        if(i==2) 
            temp_dw_m= n * delta_m.' * a{i};
            temp_dw_k= n * delta_k.' * a{i};
            sum_dw{2}(:,end-1)=temp_dw_m(:,end-1);
            sum_dw{2}(:,end)=temp_dw_k(:,end);
        end
        if i > 1
            if(i==nLayers-1)
                prev_delta=delta;
                m=real(repmat(w{2}(:,end-1).',P,1));
                k=real(repmat(w{2}(:,end).',P,1));                
                im1=2*k./(1+m.*m);
                tof=real(net{2});
                im2=(tof.*tof)-((4.*k.*k.*m.*m)./((1+m.*m)).^2);
%                  im2(im2<0)=1e-100;
%                  delta =  (delta*w{i});
%                  delta_m=delta;
%                  delta_k=delta;
                im25=(-1./(2*(sqrt(im2))));
                delta = (im25.*tof.*2).*((prev_delta*w{3}));      
                im3=2*m./((1+m.*m).^2);
                im4=2.*(m.*m.*m./((1+m.*m).^3));
                delta_m=(-2.*k.*2.*m./(1+m.*m).^2)-4.*k.*k.*(im25.*(im3-im4)).*((prev_delta*w{3}));
                delta_k=((-2./(1+m.*m))+im25.*k.*im3.*4).*((prev_delta*w{3}));
                delta(:,end)=0;%weights between bias  node of reflector layer and prev layer neurons. does not matter, can be commented
                delta_m(:,end)=0;
                delta_k(:,end)=0;
%                 delta_m(:,:)=0;delta_k(:,:)=0;
            elseif(i==2)
                delta =  (delta*w{i});
            end
        end
    end
    
    % update the prev_w, weight matrices, epoch count and mse
    for i=1:nLayers-1
        % we have the sum of the delta weights, divide through by the 
        % number of samples and add momentum * delta weight at (t-1)
        % finally, update the weight matrices
        prev_dw{i} = (sum_dw{i} ./ P) + (momentum * prev_dw{i});
        w{i} = w{i} + prev_dw{i};
    end   
    epochs = epochs + 1;
    mse = sse/(P*M); % mse = 1/P * 1/M * summed squared error
    err_log=[err_log mse];
    if(mod(epochs,100)==0)
        figure(1)
        plot(real(a{4}(1:size(train_loc,2))))
        hold on
        plot(real(D(1:size(train_loc,2))))
        hold off
        pause(0.1)
        epochs
    end
    
    

X_test=X_overall(test_loc,:);
D_test=D_overall(test_loc,:);
[P_test,N_test] = size(X_test);
[Pd_test,M_test] = size(D_test);

a_test = cell(nLayers,1);  % one activation matrix for each layer
a_test{1} = [X_test ones(P_test,1)]; % a{1} is the input + '1' for the bias node activation
    for i=1:nLayers-1
        if(i==2)           
           net{i} = a_test{i}(:,1:end-2) * (w{i}(:,1:end-2)).';
        else
           net{i} = a_test{i} * w{i}.'; % compute inputs to current layer
        end
        
        % compute activation(output of current layer, for all layers
        % exclusive the output, the last node is the bias node and
        % its activation is 1
        if i < nLayers-1 % inner layers
            if(i==1)
                %linear activation
                a_test{i+1} = [net{i}(:,1:end-2) ones(P_test,2)];
                a_test{i+1}= [a_test{i} ones(P_test,1)];
            elseif (i==2)
                %sqrt activation
                m=real(repmat(w{i}(:,end-1).',P_test,1));
                k=real(repmat(w{i}(:,end).',P_test,1));
                im1=2*k./(1+m.*m);
                tof=real(net{i});
                im2=(tof.*tof)-((4.*k.*k.*m.*m)./((1+m.*m).^2));
%                 im2(im2<0)=0;
                a_test{i+1} = [im1-sqrt(im2)];
%                 a_test{i+1} = real([net{i}]);
                a_test{i+1} = [a_test{i+1}(:,1:end-1) ones(P_test,1)];
            end
%             a_test{i+1} = [2./(1+exp(-net{i}(:,1:end-1)))-1 ones(P_test,1)];
        else             % output layers
            a_test{i+1} = 2 ./ (1 + exp(-real(net{i}))) - 1;
        end
    end
    figure(10)
    plot(a_test{4}); hold on;
    plot(D_test)
    hold off

end

% return the trained network
Network.structure = L;
Network.weights = w;
Network.epochs = epochs;
Network.mse = mse;




X_test=X_overall(test_loc,:);
D_test=D_overall(test_loc,:);
[P_test,N_test] = size(X_test);
[Pd_test,M_test] = size(D_test);

a_test = cell(nLayers,1);  % one activation matrix for each layer
a_test{1} = [X_test ones(P_test,1)]; % a{1} is the input + '1' for the bias node activation
    for i=1:nLayers-1
        if(i==2)           
           net{i} = a_test{i}(:,1:end-2) * (w{i}(:,1:end-2)).';
        else
           net{i} = a_test{i} * w{i}.'; % compute inputs to current layer
        end
        
        % compute activation(output of current layer, for all layers
        % exclusive the output, the last node is the bias node and
        % its activation is 1
        if i < nLayers-1 % inner layers
            if(i==1)
                %linear activation
                a_test{i+1} = [net{i}(:,1:end-2) ones(P_test,2)];
            elseif (i==2)
                %sqrt activation
                m=real(repmat(w{i}(:,end-1).',P_test,1));
                k=real(repmat(w{i}(:,end).',P_test,1));
                im1=2*k./(1+m.*m);
                tof=real(net{i});
                im2=(tof.*tof)-((4.*k.*k.*m.*m)./((1+m.*m).^2));
%                 im2(im2<0)=0;
                a_test{i+1} = [im1-sqrt(im2)];
%                 a_test{i+1} = real([net{i}]);
                a_test{i+1} = [a_test{i+1}(:,1:end-1) ones(P_test,1)];
            end
%             a_test{i+1} = [2./(1+exp(-net{i}(:,1:end-1)))-1 ones(P_test,1)];
        else             % output layers
            a_test{i+1} = 2 ./ (1 + exp(-real(net{i}))) - 1;
        end
    end