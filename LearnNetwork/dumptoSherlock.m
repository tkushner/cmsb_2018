
BaseName='DallaMan_pt1';

load(strcat(BaseName, '.mat'))
no_of_inputs = 37 + 37;
hidden_layer_1_neurons = 8;
hidden_layer_2_neurons = 8;
no_of_outputs = 1;
%%
% Start with the glucose values, followed by the insulin values

input_layer_weights = zeros(hidden_layer_1_neurons, no_of_inputs);
input_layer_biases = zeros(hidden_layer_1_neurons);

for i=1:hidden_layer_1_neurons
    if(i<= (hidden_layer_1_neurons/2))
        for j=1:no_of_inputs
            if(j <= no_of_inputs/2)
                input_layer_weights(i,j) = W_gluc_1(i,j);
            else
                input_layer_weights(i,j) = 0.0;
            end
        end
        input_layer_biases(i) = b_gluc_1(i);
    else
        for j=1:no_of_inputs
            if(j <= no_of_inputs/2 )
                input_layer_weights(i,j) = 0.0;
            else
                input_layer_weights(i,j) = W_insulin_1(i-(hidden_layer_1_neurons/2), j - no_of_inputs/2);
            end            
        end
        input_layer_biases(i) = b_insulin_1(i-hidden_layer_1_neurons/2);
    end
end
%%
fileID = fopen(strcat(BaseName,'_APNN.nt'),'w');
data = no_of_inputs; % No_of_inputs
fprintf(fileID,'%d \n',data);
data = 1; % No_of_outputs
fprintf(fileID,'%d \n',data);
data = 2; % No_of_hidden_layers
fprintf(fileID,'%d \n',data);
data = hidden_layer_1_neurons;
fprintf(fileID,'%d \n',data);
data = hidden_layer_2_neurons;
fprintf(fileID,'%d \n',data);

%  Input Layer data : 
for i=1:hidden_layer_1_neurons
    for j = 1:no_of_inputs
        data = input_layer_weights(i,j);
        fprintf(fileID,'%f \n',data);
    end
    data = input_layer_biases(i);
    fprintf(fileID,'%f \n',data);
end

%  Middle Layer data : 
for i=1:hidden_layer_2_neurons
    for j = 1:hidden_layer_1_neurons
        data = W_joint_1(i,j);
        fprintf(fileID,'%f \n',data);
    end
    data = b_joint_1(i);
    fprintf(fileID,'%f \n',data);
end


%  Output Layer data : 
for i=1:no_of_outputs
    for j = 1:hidden_layer_2_neurons
        data = W_joint_2(i,j);
        fprintf(fileID,'%f \n',data);
    end
    data = b_joint_2(i);
    fprintf(fileID,'%f \n',data);
end

fclose(fileID);


