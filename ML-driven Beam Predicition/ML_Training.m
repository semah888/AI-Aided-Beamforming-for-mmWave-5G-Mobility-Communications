% Load dataset
in_data = readtable('Data_General_LOS_NLOS_100itr_180_360.csv');

% Plot histogram of Beam_tx occurrences
beam_tx_values = unique(in_data.Beam_tx);
count = histcounts(categorical(in_data.Beam_tx));
bar(beam_tx_values, count);
xlabel('Beam Index at Tx');
ylabel('Number of Samples');
title('Distribution of Beam Index at Tx');

% Shuffle dataset
in_data = in_data(randperm(height(in_data)), :);

disp('Data loaded and shuffled successfully.');

% Extract input features and labels
X = table2array(in_data(:, 2:end-6)); % Assuming last column is the label
y = categorical(in_data.Beam_tx); % Convert labels to categorical for classification

% Convert categorical labels to one-hot encoding
%y = onehotencode(y, 2);

% Split dataset into training and testing sets
cv = cvpartition(size(X,1), 'HoldOut', 0.2);
X_train = X(training(cv), :);
y_train = y(training(cv), :);
X_test = X(test(cv), :);
y_test = y(test(cv), :);

disp('Dataset split into training and testing sets.');

% Define neural network architecture
numClasses = numel(unique(y_train));
layers = [ 
    featureInputLayer(size(X,2))
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(32)
    reluLayer
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

% Set training options
options = trainingOptions('adam', ...
    'MaxEpochs', 3, ...
    'MiniBatchSize', 128, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {X_test, y_test}, ...
    'Verbose', true, ...
    'Plots', 'training-progress');

% Train the network
dnn = trainNetwork(X_train, y_train, layers, options);

disp('Model training complete.');

% Evaluate performance
predictions = classify(dnn, X_test);
accuracy = mean(predictions == categorical(in_data.Beam_tx(test(cv))));
fprintf('Classification Accuracy on Test Set: %.2f%%\n', accuracy * 100);
save('dnnModel_tx.mat', 'dnn');