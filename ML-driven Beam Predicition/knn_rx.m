% Load dataset
in_data = readtable('Data_General_LOS_NLOS_100itr_180_360.csv');

% Plot histogram of Beam_tx occurrences
beam_rx_values = unique(in_data.Beam_rx);
count = histcounts(categorical(in_data.Beam_rx));

% Shuffle dataset
in_data = in_data(randperm(height(in_data)), :);
disp('Data loaded and shuffled successfully.');

% Extract input features and labels
X = table2array(in_data(:, 2:end-6)); % Adjust based on your data
y = categorical(in_data.Beam_rx);     % Convert labels to categorical for classification

% Split dataset into training and testing sets
cv = cvpartition(size(X,1), 'HoldOut', 0.2);
X_train = X(training(cv), :);
y_train = y(training(cv), :);
X_test = X(test(cv), :);
y_test = y(test(cv), :);
disp('Dataset split into training and testing sets.');

% Train a k-NN classifier
% You can adjust 'NumNeighbors' as needed (e.g., 5 is common)
knn_model = fitcknn(X_train, y_train, 'NumNeighbors', 3);
disp('k-NN model training complete.');

% Evaluate performance on the test set
predictions = predict(knn_model, X_test);
accuracy = mean(predictions == y_test);
fprintf('Classification Accuracy on Test Set: %.2f%%\n', accuracy * 100);

% If you need to use the model in another file or app, save it:
save('knn_model_rx.mat', 'knn_model');

% Later, you can load the model without retraining:
% load('knn_model.mat', 'knn_model');
% And apply it to new data:
% new_predictions = predict(knn_model, X_new);
