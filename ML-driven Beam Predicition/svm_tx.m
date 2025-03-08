% Load dataset
in_data = readtable('Data_General_LOS_NLOS_100itr_180_360.csv');

% Plot histogram of Beam_tx occurrences
beam_tx_values = unique(in_data.Beam_tx);
count = histcounts(categorical(in_data.Beam_tx));

% Shuffle dataset
in_data = in_data(randperm(height(in_data)), :);
disp('Data loaded and shuffled successfully.');

% Extract input features and labels
X = table2array(in_data(:, 2:end-6)); % Adjust columns based on your data
y = categorical(in_data.Beam_tx);     % Convert labels to categorical for classification

% Split dataset into training and testing sets
cv = cvpartition(size(X,1), 'HoldOut', 0.2);
X_train = X(training(cv), :);
y_train = y(training(cv), :);
X_test = X(test(cv), :);
y_test = y(test(cv), :);
disp('Dataset split into training and testing sets.');

% Train a multi-class SVM classifier using error-correcting output codes
% Optionally, specify the SVM kernel (e.g., 'linear' or 'rbf') using templateSVM:
template = templateSVM('KernelFunction', 'linear');
svm_model = fitcecoc(X_train, y_train, 'Learners', template);
disp('SVM model training complete.');

% Evaluate performance on the test set
predictions = predict(svm_model, X_test);
accuracy = mean(predictions == y_test);
fprintf('Classification Accuracy on Test Set: %.2f%%\n', accuracy * 100);

% Save the trained model for later use
save('svm_model_tx.mat', 'svm_model');

% To load the model in another file or app:
% load('svm_model.mat', 'svm_model');
% new_predictions = predict(svm_model, X_new);
