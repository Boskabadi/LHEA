% % Extract time and data from the timeseries object
% timeData = DI.Time;    % Extract the time vector
% valueData = DI.Data;   % Extract the data values
% 
% % Create a table for better organization
% dataTable = table(timeData, valueData, 'VariableNames', {'Time', 'Data'});
% 
% % Save the table to an Excel file
% writetable(dataTable, 'D1.xlsx');

% Extract time and data from the timeseries object DI
timeData = DIQ.Time;       % Original time points
valueData = DIQ.Data;      % Original data values

% Define the integer time intervals
integerTimes = ceil(min(timeData)):floor(max(timeData)); 

% Interpolate data at integer time intervals
newData = interp1(timeData, valueData, integerTimes, 'linear');

% Create a table for exporting
outputTable = table(integerTimes', newData, 'VariableNames', {'Time', 'Data'});

% Save the table to an Excel file
writetable(outputTable, 'DIQ.xlsx');

% Display a success message
disp('Time series data interpolated to integer time intervals and saved');
