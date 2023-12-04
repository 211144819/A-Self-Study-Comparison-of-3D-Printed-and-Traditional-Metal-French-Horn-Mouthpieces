% Specify the folder containing your audio files
folder = 'D:\MMusSamples\SustainNormalised\Nylon12CF'; % Change this to your folder path

% Search for all .wav files in the specified folder
audioFiles = dir(fullfile(folder, '*.wav'));
    
AudioFileMaxNumbers = size(audioFiles,1);

% Initialize the figure and set properties
figure('Name', 'Audio File Viewer', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 50, 1720, 900]);

% Define the number of columns and calculate the number of rows
numCols = 14;
numFiles = numel(audioFiles);
numRows = ceil(numFiles / numCols);

% Define the width and height of each subplot
subplotWidth = 1 / numCols; % Adjust as needed
subplotHeight = 0.96 / numRows; % Adjust as needed

% Create axes for each audio file and plot
for i = 1:AudioFileMaxNumbers
    % Read the audio file
    filePath = fullfile(folder, audioFiles(i).name);
    [audio, sampleRate] = audioread(filePath);
    
    % Calculate the time vector for the audio signal
    time = (0:(length(audio) - 1)) / sampleRate;
    
    % Calculate position for the current audio file subplot
    col = mod(i - 1, numCols);
    row = numRows - floor((i - 1) / numCols);
    xPosition = col * subplotWidth;
    yPosition = 0.96 - (row * subplotHeight ) - subplotHeight;
    
    % Create axes for the current audio file using subplot
    subplot('Position', [xPosition, yPosition, subplotWidth, subplotHeight]);
    
    % Plot the audio signal
    plot(time, audio);

    
     % Remove axis labels
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    
    % Put a small number above each graph with the audio file number
    text(0.1, 0.06, num2str(i), 'HorizontalAlignment', 'center', 'FontSize', 9);
    
    % Adjust axes limits for better visualization (optional)
    xlim([0, max(time)]);
    ylim([-0.07, 0.07]);
end