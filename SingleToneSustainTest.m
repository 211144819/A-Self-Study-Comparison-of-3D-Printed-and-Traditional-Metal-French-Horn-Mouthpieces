% List of material types, mouthpiece notes, and fundamental frequencies
folders = {'Metal', 'ABS', 'PLA', 'Nylon12CF', 'PA2200', 'Resin', 'PETG'};
notes = {'FHorn1', 'FHorn2', 'BbHorn1', 'BbHorn2'};
colours = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[0 1 1],[1 0 1],[0 0 0]};

% Directory for audio data and output
audioDir = 'D:\MMusSamples\SustainNormalised'; 
outputDir = 'D:\MMusSamples\SustainNormalised';

fileCounts = zeros(length(notes),length(folders));
data = cell(length(notes),length(folders));
maxDataHolder = zeros(length(notes),length(folders));
maxValueHolder = zeros(length(notes),1);
maxValue = 0;


samplerate = 48000;
window_length = 24000;
frequencies = (0:window_length-1) * samplerate / window_length;
  

for i = 1:length(notes)
    
    cla;
    
    for j = 1:length(folders)
        
        files = dir(fullfile(audioDir, folders{j}, sprintf('%s-%s-*.wav', folders{j}, notes{i})));
        fileCounts(i, j) = length(files);
        data{i, j} = struct('Amplitude', {}, 'Frequency', {}, 'ScaledAmplitude', {});
        statdata{i, j} = struct('Means', {}, 'Frequencies', {}, 'StandardDeviations', {});
        
        DataSum = zeros(window_length,1);

        for k = 1:fileCounts(i,j)
        
            tonePath = fullfile(audioDir, folders{j}, sprintf('%s-%s-%03d-SustainNormalised.wav', folders{j}, notes{i}, k));
            cd(audioDir);
            [sampledata,samplerate]=audioread(tonePath);
            
            
            hanning_window = hanning(window_length);
            
            useddata = sampledata(1:window_length);

            % Apply the window to the data
            windowed_sample = useddata(1:window_length) .* hanning_window;

            DFT = fft(windowed_sample);
            Power = abs(DFT).^2;

            %scaling factors to compensate for windowing.
            Scaling_factor = 1.63;
            Scaled_Power = Scaling_factor * Power;
            
            DataSum = DataSum + Scaled_Power;
 
        end
        DataSum = DataSum/fileCounts(i,j);
        
        
        for k = 1:length(DataSum)
            data{i, j}(k).Amplitude = DataSum(k);
            data{i, j}(k).Frequency = frequencies(k);
        end
        maxDataHolder(i,j) = max(DataSum);
        
    end
    
    maxValueHolder(i,1) = max(maxDataHolder,[],'all');
    
end

maxValue = max(maxValueHolder,[],'all');
ratio = 1/maxValue;

HarmonicCount = zeros(length(notes),length(folders));

for i = 1:length(notes)
    
    
    for j = 1:length(folders)

    cla;
    
        for k = 1:length(data{i,j})
            
        data{i,j}(k).ScaledAmplitude = data{i, j}(k).Amplitude * ratio;
       
        end
        
        for k = 1:window_length
            dataforPlot(k) = data{i,j}(k).ScaledAmplitude;
        end

        dataforPlot = dataforPlot';
        dBPower = 10 * log10(dataforPlot);

        fundamentals = {174.6, 261.6, 233.1, 349.2};
        funFreq = fundamentals{i};
        startSample = floor((funFreq/2-20));

        %exclude frequencies below fundamental frequency
        trimPower = dBPower(startSample:length(dBPower)/2);
        trimfrequencies = frequencies(startSample:length(frequencies)/2);
        
        [Amplitude, HarmonicFrequency, NumberofHarmonics] = HarmonicCalc(trimPower , trimfrequencies, funFreq);

        
        HarmonicCount(i,j) = NumberofHarmonics;
        
        for k = 1:length(Amplitude)
            statdata{i, j}(k).Means = 10^(Amplitude(k)/10);
            statdata{i, j}(k).Frequencies = HarmonicFrequency(k);
        end
        
    
        hold on;
        plot(trimfrequencies,trimPower, 'color', cell2mat(colours(3)));
        title('test')

        xlim([20,7000]);    %use [20,3000] for 110Hz fundamental, [20,5000] for 220 Hz, and [20,8000] for 440Hz.
        ylim([-100,0]);  
        ylabel('Power (dB)');
        xlabel('Frequency (Hz)');
        set(gca,'XTick',[0:500:7000]);
        set(gca,'XTickLabel', [0:500:7000]);
        set(gca,'YTick',[-100:10:0]);
        set(gca,'YTickLabel',[-100:10:0]);
        set(gcf, 'Position', [100, 50, 1000, 600]);
        
        cd(outputDir);
        saveStr = strcat(cell2mat(folders(j)),"-", cell2mat(notes(i)), ".jpg");
        figHandle = gca;
        saveas(figHandle, saveStr);

        cla(figHandle);
        
    end
    
end

maxHarmonics = min(HarmonicCount,[],'all');

%calculate the standard deviations
% Initialize array for standard deviations
SDData = zeros(length(notes),length(folders),maxHarmonics);

% Calculate standard deviations
for i = 1:length(notes)
    for j = 1:length(folders)      
        for m = 1:maxHarmonics
            sumSquaredDeviations = 0; % Initialize the sum of squared deviations
            for k = 1:fileCounts(i, j)
                
                tonePath = fullfile(audioDir, folders{j}, sprintf('%s-%s-%03d-SustainNormalised.wav', folders{j}, notes{i}, k));
                [sampledata,samplerate]=audioread(tonePath);
                hanning_window = hanning(window_length);
                useddata = sampledata(1:window_length);
                % Apply the window to the data
                windowed_sample = useddata(1:window_length) .* hanning_window;
                DFT = fft(windowed_sample);
                Power = abs(DFT).^2;
                %scaling factors to compensate for windowing.
                Scaling_factor = 1.63;
                Scaled_Power = Scaling_factor * Power;
                Power = Scaled_Power * ratio;
                
                harmonicPeak = peakfinder(Power, statdata{i, j}(m).Frequencies/2);
                deviation = Power(harmonicPeak) - statdata{i, j}(m).Means;
                sumSquaredDeviations = sumSquaredDeviations + deviation^2;
            end
            % Calculate the standard deviation for the current combination
            SDData(i, j, m) = sqrt(sumSquaredDeviations / (fileCounts(i, j) - 1));
            statdata{i, j}(m).StandardDeviations = SDData(i, j, m);
        end
    end
end

%save excels
for i = 1:length(notes)
    % Define the Excel file name for the current material
    excelFileName = fullfile(outputDir, [notes{i}, '_statdata.xlsx']);

    % Create cell array for the data
    dataExcel = cell(2 + 2 * length(folders), 2 + maxHarmonics);

    % Set harmonics labels in the data array
    dataExcel(1, 3:end) = labels;

    % Loop through folders
    for j = 1:length(folders)
        % Set folder label in the data array
        dataExcel{2 + (j - 1) * 2, 1} = folders{j};
        dataExcel{3 + (j - 1) * 2, 1} = folders{j};
        
        % Set "Mean" and "SD" labels
        dataExcel{2 + (j - 1) * 2, 2} = 'Mean';
        dataExcel{3 + (j - 1) * 2, 2} = 'SD';

        % Loop through harmonics
        for m = 1:maxHarmonics
            % Set Mean and SD values in the data array
            dataExcel{2 + (j - 1) * 2, 2 + m} = statdata{i, j}(m).Means;
            dataExcel{3 + (j - 1) * 2, 2 + m} = statdata{i, j}(m).StandardDeviations;
        end
    end

    % Write the data to an Excel file
    writecell(dataExcel, excelFileName);
end

%Z-test of 2 populations - pop 1 is always metal
for i = 1:length(notes)
   for j = 2:length(folders) 
       for m = 1:maxHarmonics
           meanDiff = statdata{i, 1}(m).Means - statdata{i, j}(m).Means;
           SDsquared1Div = (statdata{i, 1}(m).StandardDeviations)^2/fileCounts(i, j);
           SDsquared2Div = (statdata{i, j}(m).StandardDeviations)^2/fileCounts(i, j);
           RootSum = sqrt(SDsquared1Div + SDsquared2Div);
           
           ZScores(i,j-1,m) = meanDiff/RootSum;
    
 
       end
   end
end

%SAVE Z_SCORES TO EXCEL
materials = {'ABS', 'PLA', 'Nylon12CF', 'PA2200', 'Resin', 'PETG'};

% Create an Excel filename for the Z-scores
zScoreExcelFileName = fullfile(outputDir, 'ZScores.xlsx');

% Create a cell array for the Z-scores data
zScoreData = cell(length(notes) * length(materials) + 1, maxHarmonics + 3);

% Create cell array for the harmonics labels (e.g., 1, 2, 3, ...)
harmonicsLabels = cell(1, maxHarmonics);
for m = 1:maxHarmonics
    harmonicsLabels{m} = num2str(m);
end

% Assign column labels to the first row
colLabels = [{'Material', 'Note'}, harmonicsLabels];
zScoreData(1, 1:length(colLabels)) = colLabels;

% Populate the Z-scores data cell array
rowIdx = 2;
for j = 1:length(notes)
    note = notes{j};
    for i = 1:length(materials)
        material = materials{i};
        
        % Add Material and Note values
        zScoreData{rowIdx, 1} = material;
        zScoreData{rowIdx, 2} = note;
        
        % Add Z-scores for each harmonic in separate columns
        for m = 1:maxHarmonics
            zScoreData{rowIdx, m + 2} = ZScores(j, i, m);
        end
        
        rowIdx = rowIdx + 1;
    end
end

% Write the Z-scores data to an Excel file
writecell(zScoreData, zScoreExcelFileName);


function [peakIndex] = peakfinder(Power, index);

        Index = (round(index*0.92):round(index*1.08));
        maxValues = zeros(length(Index),1);
        for a = Index(1):Index(length(Index))
            maxValues(a - Index(1) + 1) = Power(a);
        end
        Max = max(maxValues, [],'all')
        MaxIndex  = find(maxValues == Max);
        peakIndex = index + (MaxIndex - find(Index == index));
end
        
        
function [Amplitude, HarmonicFrequency, NumberofHarmonics] = HarmonicCalc(trimPower , trimfrequencies, funFreq)

    MaxHarmonics = 0;
    bestp = 0;
    bestPeaks = zeros(100, 1);
    bestLocs = zeros(1, 100);
    for p = 5:50
        
        [peaks,locs] = findpeaks(trimPower,trimfrequencies, 'MinPeakProminence',p);
        goodHarmonics = 0;
        goodPeaks = zeros(length(locs), 1);
        goodLocs = zeros(1, length(locs));
        
        for f = 1:length(locs)
                if (locs(f) >= (goodHarmonics + 1) * funFreq * 0.98 && locs(f) <= (goodHarmonics + 1) * funFreq * 1.02)
                    goodHarmonics = goodHarmonics + 1;
                    goodPeaks(goodHarmonics) = peaks(f);
                    goodLocs(goodHarmonics) = locs(f);
                    
                    if locs(f) >= (goodHarmonics + 1) * funFreq * 1.05
                        break;
                    end
                    if f > MaxHarmonics
                        bestp = p;
                    end
                end              
        end 
        
        if goodHarmonics > MaxHarmonics
            MaxHarmonics = goodHarmonics - 1;
            bestPeaks = goodPeaks;
            bestLocs = goodLocs;
        end
        
    end
    
    
    Amplitude = bestPeaks(1:MaxHarmonics);
    HarmonicFrequency = bestLocs(1:MaxHarmonics);
    NumberofHarmonics = MaxHarmonics;

end
        
        
        
