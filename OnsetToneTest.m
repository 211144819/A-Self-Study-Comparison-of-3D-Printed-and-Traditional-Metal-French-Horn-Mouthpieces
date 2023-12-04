
% List of material types, mouthpiece notes, and fundamental frequencies
folders = {'Metal', 'ABS', 'PLA', 'Nylon12CF', 'PA2200', 'Resin', 'PETG'};
notes = {'FHorn1', 'FHorn2', 'BbHorn1', 'BbHorn2'};

% Directory for audio data and output
audioDir = 'E:\ProAudio\MMusSamples\OnsetNormalised'; 
outputDir = 'E:\ProAudio\MMusSamples\OnsetNormalised';

cd(audioDir);
fileCounts = zeros(length(notes),length(folders));
dataMean = cell(length(notes),length(folders));
dataSD = cell(length(notes),length(folders));
ZScores = cell(length(notes),length(folders)-1);

window_length = 4800;
overlap = 4752;

[X, Y] = getSampleDimensions(audioDir, folders, window_length, overlap);
PSum = zeros(X,Y);

%----------------------------------------------------------------------
%Find Max value for scaling
%----------------------------------------------------------------------
maxDataHolder = zeros(length(notes),length(folders));
maxValueHolder = zeros(length(notes),1);
maxValue = 0;

%----------------------------------------------------------------------
%Calculate and Average the Spectrograms
%----------------------------------------------------------------------
for i = 1:length(notes)

    for j = 1:length(folders)
    
    PSum = zeros(X,Y);
        
    cd(audioDir);
    files = dir(fullfile(audioDir, folders{j},...
        sprintf('%s-%s-*.wav', folders{j}, notes{i})));
    fileCounts(i, j) = length(files);

        for k = 1:fileCounts(i, j)

            tonePath = fullfile(audioDir, folders{j}, files(k).name);
            [sampledata,samplerate]=audioread(tonePath);

            [S,F,T] = spectrogram(sampledata,hanning(window_length),...
                overlap,24000,samplerate,'power');
            P = (abs(S).^2);

            PSum = PSum + P;
        end
        
        PAverage = PSum/fileCounts(i,j);
        
        for m = 1:size(PAverage,1)
            for n = 1:size(PAverage,2)
                dataMean{i,j}(m,n) = PAverage(m,n);
            end
        end
        
        maxDataHolder(i,j) = max(PAverage,[],'all');

    end

end

maxValue = max(maxDataHolder,[],'all');
ratio = 1/maxValue;


%----------------------------------------------------------------------
%Calculate with Scaling, the Mean and Standard Deviation
%----------------------------------------------------------------------
for i = 1:length(notes)
    
    disp(i);
    
    for j = 1:length(folders)
   
        disp(j);
    
        for m = 1:size(dataMean{i,j},1)
            
            for n = 1:size(dataMean{i,j},2)
        
                dataMean{i,j}(m,n) = dataMean{i,j}(m,n)*ratio;
            end
                 
        end
        
        sumSquaredDeviations = zeros(X,Y); % Initialize the sum of squared deviations
        cd(audioDir);
        files = dir(fullfile(audioDir, folders{j}, sprintf('%s-%s-*.wav', folders{j}, notes{i})));
        fileCounts(i, j) = length(files);
        
        for k = 1:fileCounts(i, j)
            
            tonePath = fullfile(audioDir, folders{j}, files(k).name);
            [sampledata,samplerate]=audioread(tonePath);
            [S,F,T] = spectrogram(sampledata,hanning(window_length),overlap,24000,samplerate,'power');
            P = (abs(S).^2)*ratio;
            
            for m = 1:size(dataMean{i,j},1)
                for n = 1:size(dataMean{i,j},2)
                    
                    deviation = P(m,n) - dataMean{i,j}(m,n);
                    sumSquaredDeviations(m,n) = sumSquaredDeviations(m,n) + deviation^2;
                end
            end
        end
        
            for m = 1:size(dataMean{i,j},1)
                for n = 1:size(dataMean{i,j},2)
                    dataSD{i,j}(m,n) = sqrt(sumSquaredDeviations(m,n) / (fileCounts(i, j) - 1));
                end
            end
       graphtype(i, j, X, Y, dataMean, folders, notes,  T, F, outputDir); 
    end 
end

%------------------------------------------------------------------
%  Z-TEST
%-----------------------------------------------------------------
for i = 1:length(notes)
    disp(i);
   for j = 2:length(folders) 
       disp(j);
       for m = 1:X
           for n = 1:Y
               meanDiff = dataMean{i,1}(m,n) - dataMean{i,j}(m,n);
               SDsquared1Div = dataSD{i,1}(m,n)^2/fileCounts(i, 1);
               SDsquared2Div = dataSD{i,j}(m,n)^2/fileCounts(i, j);
               RootSum = sqrt(SDsquared1Div + SDsquared2Div);

               ZScores{i,j-1}(m,n) = meanDiff/RootSum;
    
           end
       end
       graphtype2(i, j-1, X, Y, ZScores, folders, notes,  T, F, outputDir); 
   end
end




function [X, Y] = getSampleDimensions(directory, folders, window_length, overlap)
    cd(cell2mat(fullfile(directory, folders(1))));
    filesources = dir("*.wav");
    tonePath = fullfile(directory, folders{1}, filesources(1).name);
    [sampledata,samplerate]=audioread(tonePath);

    [S] = spectrogram(sampledata,hanning(window_length),overlap,24000,samplerate,'power');
    S = abs(S);
    
    X = size(S,1);
    Y = size(S,2);
end

function [] = graphtype(i,  j, X, Y, dataMean, folders, notes,  T, F, outputDir)

for m = 1:X
    for n = 1:Y

        ZP(m,n) = dataMean{i,j}(m,n);
        dbP(m,n) = 10 * log10(ZP(m,n));
      
    end
end

            imagesc( T, F, dbP); 
            set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
            xlim([0.15,0.45]);
            set(gca,'XTick',[0.15:0.05:0.45]);
            set(gca,'XTickLabel',[0:0.05:0.3]);
            set(gca, 'XMinorTick', 'on');
            ylim([0,10000]);
            xlabel('Time');
            ylabel('Frequency (Hz)');
            cb = colorbar;
            caxis([-100 0]);
            zlim([-100,  0]);
            ylabel(cb, 'Power/frequency (dB/Hz)');
            colormap hot;
            set(gcf, 'Position', [100, 50, 600, 450]);
            
            cd(outputDir);
            saveStr = strcat(cell2mat(folders(j)),"-", cell2mat(notes(i)), ".jpg");
            figHandle = gca;
            saveas(figHandle, saveStr);

            cla(figHandle);
            
end

function [] = graphtype2(i, j, X, Y, ZScores, folders, notes,  T, F, outputDir)

for m = 1:X
    for n = 1:Y

        ZP(m,n) = ZScores{i,j}(m,n);
      
    end
end

            imagesc( T, F, ZP); 
            set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
            xlim([0.15,0.45]);
            set(gca,'XTick',[0.15:0.05:0.45]);
            set(gca,'XTickLabel',[0:0.05:0.3]);
            set(gca, 'XMinorTick', 'on');
            ylim([0,10000]);
            xlabel('Time');
            ylabel('Frequency (Hz)');
            cb = colorbar;
            caxis([-9 9]);
            zlim([-9, 9]);
            ylabel(cb, 'Z-Score');
            colormap parula;
            set(gcf, 'Position', [100, 50, 1100, 400]);
            
            cd(outputDir);
            saveStr = strcat('ZScore-', cell2mat(folders(j+1)),"-", cell2mat(notes(i)), ".jpg");
            figHandle = gca;
            saveas(figHandle, saveStr);

            cla(figHandle);
            
end