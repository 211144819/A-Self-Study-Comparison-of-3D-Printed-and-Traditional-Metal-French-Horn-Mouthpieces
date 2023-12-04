
directory = 'E:\ProAudio\MMusSamples\RAW\Session1';
savedir = 'E:\ProAudio\MMusSamples\SeperatedTones';
cd(directory);
filesources = dir("*.wav");



for f = 1:size(filesources,1)
    cd(directory);
    
    NoteCount = 0;
    loadname = filesources(f).name;
    
    for session = 1:3
        
        SessionFolder = strcat("E:\ProAudio\MMusSamples\RAW\Session", num2str(session));
        SessionFolderChar = char(SessionFolder);
        cd(SessionFolderChar);

        
        fullfilename = fullfile(SessionFolderChar, loadname);
        [path,filename,ext]=fileparts(fullfilename);

        [data,samplerate]=audioread(loadname);

        samplelength = length(data);
        sumy = zeros(samplelength,1);
        avey = zeros(samplelength,1);

        volumeThresh = 0.02;
        volumeThreshAvey = 0.005;
        avSamples = samplerate/400;
        NoteOn = 0;
        SilenceTime = 1;
        SilenceSamples = samplerate*SilenceTime;

        %identify number of tones

        %create smoother averages to check for note entry and exit
        for i = 1:(samplelength-avSamples)

            for j = 1:avSamples

            sumy(i) = sumy(i) + abs(data(i+j));

            end
            avey(i) = sumy(i)/avSamples;

        end

        % Define parameters for rise and fall detection
        riseThreshold = volumeThresh;  % Adjust as needed
        fallThreshold = volumeThresh;  % Adjust as needed
        samplesForRise = 24000;  % Number of samples for rise detection
        samplesForFall = 24000;  % Number of samples for fall detection

        NoteOn = 0;
        NoteOnL = 0;
        interval = 20;

        % Identify start and end points
        for k = (samplesForRise + 1):(interval):(samplelength - samplesForFall)
            % Calculate the average value for the current and previous samples
            avgBefore = (sum(avey(k - samplesForRise:k - 1))/size(avey(k - samplesForRise:k - 1),1));
            avgCurrent = (sum(avey(k:k + samplesForFall -1))/size(avey(k:k + samplesForFall -1),1));

            % Detect note onset (rise)
            if (avgCurrent - avgBefore) >= riseThreshold
                NoteOn = 1;
                NoteOnL = k;
            end

            % Detect note offset (fall)
            if (((avgBefore - avgCurrent) >= fallThreshold) && (avey(k) < volumeThreshAvey))
                if (NoteOn == 1)
                    NoteCount = NoteCount + 1;
                    % Append start and end locations to arrays
                    NoteOn = 0;
                    
                    %save single sample to disk
                    cd(savedir);
                    filesaved = 0;
                    
                    for z = (NoteOnL-24000):(NoteOnL)
                        
                        
                        if (avey(z) >= volumeThreshAvey)
                    
                            for sample = (z - SilenceSamples):(k + SilenceSamples)

                                noteData((sample - z + SilenceSamples + 1),1) = data(sample);

                            end

                            saveStr = strcat(filename,"-", num2str(NoteCount , '%03d'), ".wav");
                            saveStrChar = char(saveStr);

                            audiowrite(saveStrChar,noteData,samplerate);

                            %reset noteData - not doing this causes tails from
                            %longer older samples
                            clear noteData;
                            filesaved = 1;
                        end
                       
                        if filesaved == 1
                            break;
                            
                        end
                    end
                    
                end
            end
        end

    end

end


