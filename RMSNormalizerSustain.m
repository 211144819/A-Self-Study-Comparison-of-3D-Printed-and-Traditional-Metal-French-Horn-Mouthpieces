directory = 'E:\ProAudio\MMusSamples\SeperatedTones';
savedir = 'E:\ProAudio\MMusSamples\SustainNormalised';
cd(directory);
filesources = dir(directory);
dirFlags = [filesources.isdir];
subFolders = filesources(dirFlags);
subFolderNames = {subFolders(3:end).name}; %obtain folder list without the extras


for DirNum = 1:size(subFolderNames,2)
    
    cd(savedir);
    mkdir(cell2mat(subFolderNames(DirNum)));

    MPDirectory = strcat(directory, '\', cell2mat(subFolderNames(DirNum)));
    cd(MPDirectory);
    filelist = dir("*.wav");
    
    AudioFileMaxNumbers = size(filelist,1);
    DesiredDecibals = -18;

    DesiredPower = db2pow(DesiredDecibals);

    PassFile = 0;
    PassedFilesNum = 0;


    for i = 1:AudioFileMaxNumbers
        filename = filelist(i).name;
        FullFilePath = fullfile(MPDirectory, filename);
        [Data,samplerate]=audioread(FullFilePath);
        Data = Data((2*samplerate):(3.5*samplerate-1));

        RootMeanSquare = rms(Data);
        Ratio = DesiredPower/RootMeanSquare;

        Datamax = max(Data);
        Datamin = min(Data);

        %check for clipping and filter out bad files
        if (Datamax >= 0.99) || (Datamin <= -0.99)
            disp('Sample has clipping')
            PassFile = 0;
        elseif (Ratio*Datamax > 1) || (Ratio*Datamin < -1)
            disp('Sample will clip if normalised')
            PassFile = 0;
        else
            PassedFilesNum = PassedFilesNum +1;
            PassFile = 1;
        end


        %Normalise File
        if (PassFile == 1)

            NewData = Ratio*Data;
            
            
            
            SaveDirectory = strcat(savedir, '\', cell2mat(subFolderNames(DirNum)));
            cd(SaveDirectory);
    
            savename = filename(1:(size(filename,2)-4));
            saveStr = strcat(savename, "-SustainNormalised.wav");
            saveStrChar = char(saveStr);

            audiowrite(saveStrChar,NewData,samplerate);
        end

    end
    
end