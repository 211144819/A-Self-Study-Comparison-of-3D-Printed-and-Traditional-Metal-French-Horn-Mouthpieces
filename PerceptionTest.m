

mainFunction();

function mainFunction()
    
    silenceBetweenTones = 0;
    numTests = 200; 

    folders = {'Metal', 'ABS', 'PLA', 'Nylon12CF', 'PA2200', 'Resin', 'PETG'};
    notes = {'FHorn1', 'FHorn2', 'BbHorn1', 'BbHorn2'};

    audioDir = 'D:\MMusSamples\FullNormalised'; 
    outputDir = 'D:\MMusSamples\FullNormalised'; 


    
    figure('Name', 'Sound Perception Test', 'NumberTitle', 'off');
    axis off;

   
    buttonA = uicontrol('Style', 'pushbutton', 'String', 'A', 'Position', [100, 50, 50, 30]);
    buttonB = uicontrol('Style', 'pushbutton', 'String', 'B', 'Position', [200, 50, 50, 30]);
    buttonC = uicontrol('Style', 'pushbutton', 'String', 'C', 'Position', [300, 50, 50, 30]);

    
    buttonPause = uicontrol('Style', 'pushbutton', 'String', 'Pause', 'Position', [400, 50, 50, 30]);

   
    buttonResume = uicontrol('Style', 'pushbutton', 'String', 'Resume', 'Position', ...
        [400, 50, 100, 50], 'Visible', 'off');

    
    sharedData = struct('userChoice', '', 'playing', false);

    
    set(buttonA, 'Callback', @(src, event) buttonCallback('A'));
    set(buttonB, 'Callback', @(src, event) buttonCallback('B'));
    set(buttonC, 'Callback', @(src, event) buttonCallback('C'));
    set(buttonPause, 'Callback', @(src, event) ...
        pauseCallback(buttonA, buttonB, buttonC, buttonPause, buttonResume));
    set(buttonResume, 'Callback', @(src, event) ...
        resumeCallback(buttonA, buttonB, buttonC, buttonPause, buttonResume));

    totalTests = zeros(length(folders), length(notes));
    correctAnswers = zeros(length(folders), length(notes));
    fileCounts = zeros(length(folders), length(notes));

    for i = 1:length(folders)
        for j = 1:length(notes)
            files = dir(fullfile(audioDir, folders{i}, sprintf('%s-%s-*.wav', folders{i}, notes{j})));
            fileCounts(i, j) = length(files);
        end
    end

    for testIndex = 1:numTests
        noteIndex = randi([1, length(notes)]);
        folderIndex = randi([2, length(folders)]); 
        order = randperm(3);
        metalTones = randperm(fileCounts(1, noteIndex), 2);
        otherTone = randperm(fileCounts(folderIndex, noteIndex), 1);

        
        sharedData.playing = true;
        for i = order
            
            if i <= 2
                toneName = fullfile(audioDir, folders{1}, ...
                    sprintf('%s-%s-%03d-FullNormalised.wav', folders{1}, notes{noteIndex}, metalTones(i)));
            else
                toneName = fullfile(audioDir, folders{folderIndex}, ...
                    sprintf('%s-%s-%03d-FullNormalised.wav', folders{folderIndex}, notes{noteIndex}, otherTone));
            end

            [tone, fs] = audioread(toneName);
            sound(tone, fs);
            sharedData.playing = true; 

            
            pause((length(tone) / fs) + silenceBetweenTones);
            sharedData.playing = false; 
        end

        while isempty(sharedData.userChoice)
            pause(0.1); 
        end

        while sharedData.playing
            pause(0.1); 
        end

        choice = sharedData.userChoice;
        sharedData.userChoice = ''; 

        if choice == char(order(3) + 64) 
            correctAnswers(folderIndex - 1, noteIndex) = correctAnswers(folderIndex - 1, noteIndex) + 1;
        end

        totalTests(folderIndex - 1, noteIndex) = totalTests(folderIndex - 1, noteIndex) + 1;
    end
    fid=fopen(fullfile(outputDir,'results2.txt'),'w');
    for i=2:length(folders)
       for j=1:length(notes)
           fprintf(fid,'Folder: %s Note: %s Tests: %d Correct: %d FileCount: %d\n',...
               folders{i}, notes{j}, totalTests(i-1,j), correctAnswers(i-1,j), fileCounts(i,j));
       end
    end
    fclose(fid);

    close;

    
     function buttonCallback(choice)
            sharedData.userChoice = choice;
        end
    

    
    function pauseCallback(buttonA, buttonB, buttonC, buttonPause, buttonResume)
        sharedData = guidata(gcf);
        set(buttonA, 'Visible', 'off');
        set(buttonB, 'Visible', 'off');
        set(buttonC, 'Visible', 'off');
        set(buttonPause, 'Visible', 'off');
        set(buttonResume, 'Visible', 'on');
        guidata(gcf, sharedData);
    end

    
    function resumeCallback(buttonA, buttonB, buttonC, buttonPause, buttonResume)
        sharedData = guidata(gcf);
        set(buttonA, 'Visible', 'on');
        set(buttonB, 'Visible', 'on');
        set(buttonC, 'Visible', 'on');
        set(buttonPause, 'Visible', 'on');
        set(buttonResume, 'Visible', 'off');
        guidata(gcf, sharedData);
    end

end
