directory = 'D:\MMusSamples\Scales\Trimmed';
savedir = 'D:\MMusSamples\Scales';
cd(directory);
filesources = dir("*.wav");

Algorithms = {'NCF', 'PEF', 'CEP', 'LHS', 'SRH'};




for j = 1:length(filesources)
    
    cd(directory);
    
    
    [data,samplerate]=audioread(filesources(j).name);
    
    fullfilename = fullfile(directory, filesources(j).name);
    [path,filename,ext]=fileparts(fullfilename);
    
    
    
    for i = 1:length(Algorithms)
 
         cd(directory);

        [f0,loc] = pitch(data,samplerate,'Range',[100,2000],'WindowLength', 2400, 'OverlapLength', 2280, 'Method',cell2mat(Algorithms(i)));
        t1 = (0:length(f0)-1)/(length(loc)/(length(data)/samplerate));
        plot(t1,f0) 
        xlabel("Time (s)")
        ylabel("Frequency (Hz)")
        grid minor
        axis tight
        ylim([150,500]);
        
        set(gcf, 'Position', [100, 50, 600, 150]);
        

        
        
        cd(savedir);
        fileName = extractBefore(filesources(j).name, ".");
        saveStr = strcat(fileName,"-", Algorithms(i), ".jpg");
        titleStr = strcat(fileName,"-", Algorithms(i));
        title(titleStr);
        figHandle = gca;
        saveas(figHandle, saveStr);

        cla(figHandle);
        clear f0 loc t1;
        
        
    end 
    
    clear data;
        
end