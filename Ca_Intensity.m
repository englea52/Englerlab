clear all
clc

mkdir('Z:\Aditya\ET & RW\11','Results');
pathout = ['Z:\Aditya\ET & RW\11' '\' 'Results'];

acq = 25; %Acquisition rate of camera in frames per second

for k = 1:10
    cd('Z:\Aditya\ET & RW\11');
    myfilename = sprintf('Stream2-%d.txt', k);
    mydata{k} = importdata(myfilename);
    
    ROI = cell2mat(mydata(k)); %Import Region of Interest intensity .txt file
    bkgrnd = importdata('Stream2-bg.txt'); %Import background intensity .txt file
    
    [minROI, minROIidx] = min(ROI); %Minimum intensity value and index of ROI
    F0 = minROI(:,2); %Background intensity of ROI (inherent cell intensity)
    B0 = bkgrnd(minROIidx(:,2),2); %Background intensity at same index as ROI (inherent background intensity)
    
    Fratio = (ROI(:,2) - bkgrnd(:,2))/(F0 - B0); %F/F0 Ratio
    time = ROI(:,1) / (acq); %Convert frame number into time (seconds)
    
    %interpolate a linear line through minimum values
    [minpeaks, minlocs] = findpeaks(-Fratio);
    int = interp1(minlocs/acq, -minpeaks, time, 'linear');
    ind1=find(~isnan(int),1,'first');
    ind2=find(~isnan(int),1,'last');
    int(1:ind1-1)=int(ind1);
    int(ind2+1:end)=int(ind2);
    
    Fratiofinal = Fratio-int+1; %ensure baseline is 1
    
    [peaks,locs,width, prominence] = findpeaks(Fratiofinal,'MinPeakProminence',.25*(max(Fratiofinal)-min(Fratiofinal)), 'WidthReference', 'halfprom');  %Locate peaks in F ratio
    TPM = 60*length(peaks)/(time(end));  %Transients per minute
    
    peak_area = trapz(time, Fratiofinal)/length(peaks); %Mean Peak Area
    mean_amp = mean(peaks);
    
    irreg_amp = std(prominence)/mean(prominence); %Irregularity of amplitude
    
    if length(locs) > 2
        %Difference in sec between each Ca peak
        d = diff(locs)/acq;
        timing_irreg = std(d)/mean(d); %Irregularity of peak timing
    else
        timing_irreg = nan;
    end
    
    Frat(:,k) = Fratiofinal;
    TPMfinal(k,:) = TPM;
    peak_areafinal(k,:) = peak_area;
    mean_ampfinal(k,:) = mean_amp;
    irreg_ampfinal(k,:) = irreg_amp;
    timing_irregfinal(k,:) = timing_irreg;
    
    a = [TPMfinal,peak_areafinal,mean_ampfinal,irreg_ampfinal,timing_irregfinal];
    
    cd(pathout)
    %export time and F ratio to excel to plot in prism
    filename = 'Stream2 F ratio.xlsx';
    y = [time, Frat];
    xlswrite(filename, y);
    
    col_header={'Transients per min','Mean Peak Area','Mean Peak Amplitude','Amplitude Irreg','Timing Irreg'}; %Row cell array (for column labels)
    xlswrite('Stream2 Analysis.xls',a,'Sheet1','A2'); %Write data
    xlswrite('Stream2 Analysis.xls',col_header,'Sheet1','A1'); %Write column header
end