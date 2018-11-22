clear all
cd('C:\Users\Aditya\Desktop\10 vs 30 kPa IF\Striation Organization\1-9\30 kPa');
for k=1:5;
    myfilename = sprintf('FITC10-%d.txt', k);
    mydata{k} = importdata(myfilename);
    y_1 = cell2mat(mydata(k)); %Import Region of Interest intensity .txt file
    y = y_1.data;
    %y = mydata{1,k}.data(:,2);
    z = fftshift(fft(y));
    rawpowerspec=abs(z).^2;
    powerspecnorm=(rawpowerspec)./(max(max(rawpowerspec)));
    pks = findpeaks(powerspecnorm);
    b = sort(pks);
    c = length(b);
    d(k) = b(c-1);
    %d = powerspecnorm(1:floor(n/2));
end
org = transpose(d)