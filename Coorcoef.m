clear all
cd('C:\Users\Aditya\Desktop\10 vs 30 kPa Calcium\1-9\txt files\30 kPa\Round 2\Results');
d = dir('*.xlsx');
data = [];
for k = 1:length(d)
    data = [data; importdata(d(k).name)];   
    data(k).Sheet1(:,1) = [];
end

for j = 1:length(d)
    b = corrcoef(data(j).Sheet1);
    blower = tril(b, -1);
    bupper = triu(b,1);
    result = blower(:, 1:end-1) + bupper(:, 2:end);
    coorelation(j) = mean2(result);
end

xlswrite('Coorelation Coefficient', transpose(coorelation));