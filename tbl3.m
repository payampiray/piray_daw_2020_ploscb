function tbl3

M1 = fit_A3;
m1 = M1(2:end,[3 4]);
m1 = cell2mat(m1);

M2 = fit_B3;
m2 = M2(2:end,[3 4]);
m2 = cell2mat(m2);

m  = [m1 m2];
m  = round(m*100)/100;
num2clip(m);

end
