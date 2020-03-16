function tbl1

[tx_sim, mx] = sim_rec;
mx = [tx_sim; mx];

m  = round(mx'*100)/100;
% % num2clip(m);

end
