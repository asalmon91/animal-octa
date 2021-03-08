function ocu_mat = subtractBackground(ocu_mat, bg)
%subtractBackground Subtracts the average background signal of the whole
%volume

if exist('bg', 'var') == 0
    bg = mean(mean(ocu_mat, 3), 2);
end
ocu_mat = ocu_mat - bg;


end