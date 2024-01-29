function [new_virmen_data] = rearrange_virmen(virmen_data)

% Dan's new virmen data has 20 rows, with the new ball voltages slotted in
% at rows 12-14.

% I will rearrange so that the new voltages are at 13-15, and the old
% voltages are the final 3 (18-20), although I could probably remove
% altogether

new_virmen_data = zeros(size(virmen_data));
new_virmen_data(1:11,:) = virmen_data(1:11,:);
new_virmen_data(12,:) = virmen_data(15,:);
new_virmen_data(13:15,:) = virmen_data(12:14,:);
new_virmen_data(16:17,:) = virmen_data(19:20,:);
new_virmen_data(18:20,:) = virmen_data(16:18,:);