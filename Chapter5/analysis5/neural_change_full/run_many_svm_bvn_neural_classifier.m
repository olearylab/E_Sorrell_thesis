function [accuracies_cell] = run_many_svm_bvn_neural_classifier(z_cell,virmen_cell,tbt_cell,nbins,kept_types,balance_types,keep_incorrect,l_r_sep)
% 02/08/2023
%%


num_mice = size(virmen_cell,1);
num_days = size(virmen_cell,2);
% vars_all = cell(num_mice,num_days);

accuracies_cell = cell(num_mice,num_days);
% stds_all = cell(num_mice,num_days);

for m = 1:num_mice
    for d = 1:num_days
        if ~isempty(virmen_cell{m,d})
            % [ss,centres] = prep_data_for_bvn_lstm(virmen_cell{m,d},tbt_cell{m,d},nbins,keep_incorrect,kept_types,balance_types);
           
            [ss,centres] = prep_data_for_neural_bvn_lstm(z_cell{m,d},virmen_cell{m,d},tbt_cell{m,d},nbins,keep_incorrect,kept_types,balance_types);
            
            
            % Additonal option to run separately on left and right trials
            if l_r_sep
                s = ss;
                kept_l = ismember(ss.origTrialType,[1,3]);
                s.binnedNeuralData2 = ss.binnedNeuralData2(kept_l,:,:);
                s.trialType = ss.trialType(kept_l);
                s.correctVec = ss.correctVec(kept_l);
                
                [accuracies_l] = svm_neural_classifier(s);
                
                s = ss;
                kept_r = ismember(ss.origTrialType,[2,4]);
                s.binnedNeuralData2 = ss.binnedNeuralData2(kept_r,:,:);
                s.trialType = ss.trialType(kept_r);
                s.correctVec = ss.correctVec(kept_r);
                [accuracies_r] = svm_neural_classifier(s);
                % average of left and right accuracies.
                accuracies_mat = (accuracies_l + accuracies_r)/2;
            else

                [accuracies_mat] = svm_neural_classifier(s);
                
            end

            accuracies_cell{m,d} = accuracies_mat;
        end
    end
end