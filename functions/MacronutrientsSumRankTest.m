function [] = MacronutrientsSumRankTest(macronutrientsMatrix,idxk)
%MACRONUTRIENTSSUMRANKTEST performs the Wilcoxon rank sum test between
%cluster 1 and cluster 2 for: Fat, Protein, Carbohydrates, Sum of Fat and
%Protein
% ATTENTION: designed for a 2-cluster solution.

arguments (Input)
    macronutrientsMatrix % #observations x 4 matrix [Fat Prot Carb Sum_Fat_Prot]
    idxk % cluster assignment
end

cluster1_set = macronutrientsMatrix(idxk == 1,:);
cluster2_set = macronutrientsMatrix(idxk == 2,:);

%% WILCOXON RANK SUM TEST
[p_set_fat, h_set_fat] = ranksum(cluster1_set(:,1),cluster2_set(:,1));
[p_set_prot, h_set_prot] = ranksum(cluster1_set(:,2),cluster2_set(:,2));

disp(' --- WILCOXON RANK SUM TEST ---')

disp(' ')
disp(' The p-values for cluster 1 vs cluster 2 are')
disp(' |    FAT    |    PROT    | ')
disp([' ' num2str(p_set_fat) '    ' num2str(p_set_prot)])
if h_set_fat == 1 && h_set_prot == 1
    disp('Null Hypothesis rejected')
end

disp(' ')
disp('-------------------------------')

end
