% Find index for each demographic group
function demog_mat_indices = find_demog_rows2(state_matrix, StateMatCols, demog_group_def, DemogTblCols)

    demog_mat_indices = ...
    find((state_matrix(:,StateMatCols.age) >= demog_group_def(DemogTblCols.age_min,1)) &...
        (state_matrix(:,StateMatCols.age) <= demog_group_def(DemogTblCols.age_max,1)) &...
         state_matrix(:,StateMatCols.race) == demog_group_def(DemogTblCols.race,1));
end