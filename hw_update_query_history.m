function query_history = hw_update_query_history(query_history,accumu_idx,current_idx,pair_rank,cons)

if length(query_history) < accumu_idx
    % initialize a new query, 0-click
    query_history{accumu_idx}.current_idx = current_idx;
    query_history{accumu_idx}.pair_rank{1} = pair_rank';
    query_history{accumu_idx}.cons{1} = []; % no constraints
    
elseif length(query_history) == accumu_idx
    % update on a old query, n-click
    query_history{accumu_idx}.pair_rank{end+1} = pair_rank';
    query_history{accumu_idx}.cons{end+1} = cons;
    
end