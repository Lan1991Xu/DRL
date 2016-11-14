function cons_set = hw_update_constraint(idx,newcons,cons_set)
% output: constraint_set(1*n) with fields:
% idx; cons;
if isempty(cons_set)
  cons_set(1).idx  = idx;
  cons_set(1).cons = newcons;
else
  idxs = [cons_set(:).idx];
  tmp  = find(idxs == idx);
  if isempty(tmp)
      num_cons = length(cons_set);
      cons_set(num_cons+1).idx  = idx;
      cons_set(num_cons+1).cons = newcons;
  else
      cons_set(tmp).cons.strong_id = union(cons_set(tmp).cons.strong_id,newcons.strong_id);
      cons_set(tmp).cons.weak_id = union(cons_set(tmp).cons.weak_id,newcons.weak_id);
      cons_set(tmp).cons.match_id = union(cons_set(tmp).cons.match_id,newcons.match_id);
  end    
end