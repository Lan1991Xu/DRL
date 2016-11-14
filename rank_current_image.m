function result = rank_current_image(all_features,all_camid,current_idx,curr_metric,incremental_method)
% give the current ranking result, based on the current metric
% show the top 14 ranking result
current_camid = all_camid(current_idx);
current_feat  = all_features(current_idx,:);
othercam_mask = find(all_camid ~= current_camid);
othercam_feat = all_features(othercam_mask,:);
%dists = pdist2(current_feat,othercam_feat,'mahalanobis',curr_metric);
%dif_ftr = repmat(current_feat,[size(othercam_feat,1),1]) - othercam_feat;
switch incremental_method
    case {'LEGO','our'}
        dif_ftr = bsxfun(@minus,current_feat,othercam_feat);
        foo = dif_ftr*curr_metric;
        dists = sum(foo.*dif_ftr,2);
        [~,rank_list] = sort(dists,'ascend');
        result        = othercam_mask(rank_list);
    case {'POP','Chunxiao'}
        dif_ftr = bsxfun(@minus,current_feat,othercam_feat);
        score_l2 = normalise_01(-sum(dif_ftr.^2,2));
        if ~isempty(curr_metric.svs)
            options = make_options('gamma_I',0.1,'gamma_A',0.1,'NN',20,'KernelParam',2.75, ...
            'Hinge', 1, 'UseBias', 0);
            K = calckernel(options, dif_ftr, dif_ftr); 
            score_lap = K(:, curr_metric.svs)*curr_metric.alpha+curr_metric.b;
            score_lap = normalise_01(score_lap);
        else
            score_lap = zeros(size(score_l2));
        end
        score_final = score_l2 * 0.2 + score_lap * 0.8;
        [~,rank_list] = sort(score_final,'descend');
        result        = othercam_mask(rank_list);
     case {'PRF'}
        dif_ftr = bsxfun(@minus,current_feat,othercam_feat);
        score_l2 = normalise_01(-sum(dif_ftr.^2,2));
        if ~isempty(curr_metric)
            pred_label_fake = ones(size(dif_ftr, 1), 1);
            [~, ~, pred_prob] = svmpredict(pred_label_fake, dif_ftr, curr_metric);
            score_re = (pred_prob(:,1));
            score_re = normalise_01(score_re);
        else
            score_re = zeros(size(score_l2));
        end
        score_final = score_l2 * 0.2 + score_re * 0.6;
        [~,rank_list] = sort(score_final,'descend');
        result        = othercam_mask(rank_list);  
     case {'NPRF'}
        dif_ftr = bsxfun(@minus,current_feat,othercam_feat);
        score_l2 = normalise_01(-sum(dif_ftr.^2,2));
        if ~isempty(curr_metric)
            pred_label_fake = ones(size(dif_ftr, 1), 1);
            [~, ~, pred_prob] = svmpredict(pred_label_fake, dif_ftr, curr_metric);
            score_re = (pred_prob(:));
            score_re = normalise_01(score_re);
        else
            score_re = zeros(size(score_l2));
        end
        score_final = score_l2 * 0.2 + score_re * 0.9;
        [~,rank_list] = sort(score_final,'descend');
        result        = othercam_mask(rank_list);     
      case {'EMR','Rocchio'}
        dif_ftr = bsxfun(@minus,current_feat,othercam_feat);
        score_l2 = normalise_01(-sum(dif_ftr.^2,2));
        if ~isempty(curr_metric)
            score_re = normalise_01(curr_metric);
        else
            score_re = zeros(size(score_l2));
        end
        score_final = score_l2 * 0.2 + score_re * 0.8;
        [~,rank_list] = sort(score_final,'descend');
        result        = othercam_mask(rank_list);   
    otherwise
        error('Unknown incremental-updating scheme.')
end

function scoren = normalise_01(score)
scoren = (score - min(score) + eps) / (max(score)-min(score)+eps); 