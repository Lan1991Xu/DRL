function varargout = GUI_online_active(varargin)
% GUI_ONLINE_ACTIVE MATLAB code for GUI_online_active.fig
%      GUI_ONLINE_ACTIVE, by itself, creates a new GUI_ONLINE_ACTIVE or raises the existing
%      singleton*.
%
%      H = GUI_ONLINE_ACTIVE returns the handle to a new GUI_ONLINE_ACTIVE or the handle to
%      the existing singleton*.
%
%      GUI_ONLINE_ACTIVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_ONLINE_ACTIVE.M with the given input arguments.
%
%      GUI_ONLINE_ACTIVE('Property','Value',...) creates a new GUI_ONLINE_ACTIVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_online_active_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_online_active_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_online_active

% Last Modified by GUIDE v2.5 03-Oct-2016 11:12:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_online_active_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_online_active_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI_online_active is made visible.
function GUI_online_active_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_online_active (see VARARGIN)

% Choose default command line output for GUI_online_active
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_online_active wait for user response (see UIRESUME)
% uiwait(handles.figure1);
rng('default');
rng(0);

% global variables declare
% epoch for saving temperal history
global epoch;              epoch = 25;

% setting panel
global dataset_name;       dataset_name   = 'viper';
global feature_name;       feature_name   = 'KCCA_feat';
global label_budget;       label_budget   = 100;
global used_budget;        used_budget    = 0;

% method panel
global incremental_method; incremental_method = 'LEGO';
global active_method;      active_method  = 'random';
global check_on_test;      check_on_test  = 0;
global query_viewid;       query_viewid   = 1;

% for whole training-set
global labelled_set;       labelled_set   = [];
global unlabelled_set;     unlabelled_set = [];
global num_trn;            num_trn        = [];

global all_images;         all_images     = [];
global all_features;       all_features   = [];
global all_pid;            all_pid        = [];
global all_camid;          all_camid      = [];
global all_clsid;          all_clsid      = [];
global affM;               affM           = [];
global nGroups;            nGroups        = 1;
global cls_centers;        cls_centers    = [];

% pre-defined order list
global predefined_list;    predefined_list = [];

% for testing set
global tst_set;            tst_set        = [];

% for current image
global current_idx;        current_idx    = [];
global pair_idx;           pair_idx       = [];
global label_flg;          label_flg      = false;
global curr_cluster;       curr_cluster   = [];
global accumu_idx;         accumu_idx     = 0;
global rank_history;       rank_history   = [];
global true_rank;          true_rank      = [];
global query_history;      query_history  = [];

%*** distance models
global curr_metric;        curr_metric    = [];

% constraint set
global constraint_set;     constraint_set = [];

% last-step history recover
global last_step_history;  last_step_history = [];

% paths
global data_path;          data_path      = './data/';
global model_path;         model_path     = './model/';
addpath(genpath('./code/'));

% check-boxes ids
global strong_checkid;     strong_checkid = [];
global weak_checkid;       weak_checkid   = [];
global match_checkid;      match_checkid  = [];


% check-boxes for convenience
global hndall;
hndall =       [handles.cand_1, handles.cand_2, handles.cand_3, ...
                handles.cand_4, handles.cand_5, handles.cand_6, ...
                handles.cand_7, handles.cand_8, handles.cand_9, ...
                handles.cand_10, handles.cand_11, handles.cand_12, ...
                handles.cand_13, handles.cand_14];

global strongchkall;
strongchkall = [handles.checkbox_strong1, handles.checkbox_strong2, handles.checkbox_strong3, ...
                handles.checkbox_strong4, handles.checkbox_strong5, handles.checkbox_strong6, ...
                handles.checkbox_strong7, handles.checkbox_strong8, handles.checkbox_strong9, ...
                handles.checkbox_strong10, handles.checkbox_strong11, handles.checkbox_strong12, ...
                handles.checkbox_strong13, handles.checkbox_strong14];

global weakchkall;
weakchkall =   [handles.checkbox_weak1, handles.checkbox_weak2, handles.checkbox_weak3, ...
                handles.checkbox_weak4, handles.checkbox_weak5, handles.checkbox_weak6, ...
                handles.checkbox_weak7, handles.checkbox_weak8, handles.checkbox_weak9, ...
                handles.checkbox_weak10, handles.checkbox_weak11, handles.checkbox_weak12, ...
                handles.checkbox_weak13, handles.checkbox_weak14];
            
global matchchkall;
matchchkall =   [handles.checkbox_match1, handles.checkbox_match2, handles.checkbox_match3, ...
                handles.checkbox_match4, handles.checkbox_match5, handles.checkbox_match6, ...
                handles.checkbox_match7, handles.checkbox_match8, handles.checkbox_match9, ...
                handles.checkbox_match10, handles.checkbox_match11, handles.checkbox_match12, ...
                handles.checkbox_match13, handles.checkbox_match14];   
            
% timers
global start_time;  start_time = [];
global finish_time; finish_time = [];

% show list
global show_order_from; show_order_from = 1; 

% --- Outputs from this function are returned to the command line.
function varargout = GUI_online_active_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_rank.
function pushbutton_rank_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global labelled_set;            global unlabelled_set;         
global active_method;           global current_idx;
%global pair_idx;                
global all_features;
global all_camid;               global curr_metric;
global rank_history;            global all_images;
global all_pid;                 global label_flg;
global strong_checkid;          global weak_checkid;
global match_checkid;           global show_order_from;
global curr_cluster;            global all_clsid;
global accumu_idx;              global incremental_method;
global true_rank;               global query_history;
global epoch;
true_rank = [];

tmp_pairidx = show_next_image(current_idx,handles);

% rank current image.
ranking_list = rank_current_image(all_features,all_camid,...
               current_idx,curr_metric{curr_cluster},incremental_method);
rank_history = ranking_list;

% set checkid to zeros, clear check boxes
strong_checkid = zeros(size(rank_history));
weak_checkid   = zeros(size(rank_history));
match_checkid  = zeros(size(rank_history));
show_order_from = 1; list_length = 14;
show_click_result(show_order_from,list_length);

% show rank result  
[true_rank,pair_rank] = update_truerank_history2(ranking_list, tmp_pairidx,true_rank);
show_truerank_history2(handles.axes_rank1,true_rank);
show_rank_result(all_images,all_pid,all_camid,ranking_list,handles);

% update query_history
query_history = hw_update_query_history(query_history,accumu_idx,current_idx,pair_rank,[]);


% --- Executes on button press in pushbutton_next.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global labelled_set;            global unlabelled_set;         
global active_method;           global current_idx;
%global pair_idx;                
global all_features;
global all_camid;               global curr_metric;
global rank_history;            global all_images;
global all_pid;                 global label_flg;
global strong_checkid;          global weak_checkid;
global match_checkid;           global show_order_from;
global curr_cluster;            global all_clsid;
global accumu_idx;              global incremental_method;
global true_rank;               global query_history;
global epoch;
true_rank = [];

% remove current-idx from unlabelled set, and add it to labelled set
if label_flg
    % manage label/unlabel sets.
    unlabelled_set = unlabelled_set(unlabelled_set~=current_idx);
    labelled_set   = [labelled_set, current_idx];
    label_flg = false;
    accumu_idx = accumu_idx + 1;
    
    %*** later: save ranking vs clicks
    % save result if epoch reached
    if mod(accumu_idx-1,epoch) == 0
       save_epoch(floor((accumu_idx-1)/epoch)); 
    end
    
    % set accumulation id window
    set(handles.edit_accumulate, 'string', sprintf('%4d',accumu_idx));
    
    %*** if incre==POP, need to declare an empty structure for next curr_metric
    if strcmp(incremental_method,'POP')
    curr_metric{accumu_idx}.svs = []; curr_metric{accumu_idx}.b   = [];
    curr_metric{accumu_idx}.w   = []; 
    end
    
else % remove the last one from query_history    
    query_history = query_history(1:end-1);
    
end

clear_rank_result(handles);


% first-by-first, need to push the current_metric and current_idx into
% storage. Since next step is to learn a new local metric.
% *** something needs to do later.
% after storing needs to initialize a new metric using either eye or
% convariance.

% second, clear laststep_rank_history

current_idx = ...
pick_next_image(labelled_set,unlabelled_set,active_method,accumu_idx);
tmp_pairidx = show_next_image(current_idx,handles);

% which cluster the image belongs?
curr_cluster = hw_which_cluster(all_clsid,current_idx,accumu_idx,incremental_method);



% clear accumulation labelling cost.

% --- Executes on button press in pushbutton_start.
function pushbutton_start_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load dataset, features, handling the indices
% pick a random image, do initial matching
global dataset_name;            global all_images;
global feature_name;            global all_features;
global labelled_set;            global all_pid;
global unlabelled_set;          global all_camid;
global data_path;               % global num_trn;
global active_method;           global current_idx;
%global pair_idx;                
global start_time;
global rank_history;            global curr_metric;
global used_budget;             global label_budget;
global strong_checkid;          global weak_checkid;
global match_checkid;           global constraint_set;     
constraint_set = [];            curr_metric = [];
global show_order_from;         global check_on_test;
global tst_set;                 global query_viewid;
global all_clsid;               global affM;
global nGroups;                 global curr_cluster;
global cls_centers;             global incremental_method;
global accumu_idx;              global query_history;

global true_rank; true_rank = [];

% reset time clock
if isfield(handles,'tmr')
  stop(handles.tmr);  delete(handles.tmr); handles = rmfield(handles,'tmr');
end

% load train-data
[all_images,all_features,all_pid,all_camid] = ...
hw_load_trndata(data_path,dataset_name,feature_name);

% load tst-data
if check_on_test
    tst_set = hw_load_tstdata(data_path,dataset_name,feature_name);
end

%*** load graph? Need to think when choose active mthods.
[all_clsid, affM, nGroups] = hw_load_graph(data_path,dataset_name,...
                                     feature_name,all_camid,nGroups);

% start timer. display in the clock window.
start_time = tic;
handles.tmr = hw_GUI_clock(handles.edit_timer,start_time,0);
guidata(hObject,handles);

% initialize metric; using I[dxd]. equivallant to l2-distance. Also
% calculate cluster centers
[curr_metric,cls_centers] = hw_initialize_metric(all_features,...
                 all_camid,query_viewid,all_clsid,nGroups,incremental_method);
             
% pick and show image.
% num_trn = size(all_features,1);
labelled_set = [];
unlabelled_set = find(all_camid == query_viewid);
accumu_idx = 1;
current_idx = ...
pick_next_image(labelled_set,unlabelled_set,active_method,accumu_idx);
tmp_pairidx = show_next_image(current_idx,handles);


% set accumulation id window
set(handles.edit_accumulate, 'string', sprintf('%4d',accumu_idx));

% which cluster the image belongs?
curr_cluster = hw_which_cluster(all_clsid,current_idx,accumu_idx,incremental_method);

% rank current image.
ranking_list = rank_current_image(all_features,all_camid,...
               current_idx,curr_metric{curr_cluster},incremental_method);
rank_history = ranking_list;

% set checkid to zeros, clear checkboxes
strong_checkid = zeros(size(rank_history));
weak_checkid   = zeros(size(rank_history));
match_checkid  = zeros(size(rank_history));
show_order_from = 1; list_length = 14;
show_click_result(show_order_from,list_length);

% show rank result  
[true_rank,pair_rank] = update_truerank_history2(ranking_list, tmp_pairidx,true_rank);
show_truerank_history2(handles.axes_rank1,true_rank);
show_rank_result(all_images,all_pid,all_camid,ranking_list,handles);

% update query_history
query_history = hw_update_query_history(query_history,accumu_idx,current_idx,pair_rank,[]);

% set budget window.
set(handles.text_usedbudget, 'string', sprintf('%3d/%3d',used_budget,label_budget));

% show test performance
if check_on_test
    %*** need to do modification later
    cmc = tst_cmc(curr_metric,tst_set,cls_centers);
    show_test_performance(handles.axes_cmc,cmc);
end
% --- Executes on button press in pushbutton_rerank.
function pushbutton_rerank_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_rerank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rank_history;       global strong_checkid; 
global weak_checkid;       global match_checkid;
global constraint_set;     global current_idx;
global curr_metric;        global all_features;
global incremental_method; global last_step_history;
global used_budget;        global label_budget;
global all_camid;          global all_images;
global all_pid;            global check_on_test;
global show_order_from;    global tst_set;
global label_flg;          global curr_cluster;
global cls_centers;        global true_rank;
global affM;               global query_history;
global accumu_idx;

if (~any(strong_checkid))&&(~any(weak_checkid))&&(~any(match_checkid))
   return; 
end    

% set current image as labelled
label_flg = true;

% push rank_history to last-step-history; - do bug check.
last_step_history.rank_history   = rank_history;
last_step_history.curr_metric    = curr_metric;
last_step_history.constraint_set = constraint_set;
last_step_history.used_budget    = used_budget;
last_step_history.true_rank      = true_rank;

% build constraint set.
cons.strong_id = rank_history(logical(strong_checkid));
cons.weak_id   = rank_history(logical(weak_checkid));
cons.match_id  = rank_history(logical(match_checkid));
constraint_set = hw_update_constraint(current_idx,cons,constraint_set);

% update metric and re-ranking.
% *** currently constraint is updated by set.
% *** maybe try update only by current labels? (i.e. instead of input cons_set, input cons)
% *** graphs?
update_options = hw_make_update_options(incremental_method,curr_metric{curr_cluster},...
cons,constraint_set,current_idx,all_features,all_camid,affM);

curr_metric{curr_cluster} = update_metric_manual(update_options);

% rank current image
ranking_list = rank_current_image(all_features,all_camid,...
               current_idx,curr_metric{curr_cluster},incremental_method);
rank_history = ranking_list;

% clear all check boxes. clear all check-ids.
strong_checkid = zeros(size(rank_history));
weak_checkid   = zeros(size(rank_history));
match_checkid  = zeros(size(rank_history));
show_order_from = 1; list_length = 14;
show_click_result(show_order_from,list_length);

% find true match rank result  
tmp_pairidx = show_next_image(current_idx,handles);
[true_rank,pair_rank] = update_truerank_history2(ranking_list, tmp_pairidx,true_rank); % will update the global true_rank
show_truerank_history2(handles.axes_rank1,true_rank);
show_rank_result(all_images,all_pid,all_camid,ranking_list,handles);

% update re-rank budget, show budget window.
used_budget = used_budget + 1;
set(handles.text_usedbudget, 'string', sprintf('%3d/%3d',used_budget,label_budget));

% update query_history
query_history = hw_update_query_history(query_history,accumu_idx,current_idx,pair_rank,cons);

% accumulate labelling cost. if thr reaches pick next image

% show test performance
if check_on_test
    %*** modeification later
    cmc = tst_cmc(curr_metric,tst_set,cls_centers);
    show_test_performance(handles.axes_cmc,cmc);
end

%% internal algorithms
%%%% begin of internal algorithms

% 1. rank current image;
% function result = rank_current_image(all_features,all_camid,current_idx,curr_metric)

% 2. pick next image;
% function idx = pick_next_image(labelled_set,unlabelled_set,active_method)

% 3. update current metric; 
% *** add graph for labelling propagation?
% function update_metric(current_metric,current_idx,constr_set,all_feats,incremental_method)

% 4. accumulate labelling cost
% function single_labelling_cost(clicks,current_image,graph)

% 5. update constraint set
% function cons_set = update_constraint(idx,newcons,cons_set)

%%%% end of internal algorithms

%% Plotting functions
% ---- Show Next Image ---------------------------------------------------%
function pair_idx = show_next_image(current_idx,handles)
% show next-to-be-labelled image
% find its pairing and show the pairing image if possible
global all_images; 
global all_pid; 
global all_camid;

current_image = all_images{current_idx};
axes(handles.current_img);
imshow(current_image);

current_pid   = all_pid(current_idx);
current_camid = all_camid(current_idx);
display(sprintf('>>> Current idx %04d, pid %04d, camid %02d. ',...
        current_idx, current_pid,current_camid));
% the below only goes with test version. delete in the final version.
% find pairing image if exist.
othercam_mask = (all_camid ~= current_camid);
samepid_mask  = (all_pid == current_pid);
pair_mask     = othercam_mask.*samepid_mask;
pair_idx      = find(pair_mask == 1);
% ---- Clear ranking result --- %
function clear_rank_result(handles)
global hndall;
list_length = 14;
for i = 1: list_length  
    axes(hndall(i));
    im = ones(160,60,3)*220;
    imshow(uint8(im));
    %display(sprintf('> Rankig-list-%02d: pid %04d, camid %02d. ',...
    %    i,pid,camid));
end

% ---- Show ranking result -----------------------------------------------%
function show_rank_result(all_images,all_pid,all_camid,ranking_list,handles)
global hndall;
list_length = 14;
global show_order_from;
for i = 1: list_length  
    idx = ranking_list(show_order_from + i - 1);
    pid = all_pid(idx);
    camid = all_camid(idx);
    im = all_images{idx};
    axes(hndall(i));
    imshow(im);
    %display(sprintf('> Rankig-list-%02d: pid %04d, camid %02d. ',...
    %    i,pid,camid));
end
set(handles.edit_showorder, 'string', sprintf('%d - %d',show_order_from,show_order_from+list_length-1));
% ---- Bounding Sim/disim samples ----------------------------------------%
function strong_check(chkflag,id)
global strong_checkid; global all_images;
global hndall; global rank_history;
global show_order_from;
if chkflag
    strong_checkid(show_order_from+id-1) = 1;
    im = all_images{rank_history(show_order_from+id-1)};
    im = draw_bounding(im, [0, 0, 1]);
    axes(hndall(id));
    imshow(im);
else
    strong_checkid(show_order_from+id-1) = 0;
    im = all_images{rank_history(show_order_from+id-1)};
    axes(hndall(id));
    imshow(im);  
end
function weak_check(chkflag,id)
global weak_checkid; global all_images;
global hndall; global rank_history;
global show_order_from;
if chkflag
    weak_checkid(show_order_from+id-1) = 1;
    im = all_images{rank_history(show_order_from+id-1)};
    im = draw_bounding(im, [1, 0, 0]);
    axes(hndall(id));
    imshow(im);
else
    weak_checkid(show_order_from+id-1) = 0;
    im = all_images{rank_history(show_order_from+id-1)};
    axes(hndall(id));
    imshow(im);  
end
function match_check(chkflag,id)
global match_checkid; global all_images;
global hndall; global rank_history;
global show_order_from;
if chkflag
    match_checkid(show_order_from+id-1) = 1;
    im = all_images{rank_history(show_order_from+id-1)};
    im = draw_bounding(im, [0, 1, 0]);
    axes(hndall(id));
    imshow(im);
else
    match_checkid(show_order_from+id-1) = 0;
    im = all_images{rank_history(show_order_from+id-1)};
    axes(hndall(id));
    imshow(im);  
end
% ---- update check boxes -----------------------------------------------%
function show_click_result(show_order_from, showtopN) 
global strongchkall; global weakchkall;
global matchchkall;  global strong_checkid;
global weak_checkid; global match_checkid;

for i = 1:showtopN
    set(strongchkall(i), 'value', strong_checkid(show_order_from+i-1));
    set(weakchkall(i), 'value', weak_checkid(show_order_from+i-1));
    set(matchchkall(i), 'value', match_checkid(show_order_from+i-1));
end
% ---- show test performance --------------------------------------------%
function show_test_performance(axes_cmc,cmc)
axes(axes_cmc);
hlines = findobj(axes_cmc,'type','line');
set(hlines,'LineWidth',1);
set(hlines,'Color',[1,0.8780,0.8]);%[0.8784,0.8784,0.8784]);
hold on;
num_rank = 50;
plot(1:num_rank ,cmc(1:num_rank),'LineWidth',2,'Color',[0.6,0.6,1]);
axis([1,num_rank , 10,90]);
xlabel('rank');
ylabel('match rate (%)');
grid on;
hold off;
% ---- show true match history ------------------------------------------%
function show_truerank_history(axes_handle,true_rank)
axes(axes_handle);
cla
if ~isempty(true_rank)
max_rank_cur = max(true_rank);
xlim([0, 10]);
%set(gca, 'XTickLabel', {'Intial','1','2','3','4','5'});
ylabel('rank'); xlabel('Round of feedback');
ylim([1, max_rank_cur+10]);
set(gca, 'ytick', 1:floor((max_rank_cur+10)/5):max_rank_cur+10);
hold(axes_handle, 'all');
if ~isempty(true_rank)
    if length(true_rank) < 2
        plot(0:length(true_rank)-1, true_rank, '-rs', 'LineWidth', 2, 'MarkerSize', 6);
    
    else
        plot(0:length(true_rank)-2, true_rank(1:end-1), '--', 'Color', [0.5 0.5 0.5], 'Marker', 's', 'LineWidth', 1.5, 'MarkerSize', 6);
        plot(length(true_rank)-2:length(true_rank)-1, true_rank(end-1:end), '-', 'Color', [1, 0, 0], 'Marker', 's', 'LineWidth', 2, 'MarkerSize', 6);
    end
end
grid on;
end

function show_truerank_history2(axes_handle,true_rank2)
axes(axes_handle);
cla
if ~isempty(true_rank2)
for k = 1:size(true_rank2,2)    
true_rank = true_rank2(:,k)
max_rank_cur = max(true_rank);
xlim([0, 10]);
%set(gca, 'XTickLabel', {'Intial','1','2','3','4','5'});
ylabel('rank'); xlabel('Round of feedback');
ylim([1, max_rank_cur+10]);
set(gca, 'ytick', 1:floor((max_rank_cur+10)/5):max_rank_cur+10);
hold(axes_handle, 'all');
if ~isempty(true_rank)
    if length(true_rank) < 2
        plot(0:length(true_rank)-1, true_rank, '-rs', 'LineWidth', 2, 'MarkerSize', 6);
    
    else
        plot(0:length(true_rank)-2, true_rank(1:end-1), '--', 'Color', [0.5 0.5 0.5], 'Marker', 's', 'LineWidth', 1.5, 'MarkerSize', 6);
        plot(length(true_rank)-2:length(true_rank)-1, true_rank(end-1:end), '-', 'Color', [1, 0, 0], 'Marker', 's', 'LineWidth', 2, 'MarkerSize', 6);
    end
end
end
grid on;
end
% ---- update & display true match ranks --------------------------------%
function [true_rank,pair_rank] = update_truerank_history(ranking_list, tmp_pairidx, true_rank)
[~,pair_rank] = intersect(ranking_list, tmp_pairidx);
if ~isempty(pair_rank)
  display(sprintf('## true match found at %04d.\n', pair_rank));
  true_rank = [true_rank,min(pair_rank)];
else
  display(sprintf('## no true match found.'));
end  

function [true_rank2,pair_rank] = update_truerank_history2(ranking_list, tmp_pairidx, true_rank2)
[~,pair_rank] = intersect(ranking_list, tmp_pairidx);
if ~isempty(pair_rank)
  display(sprintf('## true match found at %04d.\n', pair_rank));
  true_rank2 = [true_rank2; sort(pair_rank)'];
else
  display(sprintf('## no true match found.'));
end  

%% I/O functions.
% ---- save current condition. Intermediate results ---------------------%
function save_epoch(eponum)
display(sprintf('*.*saving temperal historys...epoch %02d*.*',eponum));
global epoch;  
global dataset_name;       global feature_name;  
global incremental_method; global active_method;    
global query_viewid;       global labelled_set;       
global unlabelled_set;     global nGroups;            
global query_history;      global curr_metric; 
global constraint_set;     global model_path;
% make them local instead of global
epoch_t         = epoch;        unlabelled_set_t = unlabelled_set;
labelled_set_t  = labelled_set; query_history_t  = query_history;
curr_metric_t   = curr_metric;  constraint_set_t = constraint_set;
predefined_list = labelled_set;
% save the temperal result
save([model_path,...
sprintf('epoch%02d_%s_%s_%s_%s_QVid_%d_nClster_%d_history_Time_%s.mat',...
eponum,dataset_name,feature_name,incremental_method,active_method,...
query_viewid,nGroups,datestr(now,'HH-MM-mm-dd'))],...
'eponum','epoch_t','labelled_set_t','unlabelled_set_t',...
'query_history_t','curr_metric_t','constraint_set_t',...
'predefined_list','-v7.3');

%% GUI functions.
% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton_resume.
function pushbutton_resume_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_resume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rank_history;      global curr_metric;
global constraint_set;    global used_budget;
global last_step_history; global true_rank;
global all_images;        global all_pid;
global all_camid;         global label_budget;
global query_history;     global accumu_idx;
rank_history = last_step_history.rank_history;
curr_metric  = last_step_history.curr_metric;
constraint_set = last_step_history.constraint_set;
used_budget    = last_step_history.used_budget;
true_rank      = last_step_history.true_rank;
query_history = hw_reset_query_history(query_history,accumu_idx);
show_truerank_history2(handles.axes_rank1,true_rank);
show_rank_result(all_images,all_pid,all_camid,rank_history,handles);
set(handles.text_usedbudget, 'string', sprintf('%3d/%3d',used_budget,label_budget));



% --- Executes on selection change in popupmenu_increment.
function popupmenu_increment_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_increment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_increment
global incremental_method;
contents = cellstr(get(hObject, 'String'));
incremental_method = contents{get(hObject, 'Value')};

% --- Executes during object creation, after setting all properties.
function popupmenu_increment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_active.
function popupmenu_active_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_active contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_active
global active_method;
global predefined_list;
contents = cellstr(get(hObject, 'String'));
active_method = contents{get(hObject, 'Value')};
if strcmp(active_method,'pre-defined')
[filename, pathname] = uigetfile('*.*');
    if ischar(filename) || ischar(pathname)
    tmpvar = load([pathname filename]);
    predefined_list = tmpvar.predefined_list;
    clear tmpvar;
    else
       return;
    end
end

% --- Executes during object creation, after setting all properties.
function popupmenu_active_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_dataset.
function popupmenu_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_dataset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_dataset
global dataset_name;
contents = cellstr(get(hObject, 'String'));
dataset_name = contents{get(hObject, 'Value')};

% --- Executes during object creation, after setting all properties.
function popupmenu_dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_ftrType.
function popupmenu_ftrType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ftrType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ftrType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ftrType
global feature_name;
contents = cellstr(get(hObject, 'String'));
feature_name = contents{get(hObject, 'Value')};

% --- Executes during object creation, after setting all properties.
function popupmenu_ftrType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ftrType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_budget.
function popupmenu_budget_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_budget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_budget contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_budget
global label_budget;
contents = cellstr(get(hObject, 'String'));
label_budget = str2double(contents{get(hObject, 'Value')});

% --- Executes during object creation, after setting all properties.
function popupmenu_budget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_budget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_test.
function checkbox_test_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_test
global check_on_test;
if get(hObject, 'Value')
    check_on_test = 1;
else
    check_on_test = 0;
end

% --- Executes on button press in checkbox_strong1.
function checkbox_strong1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong1
strong_check(get(hObject, 'Value'),1);

% --- Executes on button press in checkbox_strong2.
function checkbox_strong2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong2
strong_check(get(hObject, 'Value'),2);

% --- Executes on button press in checkbox_strong3.
function checkbox_strong3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong3
strong_check(get(hObject, 'Value'),3);

% --- Executes on button press in checkbox_strong4.
function checkbox_strong4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong4
strong_check(get(hObject, 'Value'),4);

% --- Executes on button press in checkbox_strong5.
function checkbox_strong5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong5
strong_check(get(hObject, 'Value'),5);

% --- Executes on button press in checkbox_strong6.
function checkbox_strong6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong6
strong_check(get(hObject, 'Value'),6);

% --- Executes on button press in checkbox_strong7.
function checkbox_strong7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong7
strong_check(get(hObject, 'Value'),7);

% --- Executes on button press in checkbox_strong8.
function checkbox_strong8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong8
strong_check(get(hObject, 'Value'),8);

% --- Executes on button press in checkbox_strong9.
function checkbox_strong9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong9
strong_check(get(hObject, 'Value'),9);

% --- Executes on button press in checkbox_strong10.
function checkbox_strong10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong10
strong_check(get(hObject, 'Value'),10);

% --- Executes on button press in checkbox_strong11.
function checkbox_strong11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong11
strong_check(get(hObject, 'Value'),11);

% --- Executes on button press in checkbox_strong12.
function checkbox_strong12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong12
strong_check(get(hObject, 'Value'),12);

% --- Executes on button press in checkbox_strong13.
function checkbox_strong13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong13
strong_check(get(hObject, 'Value'),13);

% --- Executes on button press in checkbox_strong14.
function checkbox_strong14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_strong14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_strong14
strong_check(get(hObject, 'Value'),14);

% --- Executes on button press in checkbox_weak1.
function checkbox_weak1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak1
weak_check(get(hObject, 'Value'),1);

% --- Executes on button press in checkbox_weak2.
function checkbox_weak2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak2
weak_check(get(hObject, 'Value'),2);

% --- Executes on button press in checkbox_weak3.
function checkbox_weak3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak3
weak_check(get(hObject, 'Value'),3);

% --- Executes on button press in checkbox_weak4.
function checkbox_weak4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak4
weak_check(get(hObject, 'Value'),4);

% --- Executes on button press in checkbox_weak5.
function checkbox_weak5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak5
weak_check(get(hObject, 'Value'),5);

% --- Executes on button press in checkbox_weak6.
function checkbox_weak6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak6
weak_check(get(hObject, 'Value'),6);

% --- Executes on button press in checkbox_weak7.
function checkbox_weak7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak7
weak_check(get(hObject, 'Value'),7);

% --- Executes on button press in checkbox_weak8.
function checkbox_weak8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak8
weak_check(get(hObject, 'Value'),8);

% --- Executes on button press in checkbox_weak9.
function checkbox_weak9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak9
weak_check(get(hObject, 'Value'),9);

% --- Executes on button press in checkbox_weak10.
function checkbox_weak10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak10
weak_check(get(hObject, 'Value'),10);

% --- Executes on button press in checkbox_weak11.
function checkbox_weak11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak11
weak_check(get(hObject, 'Value'),11);

% --- Executes on button press in checkbox_weak12.
function checkbox_weak12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak12
weak_check(get(hObject, 'Value'),12);

% --- Executes on button press in checkbox_weak13.
function checkbox_weak13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak13
weak_check(get(hObject, 'Value'),13);

% --- Executes on button press in checkbox_weak14.
function checkbox_weak14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_weak14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_weak14
weak_check(get(hObject, 'Value'),14);


% --- Executes on button press in checkbox_match1.
function checkbox_match1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match1
match_check(get(hObject, 'Value'),1);

% --- Executes on button press in checkbox_match2.
function checkbox_match2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match2
match_check(get(hObject, 'Value'),2);

% --- Executes on button press in checkbox_match3.
function checkbox_match3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match3
match_check(get(hObject, 'Value'),3);

% --- Executes on button press in checkbox_match4.
function checkbox_match4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match4
match_check(get(hObject, 'Value'),4);

% --- Executes on button press in checkbox_match5.
function checkbox_match5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match5
match_check(get(hObject, 'Value'),5);

% --- Executes on button press in checkbox_match6.
function checkbox_match6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match6
match_check(get(hObject, 'Value'),6);

% --- Executes on button press in checkbox_match7.
function checkbox_match7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match7
match_check(get(hObject, 'Value'),7);

% --- Executes on button press in checkbox_match8.
function checkbox_match8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match8
match_check(get(hObject, 'Value'),8);

% --- Executes on button press in checkbox_match9.
function checkbox_match9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match9
match_check(get(hObject, 'Value'),9);

% --- Executes on button press in checkbox_match10.
function checkbox_match10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match10
match_check(get(hObject, 'Value'),10);

% --- Executes on button press in checkbox_match11.
function checkbox_match11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match11
match_check(get(hObject, 'Value'),11);

% --- Executes on button press in checkbox_match12.
function checkbox_match12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match12
match_check(get(hObject, 'Value'),12);

% --- Executes on button press in checkbox_match13.
function checkbox_match13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match13
match_check(get(hObject, 'Value'),13);

% --- Executes on button press in checkbox_match14.
function checkbox_match14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_match14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_match14
match_check(get(hObject, 'Value'),14);

% --- Executes during object creation, after setting all properties.
function text_usedbudget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_usedbudget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function edit_timer_Callback(hObject, eventdata, handles)
% hObject    handle to edit_timer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_timer as text
%        str2double(get(hObject,'String')) returns contents of edit_timer as a double


% --- Executes during object creation, after setting all properties.
function edit_timer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_timer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in next_page.
function next_page_Callback(hObject, eventdata, handles)
% hObject    handle to next_page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global show_order_from;
global all_images; global all_pid; global all_camid; global rank_history;
list_length = 14; ranking_list = rank_history;
show_order_from = min(show_order_from + list_length, length(ranking_list - list_length + 1));
show_rank_result(all_images,all_pid,all_camid,ranking_list,handles);
show_click_result(show_order_from,list_length);
% --- Executes on button press in previous_page.
function previous_page_Callback(hObject, eventdata, handles)
% hObject    handle to previous_page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global show_order_from;
global all_images; global all_pid; global all_camid; global rank_history;
list_length = 14; ranking_list = rank_history;
show_order_from = max(show_order_from - list_length,1);
show_rank_result(all_images,all_pid,all_camid,ranking_list,handles);
show_click_result(show_order_from,list_length);



function edit_showorder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_showorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_showorder as text
%        str2double(get(hObject,'String')) returns contents of edit_showorder as a double


% --- Executes during object creation, after setting all properties.
function edit_showorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_showorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_nCluster.
function popupmenu_nCluster_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_nCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_nCluster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_nCluster
global nGroups; 
contents = cellstr(get(hObject, 'String'));
nGroups = str2double(contents{get(hObject, 'Value')});

% --- Executes during object creation, after setting all properties.
function popupmenu_nCluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_nCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
clearvars -global;



function edit_accumulate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_accumulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_accumulate as text
%        str2double(get(hObject,'String')) returns contents of edit_accumulate as a double


% --- Executes during object creation, after setting all properties.
function edit_accumulate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_accumulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


