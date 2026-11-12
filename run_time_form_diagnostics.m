function results=run_time_form_diagnostics(nStarts,baselineDir)
% Paper-inspired time-form diagnostics ONLY; never re-fit archived models.
% Primary-source definitions must be reviewed before this entry can run.
if nargin<1,nStarts=8;end
validateattributes(nStarts,{'numeric'},{'scalar','integer','>=',1});
root=fileparts(mfilename('fullpath'));
if nargin<2
 baselineDir=fullfile(root,'baseline-results',...
  'f05543e978e425785fdf918751185cd1a733d300','43137925939-1','joint_results');
end
paper=jsondecode(fileread(fullfile(root,'reports','paper_time_definitions.json')));
requirePaper(paper); % No arbitrary origin can silently enter a formal fit.
base=load(fullfile(baselineDir,'delay_study.mat'),'cfg','g','cs');
cfg=base.cfg;g=base.g;cs=base.cs;
assert(isequal(g,[0.031094959444983555,35143.615284616222]),'PKPD:GrowthReference','Growth reference changed.');
assert(cfg.rtol==1e-7 && cfg.atol==1e-9,'PKPD:Tolerances','Archived ODE precision changed.');
manifest=jsondecode(fileread(fullfile(fileparts(baselineDir),'manifest.json')));
assert(strcmp(manifest.source_sha,'f05543e978e425785fdf918751185cd1a733d300') && ...
 strcmp(manifest.run_id,'43137925939') && strcmp(manifest.matlab.status,'success') && ...
 strcmp(manifest.mode,'delay-full') && strcmp(manifest.archive_sha256,...
 'b39c64c22d8ebe1e9fd6bcd07d42c478c2b105e0d7cfdb1df8f8d1c2d5f593d3'),...
 'PKPD:BaselineProvenance','Unexpected archived baseline.');
out=fullfile(root,'joint_results');if ~exist(out,'dir'),mkdir(out);end
assert(~exist(fullfile(out,'time_form_inputs.mat'),'file'),'PKPD:ExistingTimeFit','Time-form outputs already exist.');
previous=get(groot,'DefaultFigureVisible');set(groot,'DefaultFigureVisible','off');
figureCleanup=onCleanup(@() set(groot,'DefaultFigureVisible',previous)); %#ok<NASGU>
diary(fullfile(out,'run.log'));diaryCleanup=onCleanup(@() diary('off')); %#ok<NASGU>
timer=tic;
fprintf('TIME-FORM ONLY: archived D, dose-beta and delay predictions are read, never re-fitted.\n');
copyfile(fullfile(root,'reports','paper_time_definitions.json'),out);
copyfile(fullfile(fileparts(baselineDir),'manifest.json'),fullfile(out,'baseline_manifest.json'));
save(fullfile(out,'time_form_inputs.mat'),'cfg','g','cs','paper','manifest','nStarts');
checks=numericalChecks(g,cfg,out); % Synthetic clocks, explicitly separate from experimental data.

origin=paper.origin_in_current_data_hours;doseStart=paper.dose_start_in_current_data_hours;
params=[];startRows=[];newPred=[];newResid=[];
lowerKD=1e-9;upperKD=1;lowerR=1e-4;upperR=100;
boundTable=table({'kd_A';'kd_B';'r'},{'h^-1';'h^-1';'h^-1'},...
 [lowerKD;lowerKD;lowerR],[upperKD;upperKD;upperR],[true;true;false],...
 'VariableNames',{'parameter','unit','positive_lower','upper','zero_also_compared'});
writetable(boundTable,fullfile(out,'time_form_bounds.csv'));
for kind=1:2
 for j=1:numel(cs)
  c=cs(j);rng(57724+100*kind+j,'twister');
  folder=fullfile(out,'time_form_checkpoints',sprintf('form_%d_dose_%d',kind,j));mkdir(folder);
  guess=0.01;lo=lowerKD;hi=upperKD;
  if kind==2,guess=[0.02,0.01];lo=[lowerKD,lowerR];hi=[upperKD,upperR];end
  fun=@(p) timeLoss(p,kind,c,g,origin,cfg);
  [p,sse,rr]=multifit(fun,guess,lo,hi,nStarts,sprintf('Time form %d dose %g',kind,c.dose),folder);
  % Exact zero is admissible even though positive searches use log coordinates.
  zero=guess;zero(1)=0;zeroSSE=fun(zero);selectedZero=zeroSSE<sse;
  if selectedZero,p=zero;sse=zeroSSE;end
  r=NaN;peak=NaN;nearRLow=false;nearRHigh=false;
  if kind==2,r=p(2);peak=1/r;nearRLow=log(r/lowerR)<0.02;nearRHigh=log(upperR/r)<0.02;end
  if kind==2 && selectedZero,r=NaN;peak=NaN;nearRLow=false;nearRHigh=false;end
  nearKDLow=p(1)==0 || log(p(1)/lowerKD)<0.02;nearKDHigh=log(upperKD/p(1))<0.02;
  t=unique([linspace(c.t(1),max(c.t),1501)';c.t;...
   max(c.t(1),min(max(c.t),origin));max(c.t(1),min(max(c.t),doseStart))]);
  if isfinite(peak) && origin+peak>=c.t(1) && origin+peak<=max(c.t),t=unique([t;origin+peak]);end
  [N,kill]=simulateTime(p,kind,t,c,g,origin,cfg);
  obs=simulateTime(p,kind,c.t,c,g,origin,cfg);res=obs-c.y;
  assert(all(N>0) && all(kill>=0),'PKPD:TimeStates','Best fit has nonpositive N or negative killing.');
  model=kind+2;
  params=[params;model,c.dose,p(1),r,peak,sse,sqrt(mean(res.^2)),...
   zeroSSE,selectedZero,nearKDLow,nearKDHigh,nearRLow,nearRHigh]; %#ok<AGROW>
  if kind==1,rr=[rr(:,1),nan(nStarts,1),rr(:,2:end)];end
  startRows=[startRows;repmat([model,c.dose],nStarts,1),(1:nStarts)',rr]; %#ok<AGROW>
  newPred=[newPred;repmat([model,c.dose],numel(t),1),t,N,kill]; %#ok<AGROW>
  newResid=[newResid;repmat([model,c.dose],numel(c.t),1),c.t,c.y,obs,res,res/max(c.y)]; %#ok<AGROW>
  save(fullfile(folder,'fit.mat'),'p','sse','rr','zeroSSE','selectedZero','c','paper','g','cfg','lo','hi');
  writetable(parameterTable(params),fullfile(out,'time_form_parameters.csv'));
  writetable(startsTable(startRows),fullfile(out,'time_form_starts.csv'));
 end
end

% Reuse the exact stored trajectories/residuals, not fresh old-model simulations.
oldPred=readtable(fullfile(baselineDir,'delay_predictions.csv'));
oldResid=readtable(fullfile(baselineDir,'delay_residuals.csv'));
oldIds=[1,2,5,6];oldPred=oldPred(ismember(oldPred.model_id,oldIds),:);
oldResid=oldResid(ismember(oldResid.model_id,oldIds),:);
predictions=[oldPred{:,{'model_id','dose_nM','time_h','N','killing_rate'}};newPred];
residuals=[oldResid{:,:};newResid];
predictionTable=array2table(predictions,'VariableNames',{'model_id','dose_nM','time_h','N','killing_rate'});
predictionTable.time_since_dose_h=predictionTable.time_h-doseStart;
predictionTable.tau_h=predictionTable.time_h-origin;
writetable(predictionTable,fullfile(out,'time_form_predictions.csv'));
writetable(array2table(residuals,'VariableNames',{'model_id','dose_nM','time_h',...
 'observed_N','predicted_N','residual','scaled_residual'}),fullfile(out,'time_form_residuals.csv'));
names={'Archived linear D';'Archived dose-specific beta';'Paper-inspired constant-rate diagnostic';...
 'Paper-inspired transient-rate diagnostic';'Archived delay stage 1';'Archived delay stage 2'};
scores=zeros(6,1);byDose=[];
for m=1:6
 scores(m)=sum(residuals(residuals(:,1)==m,7).^2);
 for j=1:numel(cs)
  rows=residuals(residuals(:,1)==m & residuals(:,2)==cs(j).dose,:);
  byDose=[byDose;m,cs(j).dose,sum(rows(:,7).^2),sqrt(mean(rows(:,6).^2))]; %#ok<AGROW>
 end
end
summary=table((1:6)',names,scores,[0;0;3;6;0;0],[2;3;0;0;2;3],...
 'VariableNames',{'model_id','model','scaled_SSE','parameters_refitted_this_comparison','previous_fit_parameters'});
writetable(summary,fullfile(out,'time_form_summary.csv'));
writetable(array2table(byDose,'VariableNames',{'model_id','dose_nM','scaled_SSE','RMSE_cells'}),...
 fullfile(out,'time_form_by_dose.csv'));
for name={'delay_parameters.csv','delay_by_dose.csv','delay_summary.csv'}
 copyfile(fullfile(baselineDir,name{1}),fullfile(out,['baseline_',name{1}]));
end
plotComparisons(predictions,residuals,cs,out);
results=struct('summary',summary,'parameters',parameterTable(params),'starts',startsTable(startRows),...
 'cfg',cfg,'growth',g,'checks',checks,'baseline',manifest,'paper',paper,'elapsed_seconds',toc(timer));
save(fullfile(out,'time_form_fit.mat'),'results','predictions','residuals','byDose');
fprintf('TIME-FORM COMPLETE: A SSE=%.10g, B SSE=%.10g, %.1f s.\n',scores(3),scores(4),toc(timer));
end

function requirePaper(paper)
assert(isfield(paper,'verified') && isequal(paper.verified,true),...
 'PKPD:UnverifiedPaper','Primary definitions are unverified; no formal A/B fit is permitted.');
for key={'origin_in_current_data_hours','dose_start_in_current_data_hours'}
 assert(isfield(paper,key{1}) && isnumeric(paper.(key{1})) && isscalar(paper.(key{1})) && ...
  isfinite(paper.(key{1})),'PKPD:PaperClock','Verified data/dose/effect clock mapping is missing.');
end
for key={'time_origin_definition','D_definition','theta_definition','initial_cell_number_definition','fitting_interval_definition'}
 k=key{1};assert(isfield(paper,k) && ischar(paper.(k)) && ~isempty(strtrim(paper.(k))),...
  'PKPD:PaperDefinition','Missing primary definition: %s',k);
 assert(isfield(paper,'source_evidence') && isfield(paper.source_evidence,k),...
  'PKPD:PaperEvidence','Missing source evidence for %s',k);
 evidence=paper.source_evidence.(k);
 for field={'source','location','excerpt'}
  assert(isfield(evidence,field{1}) && ischar(evidence.(field{1})) && ~isempty(strtrim(evidence.(field{1}))),...
   'PKPD:PaperEvidence','Record source, section/page and excerpt for %s',k);
 end
end
end

function T=parameterTable(rows)
T=array2table(rows,'VariableNames',{'model_id','dose_nM','kd_per_h','r_per_h','kernel_peak_tau_h',...
 'scaled_SSE','RMSE_cells','zero_killing_SSE','selected_zero','kd_near_lower','kd_near_upper',...
 'r_near_lower','r_near_upper'});
end

function T=startsTable(rows)
T=array2table(rows,'VariableNames',{'model_id','dose_nM','start','kd_per_h','r_per_h',...
 'scaled_SSE','exitflag','function_evaluations'});
end

function value=timeLoss(p,kind,c,g,origin,cfg)
N=simulateTime(p,kind,c.t,c,g,origin,cfg);
value=sum(((N-c.y)/max(c.y)).^2);
end

function [N,effect,transition]=simulateTime(p,kind,t,c,g,origin,cfg)
t=t(:);initial=c.t(1);assert(all(t>=initial) && c.y(1)>0);
N=zeros(size(t));N(t==initial)=c.y(1);effect=timeRate(t,p,kind,origin);
transition=struct('jump',0);z=log(c.y(1));stop=max(t);
breaks=unique([initial;max(initial,min(stop,origin));stop]);
opts=odeset('RelTol',cfg.rtol,'AbsTol',cfg.atol);
for k=1:numel(breaks)-1
 left=breaks(k);right=breaks(k+1);
 if right<=origin,rate=@(~) 0;else,rate=@(u) timeRate(u,p,kind,origin);end
 sol=ode15s(@(u,y) g(1)*(1-exp(y)/g(2))-rate(u),[left,right],z,opts);
 transition.jump=max(transition.jump,abs(sol.y(1,1)-z));
 mask=t>left & t<=right;N(mask)=exp(deval(sol,t(mask)))';z=deval(sol,right);
end
assert(all(isfinite(N)),'PKPD:TimeFinite','Nonfinite time-form solution.');
end

function value=timeRate(t,p,kind,origin)
value=zeros(size(t));tau=t-origin;active=tau>=0;
if p(1)==0,return;end
if kind==1,value(active)=p(1);
else,x=p(2)*tau(active);value(active)=p(1)*x.*exp(1-x);end
end

function checks=numericalChecks(g,cfg,out)
% Deliberately synthetic time axes; these do not set any experimental clock.
c=struct('t',[5;6;9;25;105],'y',[4000;4000;4000;4000;4000]);
t=unique([linspace(5,105,501)';9]);origin=9;
expected=g(2)./(1+(g(2)/c.y(1)-1).*exp(-g(1)*(t-c.t(1))));
rows=[];
for kind=1:2
 p=0;if kind==2,p=[0,0.1];end
 N=simulateTime(p,kind,t,c,g,origin,cfg);err=max(abs(N-expected))/max(expected);
 assert(err<1e-5,'PKPD:TimeLogistic','Zero killing does not recover logistic growth.');
 p=0.02;if kind==2,p=[0.02,0.1];end
 [N,kd,tr]=simulateTime(p,kind,t,c,g,origin,cfg);
 before=t<=origin;preError=max(abs(N(before)-expected(before)))/max(expected(before));
 assert(preError<1e-5 && all(kd(t<origin)==0) && tr.jump==0,...
  'PKPD:TimeOrigin','Pre-effect observations or state continuity mishandled.');
 assert(all(N>0) && all(kd>=0),'PKPD:TimeNonnegative','Nonpositive N or negative killing in synthetic check.');
 rows=[rows;kind,err,preError,tr.jump,min(N),min(kd)]; %#ok<AGROW>
end
peaks=[];
for r=[1e-4,0.01,100]
 kd=0.02;p=[kd,r];peak=origin+1/r;
 vals=timeRate([origin-1,origin,peak-0.01/r,peak,peak+0.01/r],p,2,origin);
 assert(vals(1)==0 && vals(2)==0 && abs(vals(4)-kd)<1e-12 && vals(4)>max(vals([3,5])),...
  'PKPD:TimePeak','B peak is not kd_B at tau=1/r.');
 peaks=[peaks;r,1/r,vals(4),abs(vals(4)-kd)]; %#ok<AGROW>
end
% A has a closed-form solution; check onset before, inside and after the window.
constantRows=[];
for onset=[-5,9,200]
 kd=0.02;[N,~,tr]=simulateTime(kd,1,t,c,g,onset,cfg);
 before=min(t,max(onset,c.t(1)))-c.t(1);
 atOnset=g(2)./(1+(g(2)/c.y(1)-1).*exp(-g(1)*before));
 dt=max(t-max(onset,c.t(1)),0);q=g(1)-kd;
 exact=1./(exp(-q*dt)./atOnset+(g(1)/g(2))/q*(1-exp(-q*dt)));
 err=max(abs(N-exact))/max(exact);
 assert(err<1e-5 && tr.jump==0,'PKPD:TimeClock','Constant-rate analytic clock check failed.');
 constantRows=[constantRows;onset,err]; %#ok<AGROW>
end
checks=struct('limits',rows,'peaks',peaks,'constant',constantRows,'synthetic_initial_time_h',5);
writetable(array2table(rows,'VariableNames',{'form','logistic_relative_error','pre_effect_relative_error',...
 'logN_state_jump','minimum_N','minimum_killing_rate'}),fullfile(out,'time_form_numerical_checks.csv'));
writetable(array2table(peaks,'VariableNames',{'r_per_h','peak_tau_h','peak_rate','peak_rate_error'}),...
 fullfile(out,'time_form_peak_checks.csv'));
writetable(array2table(constantRows,'VariableNames',{'synthetic_origin_h','relative_error'}),...
 fullfile(out,'time_form_clock_checks.csv'));
fprintf('TIME-FORM SYNTHETIC CHECKS PASSED: logistic, peak, clocks, continuity, nonnegative states.\n');
end

function plotComparisons(pred,res,cs,out)
colors=lines(6);labels={'D','Dose beta','Time A','Time B','Delay 1','Delay 2'};
f=figure('Visible','off','Position',[50,50,1400,800]);
for j=1:numel(cs)
 subplot(2,3,j);hold on;
 for m=1:6,q=pred(pred(:,1)==m & pred(:,2)==cs(j).dose,:);plot(q(:,3),q(:,4),'Color',colors(m,:));end
 plot(cs(j).t,cs(j).y,'ko','MarkerFaceColor','k');title(sprintf('%g nM',cs(j).dose));xlabel('Data time (h)');ylabel('Cells');grid on;
 if j==1,legend([labels,{'Digitized data'}],'Location','best');end
 subplot(2,3,j+3);hold on;
 for m=1:6,q=res(res(:,1)==m & res(:,2)==cs(j).dose,:);plot(q(:,3),q(:,7),'o-','Color',colors(m,:));end
 yline(0,':');xlabel('Data time (h)');ylabel('Scaled residual: prediction - data');grid on;
end
saveas(f,fullfile(out,'time_form_comparison.png'));close(f);
f=figure('Visible','off','Position',[50,50,1300,400]);
for j=1:numel(cs)
 subplot(1,3,j);hold on;
 for m=3:4,q=pred(pred(:,1)==m & pred(:,2)==cs(j).dose,:);plot(q(:,3),q(:,5),'Color',colors(m,:));end
 title(sprintf('%g nM: time-form diagnostics',cs(j).dose));xlabel('Data time (h)');ylabel('Killing rate (h^{-1})');grid on;legend('Time A','Time B');
end
saveas(f,fullfile(out,'time_form_killing.png'));close(f);
end

% Original optimizer copied unchanged to keep archived entry points untouched.
function [best,fbest,runs]=multifit(fun,guess,lo,hi,n,label,checkpointDir)
guess=guess(:)';lo=lo(:)';hi=hi(:)';d=numel(guess);
opts=optimset('Display','off','MaxIter',650,'MaxFunEvals',1800,...
 'TolX',1e-5,'TolFun',1e-8,'OutputFcn',@progress);
runs=zeros(n,d+3);fbest=Inf;best=guess;
for i=1:n
 startTimer=tic;lastProgress=0;
 fprintf('[%s] %s %d/%d START\n',datestr(now,30),label,i,n);
 if i==1,z=log(guess);
 elseif mod(i,2)==0,z=min(log(hi)-1e-5,max(log(lo)+1e-5,log(guess)+2*randn(1,d)));
 else,z=log(lo)+rand(1,d).*(log(hi)-log(lo));end
 [z,~,flag,output]=fminsearch(@(z) boundedLoss(z,fun,lo,hi),z,opts);
 p=exp(z);sse=boundedLoss(z,fun,lo,hi);
 runs(i,:)=[p,sse,flag,output.funcCount];
 if sse<fbest,fbest=sse;best=p;end
 elapsedSeconds=toc(startTimer);completedRuns=runs(1:i,:);rngState=rng;
 % One file per completed start: later failures cannot overwrite earlier starts.
 checkpoint=fullfile(checkpointDir,sprintf('%s_start_%02d.mat',...
  matlab.lang.makeValidName(label),i));
 save(checkpoint,'label','i','n','p','sse','flag','output','completedRuns',...
  'best','fbest','guess','lo','hi','elapsedSeconds','rngState');
 fprintf('[%s] %s %d/%d DONE: SSE=%.7g, exit=%d, evals=%d, %.1f s; %s\n',...
  datestr(now,30),label,i,n,sse,flag,output.funcCount,elapsedSeconds,checkpoint);
 drawnow;
end
if fbest>=1e9,error('No valid fit found for %s.',label);end
if all(runs(:,end-1)<=0),warning('%s: no converged start; inspect results.',label);end
 function stop=progress(~,optimValues,state)
  stop=false;
  if strcmp(state,'iter') && toc(startTimer)-lastProgress>=60
   lastProgress=toc(startTimer);
   fprintf('[%s] %s %d/%d running: iter=%d, evals=%d, SSE=%.7g, %.1f s\n',...
    datestr(now,30),label,i,n,optimValues.iteration,optimValues.funccount,...
    optimValues.fval,lastProgress);
  end
 end
end
function v=boundedLoss(z,fun,lo,hi)
if any(~isfinite(z)),v=1e12;return;end
violation=max(log(lo)-z,0)+max(z-log(hi),0);
if any(violation>0),v=1e9+sum(violation.^2);return;end
try,v=fun(exp(z));if ~isfinite(v),v=1e9;end
catch,v=1e9;end
end
