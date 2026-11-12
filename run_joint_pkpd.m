function results = run_joint_pkpd(nStarts,doSensitivity,doDoseBeta,doSaturation,delayMode)
% RUN_JOINT_PKPD  SUM-149PT reversible PK / direct cell-killing hypothesis.
% Open this file and click Run, or type: results = run_joint_pkpd;
% Base MATLAB only, no extra toolbox. Approximate digitized data are embedded.
% Optional quick numerical check: run_joint_pkpd(2). Default = 8 starts/model.
%
% CF' = kEF*(VE/VI)*CE - (kFE+kFB)*CF + kBF*CB
% CB' = kFB*CF - kBF*CB
% N'  = kp*N*(1-N/theta) - beta*CB*N
% NO damage-memory, NO concentration-dependent theta.
%
% This version jointly adjusts kFE,kFB,kBF,beta. kEF and volumes are fixed
% scale conventions; kp,theta are estimated from control and then fixed.
% IMPORTANT: even these four parameters are NOT individually identifiable
% from N alone. Exact equivalent parameter sets are generated below.
% Old model PK trajectories are comparison references, NOT experimental data
% and NOT an additional fitting objective. No paper-PD reproduction is claimed.
%
% Sources: article Fig4 SUM-149PT column, supplement FigS3 SUM-149PT column,
% McKenna et al. Sci Rep 7:5725 (2017), DOI 10.1038/s41598-017-05902-z.
% Sparse manual digitization, not raw data. See README_zh.md for limitations.
if nargin<1,nStarts=8;end
if nargin<2,doSensitivity=false;end
if nargin<3,doDoseBeta=false;end
if nargin<4,doSaturation=false;end
if nargin<5,delayMode='';end
if doSaturation,doDoseBeta=true;end
if ~isempty(delayMode),delayMode=validatestring(delayMode,{'quick','full'});doDoseBeta=true;end
validateattributes(nStarts,{'numeric'},{'scalar','integer','>=',1});
rng(42,'twister'); root=fileparts(mfilename('fullpath'));
out=fullfile(root,'joint_results');if ~exist(out,'dir'),mkdir(out);end
checkpointDir=fullfile(out,'checkpoints');
if ~exist(checkpointDir,'dir'),mkdir(checkpointDir);end
diary(fullfile(out,'run.log')); diary on;
logCleanup=onCleanup(@() diary('off')); %#ok<NASGU>
runTimer=tic;
fprintf('[%s] START run_joint_pkpd(%d), seed=42/twister, MATLAB %s, %s\n',...
 datestr(now,30),nStarts,version,computer);
cfg=struct('cellLine','SUM-149PT','VE',250,'VI',0.004629,'kEF',3.5e-6,...
    'oldFE',1.25,'oldFB',0.21,'rtol',1e-7,'atol',1e-9);
[C,cs]=embeddedData();
for j=1:numel(cs)
 assert(all(diff(cs(j).t)>0) && cs(j).t(1)==0 && all(cs(j).y>0));
end
fprintf('SUM-149PT, 156.25/312.5/625 nM, 6 h exposure.\n');
fprintf('APPROXIMATE DIGITIZED DATA; PK reference is simulation only.\n');
initialRng=rng;
save(fullfile(out,'run_inputs.mat'),'cfg','C','cs','nStarts','initialRng');

%% 1. Untreated growth: estimate kp and theta, then hold fixed.
gfun=@(p) sum(((logistic(p,C(:,1),C(1,2))-C(:,2))/max(C(:,2))).^2);
[g,gSSE,gRuns]=multifit(gfun,[0.03,35000],[1e-5,1e3],[1,1e6],...
 min(4,nStarts),'control',checkpointDir);

%% 2. Three comparators, all fitted to exactly the same cell-count observations.
% p is always [kFE,kFB,kBF,beta]. Different parameter spaces are mapped to p.
lo=[1e-4,1e-4,1e-7,1e-9]; hi=[10,10,10,0.1];
% A: old irreversible PK fixed; fit beta only. This is NOT original paper PD.
funA=@(x) loss([cfg.oldFE,cfg.oldFB,0,x],cs,cfg,g);
[xA,sseA,runsA]=multifit(funA,5e-4,lo(4),hi(4),min(4,nStarts),'old-PK/direct-PD',checkpointDir);
pA=[cfg.oldFE,cfg.oldFB,0,xA];
% B: irreversible model; allow both old kinetic rates and beta to refit.
funB=@(x) loss([x(1),x(2),0,x(3)],cs,cfg,g);
[xB,sseB,runsB]=multifit(funB,[cfg.oldFE,cfg.oldFB,xA],lo([1,2,4]),...
 hi([1,2,4]),nStarts,'joint irreversible',checkpointDir);pB=[xB(1:2),0,xB(3)];
% C: reversible model; jointly adjust kFE,kFB,kBF,beta.
funC=@(x) loss(x,cs,cfg,g);
[pC,sseC,runsC]=multifit(funC,[cfg.oldFE,cfg.oldFB,0.005,5e-4],lo,hi,...
 nStarts,'joint reversible',checkpointDir);
% Check a near-irreversible start as an additional nested-model sanity check.
[pExtra,sExtra,rExtra]=multifit(funC,[pB(1:2),lo(3)*1.01,pB(4)],lo,hi,1,...
 'reversible near zero',checkpointDir);runsC=[runsC;rExtra];
if sExtra<sseC,pC=pExtra;sseC=sExtra;end
if sseC>sseB+1e-4
 warning('Reversible fit worse than irreversible. Check convergence/search bounds.');
end

%% 2D. Fixed-reference PK; fit only shared kBF and beta.
% Exact growth doubles from full run 43100077010 (commit 39cf75b).
% Reference PK values remain the original cfg values, not experimental estimates.
gD=[0.031094959444983555,35143.615284616222];
fixedNames={'kEF';'kFE';'kFB';'VI';'VE';'kp';'theta'};
fixedValues=[cfg.kEF;cfg.oldFE;cfg.oldFB;cfg.VI;cfg.VE;gD'];
fixedUnits={'h^-1';'h^-1';'h^-1';'uL';'uL';'h^-1';'cells'};
fixedSources=[repmat({'6325d45:run_joint_pkpd.m cfg'},5,1);...
 repmat({'full run 43100077010 control fit, commit 39cf75b'},2,1)];
fixedTable=table(fixedNames,fixedValues,fixedUnits,fixedSources,...
 'VariableNames',{'parameter','fixed_value','unit','source'});
writetable(fixedTable,fullfile(out,'model_d_fixed_parameters.csv'));
funD=@(x) loss([cfg.oldFE,cfg.oldFB,x(1),x(2)],cs,cfg,gD);
[xD,sseD,runsD]=multifit(funD,[0.005,5e-4],lo([3,4]),hi([3,4]),...
 nStarts,'fixed reference D',checkpointDir);
pD=[cfg.oldFE,cfg.oldFB,xD];
writeRuns(runsD,{'kBF','beta'},out,'model_d');
save(fullfile(out,'model_d_fit.mat'),'pD','sseD','runsD','gD','fixedTable','cfg','C','cs');
profileD=profileFixedPK(pD,sseD,cs,cfg,gD,lo,hi,out);

%% 3. Exact non-identifiability demonstration, not a confidence interval.
% E=beta*CB satisfies E''+a*E'+b*E=c*CE for zero initial intracellular drug:
% a=kFE+kFB+kBF; b=kFE*kBF; c=beta*kEF*(VE/VI)*kFB.
% Keeping (a,b,c) fixed preserves E(t), hence N(t), for ANY prescribed input.
% Thus fixing kEF alone does not identify all remaining four parameters.
abc=combinations(pC,cfg);
eq=equivalentFamily(pC,lo,hi); % pC plus up to four distinct representatives
for j=1:size(eq,1)
 assert(max(abs(combinations(eq(j,:),cfg)./abc-1))<1e-7);
end

%% 4. Parameter and fit summaries. 'Representative' does NOT mean unique.
model={'Old PK, direct PD';'Joint irreversible';'Joint reversible';'Fixed reference D'};
P=[pA;pB;pC;pD];modelGrowth=[g;g;g;gD];
summary=table(model,P(:,1),P(:,2),P(:,3),P(:,4),[sseA;sseB;sseC;sseD],...
 'VariableNames',{'model','kFE','kFB','kBF','beta','scaled_SSE'});
disp(summary);
growth=table(g(1),g(2),gSSE,'VariableNames',{'kp','theta','scaled_SSE'});
writetable(summary,fullfile(out,'representative_parameters.csv'));
writetable(growth,fullfile(out,'growth_parameters.csv'));
writetable(array2table(eq,'VariableNames',{'kFE','kFB','kBF','beta'}),...
 fullfile(out,'equivalent_parameter_sets.csv'));
writetable(array2table(abc,'VariableNames',{'a_rate_sum','b_rate_product','c_effect_gain'}),...
 fullfile(out,'reversible_combinations.csv'));
writeRuns(runsA,{'beta'},out,'old_pk');
writeRuns(runsB,{'kFE','kFB','beta'},out,'joint_irreversible');
writeRuns(runsC,{'kFE','kFB','kBF','beta'},out,'joint_reversible');
% The bounds below are computational exploration limits, not biological priors.
writetable(table({'kFE';'kFB';'kBF';'beta'},lo',hi',...
 'VariableNames',{'parameter','lower','upper'}),fullfile(out,'search_bounds.csv'));
nearBound=any((log(pC)-log(lo))<0.02 | (log(hi)-log(pC))<0.02);
if nearBound,warning('At least one representative parameter is near a bound.');end

%% 5. Cell fits and residuals; no paper fit line is used as target.
f=figure('Visible','off','Color','w','Position',[80,80,1100,760]);
subplot(2,2,1);tt=linspace(0,max(C(:,1)),300)';
plot(C(:,1),C(:,2),'ko',tt,logistic(g,tt,C(1,2)),'k-','LineWidth',1.4);
title('SUM-149PT untreated');xlabel('Hours');ylabel('Cell count');grid on;
cols=[0.5,0.5,0.5;0.85,0.25,0.15;0.1,0.35,0.85;0.1,0.6,0.3];styles={':','--','-','-.'};
predRows=[];resRows=[];byDose=[];
for j=1:numel(cs)
 c=cs(j);tt=unique([linspace(0,max(c.t),350)';c.wash]);
 subplot(2,2,j+1);plot(c.t,c.y,'ko','MarkerSize',5);hold on;
 for m=1:4
  Y=simulate(P(m,:),tt,c,cfg,modelGrowth(m,:));
  plot(tt,Y(:,3),styles{m},'Color',cols(m,:),'LineWidth',1.5);
  predRows=[predRows;[repmat([m,c.dose],numel(tt),1),tt,Y]]; %#ok<AGROW>
  Yobs=simulate(P(m,:),c.t,c,cfg,modelGrowth(m,:));rr=Yobs(:,3)-c.y;
  resRows=[resRows;[repmat([m,c.dose],numel(c.t),1),c.t,c.y,Yobs(:,3),rr,rr/max(c.y)]]; %#ok<AGROW>
  byDose=[byDose;[m,c.dose,sum((rr/max(c.y)).^2),sqrt(mean(rr.^2))]]; %#ok<AGROW>
 end
 title(sprintf('SUM-149PT, %g nM, 6 h',c.dose));xlabel('Hours from treatment start');ylabel('Cell count');grid on;
 legend({'Digitized points','Old PK/direct PD','Joint irreversible','Joint reversible','Fixed reference D'},'Location','best');
end
saveas(f,fullfile(out,'01_cell_fits.png'));close(f);
fprintf('[%s] Saved 01_cell_fits.png (total %.1f s)\n',datestr(now,30),toc(runTimer));
writetable(array2table(predRows,'VariableNames',{'model_id','dose_nM','time_h','CF','CB','N'}),...
 fullfile(out,'predictions.csv'));
writetable(array2table(resRows,'VariableNames',{'model_id','dose_nM','time_h','observed_N','predicted_N','residual','scaled_residual'}),...
 fullfile(out,'residuals.csv'));
writetable(array2table(byDose,'VariableNames',{'model_id','dose_nM','scaled_SSE','RMSE_cells'}),...
 fullfile(out,'fit_by_dose.csv'));
f=figure('Visible','off','Color','w','Position',[80,80,1100,360]);
for j=1:3
 subplot(1,3,j);hold on;
 for m=1:4
  Q=resRows(resRows(:,1)==m & resRows(:,2)==cs(j).dose,:);
  plot(Q(:,3),Q(:,7),'o-','Color',cols(m,:));
 end
 plot([0,max(cs(j).t)],[0,0],'k:');title(sprintf('%g nM',cs(j).dose));
 xlabel('Hours');ylabel('Scaled residual: prediction - data');grid on;
end
saveas(f,fullfile(out,'02_residuals.png'));close(f);
fprintf('[%s] Saved 02_residuals.png (total %.1f s)\n',datestr(now,30),toc(runTimer));

%% 6. PK predictions versus old simulated PK. These are NOT PK observations.
f=figure('Visible','off','Color','w','Position',[80,80,1100,650]);
for j=1:3
 c=cs(j);tt=unique([linspace(0,24,120)';linspace(24,600,240)']);
 old=simulate(pA,tt,c,cfg,g);new=simulate(pC,tt,c,cfg,g);
 subplot(2,3,j);plot(tt,sum(old(:,1:2),2),'k--',tt,sum(new(:,1:2),2),'b-');
 title(sprintf('%g nM: total DOX',c.dose));xlabel('Hours');ylabel('CF+CB (nM)');grid on;
 legend({'Old simulated PK','One joint-fit representative'},'Location','best');
 subplot(2,3,j+3);plot(tt,old(:,2),'k--',tt,new(:,2),'b-');
 title(sprintf('%g nM: bound DOX',c.dose));xlabel('Hours');ylabel('CB (nM)');grid on;
end
saveas(f,fullfile(out,'03_PK_reference.png'));close(f);
fprintf('[%s] Saved 03_PK_reference.png (total %.1f s)\n',datestr(now,30),toc(runTimer));

%% 7. Distinct PK parameter sets can predict the same cell counts.
c=cs(3);tt=unique([linspace(0,24,120)';linspace(24,600,240)']);
f=figure('Visible','off','Color','w','Position',[80,80,1100,360]);
labels=cell(size(eq,1),1);ref=simulate(eq(1,:),tt,c,cfg,g);maxDelta=0;
for j=1:size(eq,1)
 Y=simulate(eq(j,:),tt,c,cfg,g);labels{j}=sprintf('set %d: kBF=%.3g',j,eq(j,3));
 maxDelta=max(maxDelta,max(abs(Y(:,3)-ref(:,3)))/max(ref(:,3)));
 subplot(1,3,1);hold on;plot(tt,sum(Y(:,1:2),2),'LineWidth',1.2);
 subplot(1,3,2);hold on;plot(tt,eq(j,4)*Y(:,2),'LineWidth',1.2);
 subplot(1,3,3);hold on;plot(tt,Y(:,3),'LineWidth',1.2);
end
subplot(1,3,1);title('Different intracellular DOX');ylabel('CF+CB (nM)');xlabel('Hours');legend(labels,'Location','best');grid on;
subplot(1,3,2);title('Same killing rate');ylabel('beta CB (h^{-1})');xlabel('Hours');grid on;
subplot(1,3,3);title('Same cell-count prediction');ylabel('N');xlabel('Hours');grid on;
saveas(f,fullfile(out,'04_equivalent_parameters.png'));close(f);
fprintf('[%s] Saved 04_equivalent_parameters.png (total %.1f s)\n',datestr(now,30),toc(runTimer));
if maxDelta>1e-3,warning('Equivalent-set numerical discrepancy: tighten tolerances.');end
fprintf('Equivalent-set max relative N difference (numerical): %.3g\n',maxDelta);
fprintf('These individual PK rates are representatives, NOT unique estimates.\n');
fprintf('Old simulated PK was NOT included in the loss. No held-out validation.\n');
results=struct('config',cfg,'growth',g,'summary',summary,'equivalent_sets',eq,...
 'combinations',abc,'control_runs',gRuns,'max_equivalent_N_difference',maxDelta,...
 'elapsed_seconds',toc(runTimer),'seed',42,'matlab_version',version,...
 'model_d_growth',gD,'model_d_fixed',fixedTable,'model_d_profile',profileD);
if doSensitivity
 results.model_d_sensitivity=sensitivityFixedPK(cs,cfg,gD,lo,hi,nStarts,out);
 results.elapsed_seconds=toc(runTimer);
end
if doDoseBeta
 results.model_d_dose_beta=doseBetaDiagnostic(cs,cfg,gD,lo(4),hi(4),nStarts,out);
 results.elapsed_seconds=toc(runTimer);
end
if doSaturation
 results.model_d_saturation=saturationDiagnostic(cs,cfg,gD,lo(4),hi(4),...
  nStarts,out,results.model_d_dose_beta);
 results.elapsed_seconds=toc(runTimer);
end
if ~isempty(delayMode)
 results.model_d_delay=delayStudy(cs,cfg,gD,lo,hi,nStarts,out,results.model_d_dose_beta,delayMode);
 results.elapsed_seconds=toc(runTimer);
end
save(fullfile(out,'joint_fit.mat'),'results','C','cs','runsA','runsB','runsC','runsD');
fprintf('[%s] COMPLETE in %.1f s. Results saved in %s\n',datestr(now,30),toc(runTimer),out);
end

%% Optimizer: positive parameters, log coordinates, explicit exploratory bounds.
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
function v=loss(p,cs,cfg,g)
v=0;
for j=1:numel(cs)
 Y=simulate(p,cs(j).t,cs(j),cfg,g);
 v=v+sum(((Y(:,3)-cs(j).y)/max(cs(j).y)).^2);
end
end
function Y=simulate(p,t,c,cfg,g)
t=t(:);assert(all(t>=0));Y=zeros(numel(t),3);
y0=[0;0;log(c.y(1))];stop=max(t);
if stop==0,Y(:,3)=c.y(1);return;end
opts=odeset('RelTol',cfg.rtol,'AbsTol',cfg.atol);
end1=min(c.wash,stop);
s1=ode15s(@(~,y) rhs(y,c.dose,p,cfg,g),[0,end1],y0,opts);
m=t<=end1;v=deval(s1,t(m));Y(m,:)=v';
if any(~m)
 s2=ode15s(@(~,y) rhs(y,0,p,cfg,g),[end1,stop],deval(s1,end1),opts);
 v=deval(s2,t(~m));Y(~m,:)=v';
end
Y(:,3)=exp(Y(:,3));
if any(~isfinite(Y(:))),error('Nonfinite simulation');end
end
function dy=rhs(y,ce,p,cfg,g)
% Physical states CF, CB and log(N). No damage state.
dy=[cfg.kEF*(cfg.VE/cfg.VI)*ce-(p(1)+p(2))*y(1)+p(3)*y(2);...
 p(2)*y(1)-p(3)*y(2);g(1)*(1-exp(y(3))/g(2))-p(4)*y(2)];
end
function y=logistic(g,t,n0)
y=g(2)./(1+(g(2)/n0-1).*exp(-g(1)*t));
end
function abc=combinations(p,cfg)
abc=[sum(p(1:3)),p(1)*p(3),p(4)*cfg.kEF*(cfg.VE/cfg.VI)*p(2)];
end
function eq=equivalentFamily(p,lo,hi)
a=sum(p(1:3));b=p(1)*p(3);gain=p(4)*p(2);
% Valid positive kBF lies between roots of x^2-a*x+b=0.
fast=(a+sqrt(max(0,a*a-4*b)))/2;slow=b/fast;
xs=logspace(log10(max(lo(3),slow*(1+1e-5))),...
 log10(min(hi(3),fast*(1-1e-5))),500)';
f=b./xs;fb=a-xs-f;beta=gain./fb;allp=[f,fb,xs,beta];
valid=all(isfinite(allp),2)&all(allp>=lo,2)&all(allp<=hi,2);
allp=allp(valid,:);eq=p;
if ~isempty(allp)
 ids=unique(round(linspace(1,size(allp,1),4)));eq=[p;allp(ids,:)];
end
end
function writeRuns(R,names,out,label)
names=[names,{'scaled_SSE','exitflag','function_evaluations'}];
writetable(array2table(R,'VariableNames',names),fullfile(out,[label,'_multistart.csv']));
end

function [C,cs]=embeddedData()
% Manual black-marker pixel centers from extracted PDF figures, not red fits.
% FigS3 original image: x0=210, x200=308; y0=558, y5000=504.
q=[210,519;259,262;271,267;284,219;295,170;307,189;318,170;330,171;341,197;354,179];
C=[(q(:,1)-210)*200/98,(558-q(:,2))*5000/54];
% Fig4 SUM-149PT: x0=172, x200~254.5; y0 differs by row; 10000=82.5px.
q1=[172,351;193,334;213,310;233,263;253,217;273,185;293,176;311,194;332,197;352,183];
q2=[172,733;193,732;213,725;233,716;253,699;273,602;293,568;304,550;321,549;343,571];
q3=[172,1123;193,1127;213,1129;233,1132;253,1131;273,1123;293,1092;304,1050;312,1015;371,936;383,957;392,928];
Q={q1,q2,q3};dose=[156.25,312.5,625];zero=[381,771,1160];
cs=struct([]);
for j=1:3
 q=Q{j};cs(j).dose=dose(j);cs(j).wash=6;
 cs(j).t=(q(:,1)-172)*200/82.5;cs(j).y=(zero(j)-q(:,2))*10000/82.5;
end
end

function profile=profileFixedPK(pD,sseD,cs,cfg,g,lo,hi,out)
% Conditional SSE profile: re-optimize beta at EVERY fixed kBF.
% kBF=0 is an explicit irreversible reference, outside the positive fit bounds.
profileDir=fullfile(out,'profile_checkpoints');
if ~exist(profileDir,'dir'),mkdir(profileDir);end
xs=unique([0;logspace(log10(lo(3)),log10(hi(3)),81)';...
 logspace(log10(max(lo(3),pD(3)/3)),log10(min(hi(3),pD(3)*3)),41)';pD(3)]);
rows=zeros(0,8);allCandidates=cell(numel(xs),1);
for j=1:numel(xs)
 timer=tic;
 fun=@(beta) loss([cfg.oldFE,cfg.oldFB,xs(j),beta],cs,cfg,g);
 [beta,sse,flag,evals,candidates]=profileBeta(fun,lo(4),hi(4),pD(4));
 nearBeta=log(beta/lo(4))<0.02 || log(hi(4)/beta)<0.02;
 rows(end+1,:)=[xs(j),beta,sse,flag,evals,nearBeta,toc(timer),j]; %#ok<AGROW>
 allCandidates{j}=candidates;
 point=rows(end,:);kBF=xs(j);
 save(fullfile(profileDir,sprintf('kbf_%03d.mat',j)),...
  'kBF','beta','sse','flag','evals','candidates','point');
 partial=array2table(rows,'VariableNames',{'kBF','beta','scaled_SSE',...
  'exitflag','function_evaluations','beta_near_bound','elapsed_seconds','point_id'});
 writetable(partial,fullfile(out,'model_d_kbf_profile.csv'));
 fprintf('[%s] D profile %d/%d: kBF=%.7g, beta=%.7g, SSE=%.9g, flag=%d, %.1f s\n',...
  datestr(now,30),j,numel(xs),kBF,beta,sse,flag,toc(timer));
end
profile=partial;profile.delta_SSE=profile.scaled_SSE-min(profile.scaled_SSE);
profile.relative_increase=profile.delta_SSE/min(profile.scaled_SSE);
writetable(profile,fullfile(out,'model_d_kbf_profile.csv'));
save(fullfile(out,'model_d_kbf_profile.mat'),'profile','allCandidates','pD','sseD','g','cfg');
if abs(min(profile.scaled_SSE)-sseD)>1e-4
 warning('PKPD:ProfileMismatch','Profile and joint D optimum differ; inspect optimization.');
end
f=figure('Visible','off','Color','w','Position',[80,80,1100,400]);
positive=profile.kBF>0;
subplot(1,2,1);semilogx(profile.kBF(positive),profile.scaled_SSE(positive),'o-');
hold on;plot(pD(3),sseD,'rp','MarkerSize',12);grid on;
xlabel('Fixed kBF (h^{-1})');ylabel('Profile scaled SSE');
title('D: beta re-optimized; fixed-reference conditional profile');
subplot(1,2,2);loglog(profile.kBF(positive),profile.beta(positive),'o-');grid on;
xlabel('Fixed kBF (h^{-1})');ylabel('Optimized beta');
title('No likelihood calibration or confidence interval');
saveas(f,fullfile(out,'05_model_d_kbf_profile.png'));close(f);
end

function [bestBeta,bestSSE,bestFlag,evals,candidates]=profileBeta(fun,lower,upper,guess)
% Scan log(beta), refine each detected local minimum, and check endpoints.
% The objective and ODE tolerances are unchanged; no Optimization Toolbox.
zs=unique([linspace(log(lower),log(upper),25),log(guess)]);
evals=0;vs=zeros(size(zs));
for k=1:numel(zs),vs(k)=objective(zs(k));end
[bestSSE,idx]=min(vs);bestBeta=exp(zs(idx));bestFlag=0;
if idx==1 || idx==numel(zs),bestFlag=2;end % Explicitly tested bound, not optimizer flag.
candidates=[zs(:),vs(:),zeros(numel(zs),1)];
local=find(vs(2:end-1)<=vs(1:end-2) & vs(2:end-1)<=vs(3:end) & ...
 (vs(2:end-1)<vs(1:end-2) | vs(2:end-1)<vs(3:end)))+1;
ids=unique([idx,local]);
opts=optimset('Display','off','MaxIter',650,'MaxFunEvals',1800,'TolX',1e-5);
for k=ids
 left=max(1,k-1);right=min(numel(zs),k+1);
 [z,v,flag]=fminbnd(@objective,zs(left),zs(right),opts);
 candidates(end+1,:)=[z,v,flag]; %#ok<AGROW>
 if v<bestSSE,bestSSE=v;bestBeta=exp(z);bestFlag=flag;end
end
if ~isfinite(bestSSE) || bestSSE>=1e9,error('PKPD:ProfileFit','No valid beta fit.');end
 function value=objective(z)
  evals=evals+1;
  value=boundedLoss(z,fun,lower,upper);
 end
end

function summary=sensitivityFixedPK(cs,cfg,g,lo,hi,nStarts,out)
% One-at-a-time exploratory reference perturbations, not uncertainty intervals.
% Each of the nine cases is re-fitted, including all three 1x reference cases.
originalRng=rng;rngCleanup=onCleanup(@() rng(originalRng)); %#ok<NASGU>
paramNames={'kEF','kFE','kFB'};multipliers=[0.5,1,2];
caseNames=cell(9,1);summaryRows=zeros(9,13);allRuns=[];predictions=[];residuals=[];byDose=[];
caseId=0;
for field=1:3
 for factor=multipliers
  caseId=caseId+1;caseNames{caseId}=paramNames{field};
  workCfg=cfg;fe=cfg.oldFE;fb=cfg.oldFB;
  if field==1,workCfg.kEF=cfg.kEF*factor;
  elseif field==2,fe=cfg.oldFE*factor;
  else,fb=cfg.oldFB*factor;end
  label=sprintf('D sensitivity %s x%g',paramNames{field},factor);
  checkpointDir=fullfile(out,'sensitivity_checkpoints',sprintf('case_%02d',caseId));
  if ~exist(checkpointDir,'dir'),mkdir(checkpointDir);end
  % Identical starting points for every scenario; same optimizer and bounds.
  rng(314159,'twister');
  fun=@(x) loss([fe,fb,x(1),x(2)],cs,workCfg,g);
  timer=tic;
  [x,sse,runs]=multifit(fun,[0.005,5e-4],lo([3,4]),hi([3,4]),...
   nStarts,label,checkpointDir);
  nearLower=log(x./lo([3,4]))<0.02;nearUpper=log(hi([3,4])./x)<0.02;
  summaryRows(caseId,:)=[caseId,factor,workCfg.kEF,fe,fb,x,sse,...
   nearLower(1),nearUpper(1),nearLower(2),nearUpper(2),toc(timer)];
  allRuns=[allRuns;[repmat(caseId,nStarts,1),(1:nStarts)',runs]]; %#ok<AGROW>
  p=[fe,fb,x];
  save(fullfile(checkpointDir,'case_fit.mat'),'p','sse','runs','workCfg','g',...
   'cs','lo','hi','factor','field','caseId');
  for doseId=1:numel(cs)
   c=cs(doseId);
   % Resolve the transient around the unchanged 6 h washout and the long tail.
   t=unique([linspace(0,12,241)';linspace(12,max(c.t),501)';c.wash;c.t]);
   Y=simulate(p,t,c,workCfg,g);effect=x(2)*Y(:,2);
   predictions=[predictions;[repmat([caseId,c.dose],numel(t),1),t,Y,effect]]; %#ok<AGROW>
   Yo=simulate(p,c.t,c,workCfg,g);rr=Yo(:,3)-c.y;
   residuals=[residuals;[repmat([caseId,c.dose],numel(c.t),1),c.t,c.y,Yo(:,3),rr,rr/max(c.y)]]; %#ok<AGROW>
   byDose=[byDose;[caseId,c.dose,sum((rr/max(c.y)).^2),sqrt(mean(rr.^2))]]; %#ok<AGROW>
  end
  summary=array2table(summaryRows(1:caseId,:),'VariableNames',{'case_id','multiplier',...
   'fixed_kEF','fixed_kFE','fixed_kFB','kBF','beta','scaled_SSE',...
   'kBF_near_lower','kBF_near_upper','beta_near_lower','beta_near_upper','elapsed_seconds'});
  summary.varied_parameter=caseNames(1:caseId);
  writetable(summary,fullfile(out,'model_d_sensitivity_summary.csv'));
  fprintf('[%s] Sensitivity case %d/9 COMPLETE: %s, SSE=%.10g\n',...
   datestr(now,30),caseId,label,sse);
 end
end
writetable(array2table(allRuns,'VariableNames',{'case_id','start','kBF','beta',...
 'scaled_SSE','exitflag','function_evaluations'}),fullfile(out,'model_d_sensitivity_starts.csv'));
writetable(array2table(predictions,'VariableNames',{'case_id','dose_nM','time_h','CF','CB','N','beta_CB'}),...
 fullfile(out,'model_d_sensitivity_predictions.csv'));
writetable(array2table(residuals,'VariableNames',{'case_id','dose_nM','time_h',...
 'observed_N','predicted_N','residual','scaled_residual'}),fullfile(out,'model_d_sensitivity_residuals.csv'));
writetable(array2table(byDose,'VariableNames',{'case_id','dose_nM','scaled_SSE','RMSE_cells'}),...
 fullfile(out,'model_d_sensitivity_by_dose.csv'));
save(fullfile(out,'model_d_sensitivity.mat'),'summary','allRuns','predictions',...
 'residuals','byDose','cfg','g','cs','lo','hi','nStarts');
end

function summary=doseBetaDiagnostic(cs,cfg,g,lower,upper,nStarts,out)
% Diagnostic only: freeze the archived reference D, fit one beta per dose.
% Exact doubles from run 43101764836 model_d_fit.mat, commit 950b8e5.
kBF=0.0049373889028291049;sharedBeta=0.00048479871747565857;
sourceRun='43101764836';sourceCommit='950b8e52466c0bc98cff98ce759eea70c5b1929b';
fixed=table(cfg.kEF,cfg.oldFE,cfg.oldFB,cfg.VI,cfg.VE,g(1),g(2),kBF,sharedBeta,...
 'VariableNames',{'kEF','kFE','kFB','VI','VE','kp','theta','kBF','shared_beta'});
writetable(fixed,fullfile(out,'model_d_dose_beta_fixed.csv'));
originalRng=rng;rngCleanup=onCleanup(@() rng(originalRng)); %#ok<NASGU>
rows=[];allRuns=[];predictions=[];residuals=[];
for j=1:numel(cs)
 c=cs(j);timer=tic;
 checkpointDir=fullfile(out,'dose_beta_checkpoints',sprintf('dose_%02d',j));
 if ~exist(checkpointDir,'dir'),mkdir(checkpointDir);end
 rng(161803,'twister');
 fun=@(beta) loss([cfg.oldFE,cfg.oldFB,kBF,beta],c,cfg,g);
 [beta,sse,runs]=multifit(fun,sharedBeta,lower,upper,nStarts,...
  sprintf('D dose beta %g nM',c.dose),checkpointDir);
 % Independent one-dimensional scan/refinement checks the multistart optimum.
 % It uses the existing profileBeta helper, unchanged bounds and ODE accuracy.
 [checkBeta,checkSSE,checkFlag,checkEvals,candidates]=profileBeta(fun,lower,upper,sharedBeta);
 % Audit the shallow local region without changing solver/optimizer accuracy.
 auditBetas=unique([linspace(max(lower,min(beta,checkBeta)*0.998),...
  min(upper,max(beta,checkBeta)*1.002),101),beta,checkBeta]);
 auditSSE=arrayfun(fun,auditBetas);
 localAudit=table(auditBetas',auditSSE','VariableNames',{'beta','scaled_SSE'});
 writetable(localAudit,fullfile(checkpointDir,'local_objective.csv'));
 [auditMin,auditIndex]=min(auditSSE);auditBeta=auditBetas(auditIndex);
 sharedSSE=fun(sharedBeta);
 nearLower=log(beta/lower)<0.02;nearUpper=log(upper/beta)<0.02;
 rows=[rows;[c.dose,sharedBeta,beta,sharedSSE,sse,nearLower,nearUpper,...
  checkBeta,checkSSE,checkFlag,checkEvals,auditBeta,auditMin,toc(timer)]]; %#ok<AGROW>
 allRuns=[allRuns;[repmat(c.dose,nStarts,1),(1:nStarts)',runs]]; %#ok<AGROW>
 save(fullfile(checkpointDir,'dose_fit.mat'),'beta','sse','runs','kBF','sharedBeta',...
  'sharedSSE','c','cfg','g','lower','upper','checkBeta','checkSSE','checkFlag',...
  'checkEvals','candidates','localAudit','sourceRun','sourceCommit');
 t=unique([linspace(0,12,241)';linspace(12,max(c.t),1001)';c.wash;c.t]);
 for modelId=1:2
  currentBeta=sharedBeta;if modelId==2,currentBeta=beta;end
  p=[cfg.oldFE,cfg.oldFB,kBF,currentBeta];
  Y=simulate(p,t,c,cfg,g);
  predictions=[predictions;[repmat([modelId,c.dose],numel(t),1),t,Y,...
   currentBeta*Y(:,2)]]; %#ok<AGROW>
  Yo=simulate(p,c.t,c,cfg,g);rr=Yo(:,3)-c.y;
  residuals=[residuals;[repmat([modelId,c.dose],numel(c.t),1),c.t,c.y,...
   Yo(:,3),rr,rr/max(c.y)]]; %#ok<AGROW>
 end
 summary=array2table(rows,'VariableNames',{'dose_nM','shared_beta','dose_beta',...
  'shared_scaled_SSE','dose_scaled_SSE','beta_near_lower','beta_near_upper',...
  'scan_beta','scan_scaled_SSE','scan_exitflag','scan_evaluations',...
  'local_audit_beta','local_audit_scaled_SSE','elapsed_seconds'});
 writetable(summary,fullfile(out,'model_d_dose_beta_summary.csv'));
 fprintf('[%s] D dose beta %g nM COMPLETE: beta=%.10g, shared SSE=%.10g, dose SSE=%.10g\n',...
  datestr(now,30),c.dose,beta,sharedSSE,sse);
end
writetable(array2table(allRuns,'VariableNames',{'dose_nM','start','beta',...
 'scaled_SSE','exitflag','function_evaluations'}),fullfile(out,'model_d_dose_beta_starts.csv'));
writetable(array2table(predictions,'VariableNames',{'model_id','dose_nM','time_h','CF','CB','N','beta_CB'}),...
 fullfile(out,'model_d_dose_beta_predictions.csv'));
writetable(array2table(residuals,'VariableNames',{'model_id','dose_nM','time_h',...
 'observed_N','predicted_N','residual','scaled_residual'}),fullfile(out,'model_d_dose_beta_residuals.csv'));
save(fullfile(out,'model_d_dose_beta.mat'),'summary','fixed','allRuns','predictions',...
 'residuals','cfg','g','cs','lower','upper','nStarts','kBF','sharedBeta','sourceRun','sourceCommit');
f=figure('Visible','off','Color','w','Position',[80,80,1200,720]);
colors=[0.1,0.6,0.3;0.75,0.2,0.15];styles={'--','-'};
for j=1:numel(cs)
 c=cs(j);
 subplot(2,3,j);plot(c.t,c.y,'ko','DisplayName','Digitized data');hold on;
 for modelId=1:2
  Q=predictions(predictions(:,1)==modelId & predictions(:,2)==c.dose,:);
  plot(Q(:,3),Q(:,6),styles{modelId},'Color',colors(modelId,:),'LineWidth',1.5);
 end
 title(sprintf('%g nM',c.dose));xlabel('Hours');ylabel('Cell count');grid on;
 legend({'Digitized data','Shared beta D','Dose-specific beta diagnostic'},'Location','best');
 subplot(2,3,j+3);hold on;
 for modelId=1:2
  Q=residuals(residuals(:,1)==modelId & residuals(:,2)==c.dose,:);
  plot(Q(:,3),Q(:,7),['o',styles{modelId}],'Color',colors(modelId,:),'LineWidth',1.2);
 end
 plot([0,max(c.t)],[0,0],'k:');xlabel('Hours');ylabel('Scaled residual: prediction - data');grid on;
end
saveas(f,fullfile(out,'06_model_d_dose_beta.png'));close(f);
end

function summary=saturationDiagnostic(cs,cfg,g,betaLower,betaUpper,nStarts,out,doseBeta)
% Shared saturating PD, with all archived D PK and growth values fixed.
kBF=0.0049373889028291049;sharedBeta=0.00048479871747565857;
lower=[betaLower,1e-3];upper=[betaUpper,1e8]; % beta: 1/(nM h), Ks: nM
originalRng=rng;rngCleanup=onCleanup(@() rng(originalRng)); %#ok<NASGU>
rng(271828,'twister');timer=tic;
fitDir=fullfile(out,'saturation_checkpoints');mkdir(fitDir);
fun=@(x) saturationLoss(x,cs,cfg,g,kBF);
[x,sse,runs]=multifit(fun,[sharedBeta,50],lower,upper,nStarts,...
 'D shared saturating PD',fitDir);
% A near-linear start checks nesting without treating large Ks as identified.
[extra,extraSSE,extraRuns]=multifit(fun,[sharedBeta,1e7],lower,upper,1,...
 'D saturation near linear',fitDir);
runs=[runs;extraRuns];
if extraSSE<sse,x=extra;sse=extraSSE;end
jointSeconds=toc(timer);
nearLower=log(x./lower)<0.02;nearUpper=log(upper./x)<0.02;
parameters=table(x(1),x(2),x(1)*x(2),sse,nearLower(1),nearUpper(1),...
 nearLower(2),nearUpper(2),jointSeconds,'VariableNames',{'beta','Ks_nM',...
 'max_killing_rate_per_h','scaled_SSE','beta_near_lower','beta_near_upper',...
 'Ks_near_lower','Ks_near_upper','joint_fit_seconds'});
writetable(parameters,fullfile(out,'model_d_saturation_parameters.csv'));
writetable(array2table(runs,'VariableNames',{'beta','Ks_nM','scaled_SSE',...
 'exitflag','function_evaluations'}),fullfile(out,'model_d_saturation_starts.csv'));
writetable(table({'beta';'Ks'},{'nM^-1 h^-1';'nM'},lower',upper',...
 'VariableNames',{'parameter','unit','lower','upper'}),fullfile(out,'model_d_saturation_bounds.csv'));
save(fullfile(out,'model_d_saturation_fit.mat'),'x','sse','runs','parameters',...
 'cfg','g','cs','kBF','sharedBeta','lower','upper','jointSeconds');

% Profile Ks over the entire exploratory range; re-optimize beta at each point.
ksValues=unique([logspace(log10(lower(2)),log10(upper(2)),81)';...
 logspace(log10(max(lower(2),x(2)/4)),log10(min(upper(2),x(2)*4)),41)';x(2);Inf]);
profileDir=fullfile(out,'saturation_profile_checkpoints');mkdir(profileDir);
profileRows=[];
for j=1:numel(ksValues)
 pointTimer=tic;Ks=ksValues(j);
 objective=@(beta) saturationLoss([beta,Ks],cs,cfg,g,kBF);
 [beta,value,flag,evals,candidates]=profileBeta(objective,betaLower,betaUpper,x(1));
 boundLow=log(beta/betaLower)<0.02;boundHigh=log(betaUpper/beta)<0.02;
 profileRows=[profileRows;[Ks,beta,value,flag,evals,boundLow,boundHigh,toc(pointTimer)]]; %#ok<AGROW>
 save(fullfile(profileDir,sprintf('ks_%03d.mat',j)),...
  'Ks','beta','value','flag','evals','candidates');
 profile=array2table(profileRows,'VariableNames',{'Ks_nM','beta','scaled_SSE',...
  'exitflag','function_evaluations','beta_near_lower','beta_near_upper','elapsed_seconds'});
 writetable(profile,fullfile(out,'model_d_saturation_profile.csv'));
 fprintf('[%s] D saturation profile %d/%d: Ks=%.7g, beta=%.7g, SSE=%.10g, %.1f s\n',...
  datestr(now,30),j,numel(ksValues),Ks,beta,value,toc(pointTimer));
end
profile.delta_SSE=profile.scaled_SSE-min(profile.scaled_SSE);
writetable(profile,fullfile(out,'model_d_saturation_profile.csv'));
save(fullfile(out,'model_d_saturation_profile.mat'),'profile','x','sse','cfg','g','kBF');

% Three comparators, always using the same data, scaling and fixed kBF.
predictions=[];residuals=[];byDose=[];limits=[];
for j=1:numel(cs)
 c=cs(j);t=unique([linspace(0,12,241)';linspace(12,max(c.t),1001)';c.wash;c.t]);
 pRef=[cfg.oldFE,cfg.oldFB,kBF,sharedBeta];
 reference=simulate(pRef,t,c,cfg,g);referenceEffect=sharedBeta*reference(:,2);
 for modelId=1:3
  beta=sharedBeta;Ks=Inf;
  if modelId==2,beta=doseBeta.dose_beta(j);end
  if modelId==3,beta=x(1);Ks=x(2);end
  p=[cfg.oldFE,cfg.oldFB,kBF,beta];
  if modelId==3
   Y=simulateSaturation(p,Ks,t,c,cfg,g);Yo=simulateSaturation(p,Ks,c.t,c,cfg,g);
  else
   Y=simulate(p,t,c,cfg,g);Yo=simulate(p,c.t,c,cfg,g);
  end
  effect=beta*Y(:,2)./(1+Y(:,2)/Ks);rr=Yo(:,3)-c.y;
  predictions=[predictions;[repmat([modelId,c.dose],numel(t),1),t,Y,effect]]; %#ok<AGROW>
  residuals=[residuals;[repmat([modelId,c.dose],numel(c.t),1),c.t,c.y,Yo(:,3),rr,rr/max(c.y)]]; %#ok<AGROW>
  byDose=[byDose;[modelId,c.dose,beta,Ks,sum((rr/max(c.y)).^2),sqrt(mean(rr.^2))]]; %#ok<AGROW>
 end
 % Hold beta at the original shared value to isolate the large-Ks limit.
 for Ks=[1e3,1e5,1e8,Inf]
  Y=simulateSaturation(pRef,Ks,t,c,cfg,g);
  effect=sharedBeta*Y(:,2)./(1+Y(:,2)/Ks);
  Yo=simulateSaturation(pRef,Ks,c.t,c,cfg,g);
  limits=[limits;[Ks,c.dose,sharedBeta,max(abs(Y(:,3)-reference(:,3)))/max(reference(:,3)),...
   max(abs(effect-referenceEffect))/max(referenceEffect),sum(((Yo(:,3)-c.y)/max(c.y)).^2)]]; %#ok<AGROW>
 end
end
modelId=(1:3)';modelNames={'Linear shared beta D';'Dose-specific beta diagnostic';'Shared saturating beta Ks'};
totalSSE=arrayfun(@(m) sum(byDose(byDose(:,1)==m,5)),modelId);
summary=table(modelId,modelNames,totalSSE,[1;3;2],...
 'VariableNames',{'model_id','model','scaled_SSE','conditional_free_PD_parameters'});
writetable(summary,fullfile(out,'model_d_saturation_summary.csv'));
writetable(array2table(byDose,'VariableNames',{'model_id','dose_nM','beta','Ks_nM','scaled_SSE','RMSE_cells'}),...
 fullfile(out,'model_d_saturation_by_dose.csv'));
writetable(array2table(predictions,'VariableNames',{'model_id','dose_nM','time_h','CF','CB','N','killing_rate'}),...
 fullfile(out,'model_d_saturation_predictions.csv'));
writetable(array2table(residuals,'VariableNames',{'model_id','dose_nM','time_h','observed_N',...
 'predicted_N','residual','scaled_residual'}),fullfile(out,'model_d_saturation_residuals.csv'));
writetable(array2table(limits,'VariableNames',{'Ks_nM','dose_nM','fixed_beta',...
 'max_N_relative_difference','max_killing_relative_difference','scaled_SSE'}),...
 fullfile(out,'model_d_saturation_linear_limit.csv'));
save(fullfile(out,'model_d_saturation.mat'),'summary','parameters','profile','runs',...
 'byDose','predictions','residuals','limits','cfg','g','cs','kBF','sharedBeta','lower','upper');
f=figure('Visible','off','Color','w','Position',[80,80,1200,720]);
colors=[0.1,0.6,0.3;0.75,0.2,0.15;0.1,0.3,0.85];styles={'--',':','-'};
for j=1:numel(cs)
 c=cs(j);subplot(2,3,j);plot(c.t,c.y,'ko');hold on;
 for m=1:3
  Q=predictions(predictions(:,1)==m & predictions(:,2)==c.dose,:);
  plot(Q(:,3),Q(:,6),styles{m},'Color',colors(m,:),'LineWidth',1.5);
 end
 title(sprintf('%g nM',c.dose));xlabel('Hours');ylabel('Cell count');grid on;
 legend([{'Digitized data'};modelNames],'Location','best');
 subplot(2,3,j+3);hold on;
 for m=1:3
  Q=residuals(residuals(:,1)==m & residuals(:,2)==c.dose,:);
  plot(Q(:,3),Q(:,7),['o',styles{m}],'Color',colors(m,:));
 end
 plot([0,max(c.t)],[0,0],'k:');xlabel('Hours');ylabel('Scaled residual: prediction - data');grid on;
end
saveas(f,fullfile(out,'07_model_d_saturation.png'));close(f);
f=figure('Visible','off','Color','w','Position',[80,80,1100,400]);finite=isfinite(profile.Ks_nM);
subplot(1,2,1);semilogx(profile.Ks_nM(finite),profile.scaled_SSE(finite),'o-');hold on;
yline(profile.scaled_SSE(~finite),'k--','Linear limit');plot(x(2),sse,'rp','MarkerSize',12);
xlabel('Fixed Ks (nM)');ylabel('Profile scaled SSE');grid on;title('Beta re-optimized at each Ks');
subplot(1,2,2);loglog(profile.Ks_nM(finite),profile.beta(finite),'o-');
xlabel('Fixed Ks (nM)');ylabel('Optimized beta (nM^{-1} h^{-1})');grid on;title('Conditional profile; not a confidence interval');
saveas(f,fullfile(out,'08_model_d_saturation_profile.png'));close(f);
fprintf('[%s] D saturation COMPLETE: beta=%.10g, Ks=%.10g, SSE=%.10g, %.1f s\n',...
 datestr(now,30),x(1),x(2),sse,toc(timer));
end

function v=saturationLoss(x,cs,cfg,g,kBF)
v=0;p=[cfg.oldFE,cfg.oldFB,kBF,x(1)];
for j=1:numel(cs)
 Y=simulateSaturation(p,x(2),cs(j).t,cs(j),cfg,g);
 v=v+sum(((Y(:,3)-cs(j).y)/max(cs(j).y)).^2);
end
end

function Y=simulateSaturation(p,Ks,t,c,cfg,g)
% Same three states, integrator, tolerances, initial conditions and washout.
t=t(:);assert(all(t>=0) && Ks>0);Y=zeros(numel(t),3);
y0=[0;0;log(c.y(1))];stop=max(t);
if stop==0,Y(:,3)=c.y(1);return;end
opts=odeset('RelTol',cfg.rtol,'AbsTol',cfg.atol);
end1=min(c.wash,stop);
s1=ode15s(@(~,y) rhsSaturation(y,c.dose,p,Ks,cfg,g),[0,end1],y0,opts);
m=t<=end1;v=deval(s1,t(m));Y(m,:)=v';
if any(~m)
 s2=ode15s(@(~,y) rhsSaturation(y,0,p,Ks,cfg,g),[end1,stop],deval(s1,end1),opts);
 v=deval(s2,t(~m));Y(~m,:)=v';
end
Y(:,3)=exp(Y(:,3));
if any(~isfinite(Y(:))),error('Nonfinite saturation simulation');end
end

function dy=rhsSaturation(y,ce,p,Ks,cfg,g)
dy=rhs(y,ce,p,cfg,g); % The two PK derivatives are exactly the original ones.
dy(3)=g(1)*(1-exp(y(3))/g(2))-p(4)*y(2)/(1+y(2)/Ks);
end

function summary=delayStudy(cs,cfg,g,oldLower,oldUpper,nStarts,out,doseBeta,mode)
% Two serial effect-transit states; no molecular interpretation is imposed.
kBFRef=0.0049373889028291049;betaRef=0.00048479871747565857;
lower=[oldLower(3),oldLower(4),1e-4];upper=[oldUpper(3),oldUpper(4),100];
root=fileparts(mfilename('fullpath'));
paper=jsondecode(fileread(fullfile(root,'reports','paper_time_definitions.json')));
if paper.verified
 assert(isnumeric(paper.origin_in_current_data_hours) && isscalar(paper.origin_in_current_data_hours) && ...
  isfinite(paper.origin_in_current_data_hours),'PKPD:PaperOrigin','Verified paper origin is missing.');
end
fid=fopen(fullfile(out,'paper_time_definitions.json'),'w');fprintf(fid,'%s\n',jsonencode(paper));fclose(fid);
originalRng=rng;cleanup=onCleanup(@() rng(originalRng)); %#ok<NASGU>
checks=delayNumericalChecks(cs,cfg,g,kBFRef,betaRef,out);
fitDir=fullfile(out,'delay_checkpoints');mkdir(fitDir);
rng(57721,'twister');
f1=@(x) delayLoss([kBFRef,x],cs,cfg,g);
[x1,sse1,runs1]=multifit(f1,[betaRef,0.02],lower(2:3),upper(2:3),nStarts,'Delay stage 1',fitDir);
[xe,se,re]=multifit(f1,[betaRef,upper(3)],lower(2:3),upper(2:3),1,'Delay stage 1 fast',fitDir);
runs1=[runs1;re];if se<sse1,x1=xe;sse1=se;end
p1=[kBFRef,x1];
save(fullfile(out,'delay_stage1.mat'),'p1','sse1','runs1','cfg','g','cs','lower','upper');
rng(57722,'twister');
f2=@(x) delayLoss(x,cs,cfg,g);
[p2,sse2,runs2]=multifit(f2,p1,lower,upper,nStarts,'Delay stage 2',fitDir);
[xe,se,re]=multifit(f2,[kBFRef,betaRef,upper(3)],lower,upper,1,'Delay stage 2 fast',fitDir);
runs2=[runs2;re];if se<sse2,p2=xe;sse2=se;end
profile=table();
if strcmp(mode,'full')
 profileDir=fullfile(out,'delay_r_profile_checkpoints');mkdir(profileDir);
 rs=unique([logspace(log10(lower(3)),log10(upper(3)),21)';...
  logspace(log10(max(lower(3),p2(3)/4)),log10(min(upper(3),p2(3)*4)),11)';p2(3)]);
 profileRows=[];
 for j=1:numel(rs)
  r=rs(j);rng(57723,'twister');pointTimer=tic;
  fun=@(x) delayLoss([x,r],cs,cfg,g);
  folder=fullfile(profileDir,sprintf('point_%03d',j));mkdir(folder);
  [xp,sp,rp]=multifit(fun,p2(1:2),lower(1:2),upper(1:2),3,...
   sprintf('Delay r profile %03d',j),folder);
  near=log(xp./lower(1:2))<0.02 | log(upper(1:2)./xp)<0.02;
  profileRows=[profileRows;[r,xp,sp,near,toc(pointTimer)]]; %#ok<AGROW>
  profile=array2table(profileRows,'VariableNames',{'r_per_h','kBF','alpha','scaled_SSE',...
   'kBF_near_bound','alpha_near_bound','elapsed_seconds'});
  writetable(profile,fullfile(out,'delay_stage2_r_profile.csv'));
  save(fullfile(folder,'point.mat'),'r','xp','sp','rp');
 end
 [sp,idx]=min(profile.scaled_SSE);
 if sp<sse2-1e-7
  guess=[profile.kBF(idx),profile.alpha(idx),profile.r_per_h(idx)];
  [xe,se,re]=multifit(f2,guess,lower,upper,1,'Delay stage 2 profile polish',fitDir);
  runs2=[runs2;re];if se<sse2,p2=xe;sse2=se;end
 end
end
save(fullfile(out,'delay_stage2.mat'),'p2','sse2','runs2','profile','cfg','g','cs','lower','upper');
writetable(array2table(runs1,'VariableNames',{'alpha','r_per_h','scaled_SSE','exitflag','function_evaluations'}),...
 fullfile(out,'delay_stage1_starts.csv'));
writetable(array2table(runs2,'VariableNames',{'kBF','alpha','r_per_h','scaled_SSE','exitflag','function_evaluations'}),...
 fullfile(out,'delay_stage2_starts.csv'));
P=[p1;p2];nearLow=log(P./lower)<0.02;nearHigh=log(upper./P)<0.02;
params=array2table([(1:2)',P,[sse1;sse2],nearLow,nearHigh],...
 'VariableNames',{'stage','kBF','alpha','r_per_h','scaled_SSE','kBF_near_lower','alpha_near_lower',...
 'r_near_lower','kBF_near_upper','alpha_near_upper','r_near_upper'});
writetable(params,fullfile(out,'delay_parameters.csv'));
writetable(table({'kBF';'alpha';'r'},{'h^-1';'nM^-1 h^-1';'h^-1'},lower',upper',...
 'VariableNames',{'parameter','unit','lower','upper'}),fullfile(out,'delay_bounds.csv'));

% Time-form diagnostics are explicitly gated on verified primary definitions.
timeA=nan(numel(cs),1);timeB=nan(numel(cs),2);timeRuns=[];
if paper.verified
 for j=1:numel(cs)
  c=cs(j);rng(57724,'twister');
  for kind=1:2
   guess=0.01;loT=1e-9;hiT=1;
   if kind==2,guess=[0.02,0.01];loT=[1e-9,lower(3)];hiT=[1,upper(3)];end
   folder=fullfile(out,'time_form_checkpoints',sprintf('form_%d_dose_%d',kind,j));mkdir(folder);
   objective=@(p) timeFormLoss(p,kind,c,g,paper.origin_in_current_data_hours,cfg);
   [p,sse,rr]=multifit(objective,guess,loT,hiT,nStarts,...
    sprintf('Time form %d dose %g',kind,c.dose),folder);
   if kind==1,timeA(j)=p;rr=[rr(:,1),nan(nStarts,1),rr(:,2:end)];else,timeB(j,:)=p;end
   timeRuns=[timeRuns;[repmat([kind,c.dose],nStarts,1),(1:nStarts)',rr]]; %#ok<AGROW>
   save(fullfile(folder,'fit.mat'),'p','sse','rr','c','paper','g','cfg','loT','hiT');
  end
 end
 writetable(array2table(timeRuns,'VariableNames',{'form','dose_nM','start','kd_per_h','r_per_h',...
  'scaled_SSE','exitflag','function_evaluations'}),fullfile(out,'time_form_starts.csv'));
else
 fprintf('TIME-FORM A/B NOT RUN: primary-paper time definitions are unverified.\n');
end

predictions=[];residuals=[];byDose=[];
active=[1,2,5,6];if paper.verified,active=1:6;end
for j=1:numel(cs)
 c=cs(j);t=unique([linspace(0,12,241)';linspace(12,max(c.t),1001)';c.wash;c.t]);
 for model=active
  kbf=kBFRef;amplitude=betaRef;r=Inf;
  if model<=2
   if model==2,amplitude=doseBeta.dose_beta(j);end
   Y3=simulate([cfg.oldFE,cfg.oldFB,kbf,amplitude],t,c,cfg,g);
   Yo3=simulate([cfg.oldFE,cfg.oldFB,kbf,amplitude],c.t,c,cfg,g);
   Y=[Y3(:,1:2),nan(numel(t),2),Y3(:,3)];No=Yo3(:,3);effect=amplitude*Y3(:,2);
  elseif model<=4
   kind=model-2;kbf=NaN;
   if kind==1,p=timeA(j);else,p=timeB(j,:);r=p(2);end
   amplitude=p(1);origin=paper.origin_in_current_data_hours;
   [N,effect]=simulateTimeForm(p,kind,t,c,g,origin,cfg);
   No=simulateTimeForm(p,kind,c.t,c,g,origin,cfg);Y=[nan(numel(t),4),N];
  else
   p=P(model-4,:);kbf=p(1);amplitude=p(2);r=p(3);
   Y=simulateDelay(p,t,c,cfg,g);Yo=simulateDelay(p,c.t,c,cfg,g);No=Yo(:,5);effect=p(2)*Y(:,4);
  end
  rr=No-c.y;
  predictions=[predictions;[repmat([model,c.dose],numel(t),1),t,Y,effect]]; %#ok<AGROW>
  residuals=[residuals;[repmat([model,c.dose],numel(c.t),1),c.t,c.y,No,rr,rr/max(c.y)]]; %#ok<AGROW>
  byDose=[byDose;[model,c.dose,kbf,amplitude,r,sum((rr/max(c.y)).^2),sqrt(mean(rr.^2))]]; %#ok<AGROW>
 end
end
names={'Archived linear D';'Dose-specific beta diagnostic';'Paper-inspired constant-rate diagnostic';...
 'Paper-inspired transient-rate diagnostic';'CB delay stage 1 fixed kBF';'CB delay stage 2 shared kBF'};
scores=nan(6,1);states=repmat({'completed'},6,1);
for m=active,scores(m)=sum(byDose(byDose(:,1)==m,6));end
if ~paper.verified,states(3:4)={'not_run_unverified_paper';'not_run_unverified_paper'};end
summary=table((1:6)',names,states,scores,[0;3;3;6;2;3],...
 'VariableNames',{'model_id','model','status','scaled_SSE','parameters_refitted_this_comparison'});
writetable(summary,fullfile(out,'delay_summary.csv'));
writetable(array2table(byDose,'VariableNames',{'model_id','dose_nM','kBF','amplitude','r_per_h',...
 'scaled_SSE','RMSE_cells'}),fullfile(out,'delay_by_dose.csv'));
writetable(array2table(predictions,'VariableNames',{'model_id','dose_nM','time_h','CF','CB','E1','E2','N','killing_rate'}),...
 fullfile(out,'delay_predictions.csv'));
writetable(array2table(residuals,'VariableNames',{'model_id','dose_nM','time_h','observed_N',...
 'predicted_N','residual','scaled_residual'}),fullfile(out,'delay_residuals.csv'));
save(fullfile(out,'delay_study.mat'),'summary','params','predictions','residuals','byDose',...
 'p1','p2','sse1','sse2','runs1','runs2','profile','checks','paper','timeA','timeB','timeRuns',...
 'cfg','g','cs','kBFRef','betaRef','lower','upper','mode');
end

function v=delayLoss(p,cs,cfg,g)
v=0;
for j=1:numel(cs)
 Y=simulateDelay(p,cs(j).t,cs(j),cfg,g);
 v=v+sum(((Y(:,5)-cs(j).y)/max(cs(j).y)).^2);
end
end

function [Y,transition]=simulateDelay(p,t,c,cfg,g)
t=t(:);assert(all(t>=0));Y=zeros(numel(t),5);transition=struct();
y0=[0;0;0;0;log(c.y(1))];stop=max(t);
if stop==0,Y(:,5)=c.y(1);return;end
opts=odeset('RelTol',cfg.rtol,'AbsTol',cfg.atol);end1=min(c.wash,stop);
s1=ode15s(@(~,y) rhsDelay(y,c.dose,p,cfg,g),[0,end1],y0,opts);
m=t<=end1;v=deval(s1,t(m));Y(m,:)=v';
if any(~m)
 before=deval(s1,end1);
 s2=ode15s(@(~,y) rhsDelay(y,0,p,cfg,g),[end1,stop],before,opts);
 v=deval(s2,t(~m));Y(~m,:)=v';
 transition.left=before;transition.right=s2.y(:,1);
 transition.derivative_jump=rhsDelay(before,0,p,cfg,g)-rhsDelay(before,c.dose,p,cfg,g);
end
Y(:,5)=exp(Y(:,5));
if any(~isfinite(Y(:))),error('Nonfinite delay simulation');end
end

function dy=rhsDelay(y,ce,p,cfg,g)
% Physical states CF, CB, E1, E2, log(N); E1/E2 are effect-transit states.
original=rhs([y(1:2);y(5)],ce,[cfg.oldFE,cfg.oldFB,p(1),p(2)],cfg,g);
dy=[original(1:2);p(3)*(y(2)-y(3));p(3)*(y(3)-y(4));...
 g(1)*(1-exp(y(5))/g(2))-p(2)*y(4)];
end

function checks=delayNumericalChecks(cs,cfg,g,kBF,beta,out)
rows=[];washRows=[];
for j=1:numel(cs)
 c=cs(j);t=unique([linspace(0,12,241)';linspace(12,max(c.t),501)';c.t]);
 Y=simulateDelay([kBF,0,0.01],t,c,cfg,g);expected=logistic(g,t,c.y(1));
 logisticError=max(abs(Y(:,5)-expected))/max(expected);
 assert(logisticError<1e-5,'PKPD:DelayLogistic','alpha=0 does not recover logistic growth.');
 ref=simulate([cfg.oldFE,cfg.oldFB,kBF,beta],t,c,cfg,g);
 for r=[1,100,1e4]
  [Y,tr]=simulateDelay([kBF,beta,r],t,c,cfg,g);
  nError=max(abs(Y(:,5)-ref(:,3)))/max(ref(:,3));
  eError=max(abs(Y(:,4)-ref(:,2)))/max(ref(:,2));
  assert(all(Y(:)>=0),'PKPD:DelayNonnegative','Negative state in numerical check.');
  jumpError=max(abs(tr.left-tr.right));
  derivativeExpected=[-cfg.kEF*(cfg.VE/cfg.VI)*c.dose;0;0;0;0];
  derivativeError=max(abs(tr.derivative_jump-derivativeExpected));
  assert(jumpError==0 && derivativeError<1e-10,'PKPD:DelayWashout','Washout reset or RHS discontinuity is incorrect.');
  if r==1e4
   assert(nError<1e-4 && eError<1e-4,'PKPD:DelayFastLimit','Large-r limit disagrees with instantaneous effect.');
  end
  rows=[rows;[c.dose,r,logisticError,nError,eError,jumpError,derivativeError,min(Y(:))]]; %#ok<AGROW>
 end
 % Independent matrix-exponential check of the linear PK + transit subsystem.
 p=[kBF,beta,0.01];tt=unique([0;0.05;1;c.wash;6.1;24;100;max(c.t)]);
 Y=simulateDelay(p,tt,c,cfg,g);A=[-(cfg.oldFE+cfg.oldFB),kBF,0,0;...
  cfg.oldFB,-kBF,0,0;0,p(3),-p(3),0;0,0,p(3),-p(3)];
 b=[cfg.kEF*(cfg.VE/cfg.VI)*c.dose;0;0;0];M=[A,b;zeros(1,5)];
 atWash=expm(M*c.wash)*[zeros(4,1);1];exact=zeros(numel(tt),4);
 for k=1:numel(tt)
  if tt(k)<=c.wash,z=expm(M*tt(k))*[zeros(4,1);1];z=z(1:4);
  else,z=expm(A*(tt(k)-c.wash))*atWash(1:4);end
  exact(k,:)=z';
 end
 exactError=max(abs(Y(:,1:4)-exact),[],1)./max(exact,[],1);
 assert(all(exactError<1e-5),'PKPD:DelayLinearStates','PK/transit ODEs disagree with matrix exponential.');
 washRows=[washRows;[c.dose,exactError]]; %#ok<AGROW>
end
checks=struct('limits',rows,'matrix_exponential',washRows);
writetable(array2table(rows,'VariableNames',{'dose_nM','r_per_h','logistic_relative_error',...
 'N_vs_instant_relative_error','E2_vs_CB_relative_error','washout_state_jump',...
 'washout_derivative_error','minimum_state'}),fullfile(out,'delay_numerical_checks.csv'));
writetable(array2table(washRows,'VariableNames',{'dose_nM','CF_relative_error','CB_relative_error',...
 'E1_relative_error','E2_relative_error'}),fullfile(out,'delay_matrix_exponential_checks.csv'));
save(fullfile(out,'delay_numerical_checks.mat'),'checks');
fprintf('DELAY NUMERICAL CHECKS PASSED: logistic, fast limit, continuous washout, nonnegative states, matrix exponential.\n');
end

function v=timeFormLoss(p,kind,c,g,origin,cfg)
N=simulateTimeForm(p,kind,c.t,c,g,origin,cfg);
v=sum(((N-c.y)/max(c.y)).^2);
end

function [N,effect]=simulateTimeForm(p,kind,t,c,g,origin,cfg)
t=t(:);N=zeros(size(t));effect=arrayfun(@(u) timeFormRate(u,p,kind,origin),t);
breaks=unique([0;max(0,min(max(t),origin));max(t)]);z=log(c.y(1));
N(t==0)=c.y(1);opts=odeset('RelTol',cfg.rtol,'AbsTol',cfg.atol);
for k=1:numel(breaks)-1
 left=breaks(k);right=breaks(k+1);
 % Split at verified effect onset; before it all killing is exactly zero.
 if right<=origin,rate=@(~) 0;else,rate=@(u) timeFormRate(u,p,kind,origin);end
 sol=ode15s(@(u,y) g(1)*(1-exp(y)/g(2))-rate(u),[left,right],z,opts);
 mask=t>left & t<=right;N(mask)=exp(deval(sol,t(mask)))';z=deval(sol,right);
end
if any(~isfinite(N)),error('Nonfinite time-form solution');end
end

function value=timeFormRate(t,p,kind,origin)
tau=t-origin;
if tau<0,value=0;elseif kind==1,value=p(1);else,value=p(1)*p(2)*tau*exp(1-p(2)*tau);end
end
