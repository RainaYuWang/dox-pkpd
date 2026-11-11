function results = run_joint_pkpd(nStarts,doSensitivity)
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
