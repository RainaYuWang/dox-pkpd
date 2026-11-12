function results=run_saturation_link(mode)
% Shared instantaneous nonlinear PK-PD; no damage states or changed growth.
% Full-run release: identical numerical code to successful quick 43160254529.
mode=validatestring(mode,{'sat-quick','sat-full'});
root=fileparts(mfilename('fullpath'));out=fullfile(root,'joint_results');
assert(~exist(out,'dir'),'Move previous joint_results before running.');mkdir(out);
timer=tic;status=struct('status','running','mode',mode,'old_models_refitted',false,...
 'matlab_version',version,'commit',getenv('GITHUB_SHA'),'run_id',getenv('GITHUB_RUN_ID'));
writeStatus(out,status);diary(fullfile(out,'run.log'));cleanup=onCleanup(@() diary('off'));
try
 baseDir=fullfile(root,'baseline-results','f05543e978e425785fdf918751185cd1a733d300','43137925939-1','joint_results');
 base=load(fullfile(baseDir,'delay_study.mat'),'cfg','g','cs');cfg=base.cfg;g=base.g;cs=base.cs;
 assert(isequal(g,[0.031094959444983555,35143.615284616222]));
 assert(cfg.rtol==1e-7 && cfg.atol==1e-9);
 copyfile(fullfile(fileparts(baseDir),'manifest.json'),fullfile(out,'baseline_manifest.json'));
 kref=0.0049373889028291049;bref=0.00048479871747565857;
 lo=[1e-7,1e-9,1e-3];hi=[10,0.1,1e8]; % kBF, beta, Ks: retain existing bounds.
 n=2;if strcmp(mode,'sat-full'),n=8;end
 save(fullfile(out,'saturation_inputs.mat'),'cfg','g','cs','lo','hi','kref','bref','mode');
 writetable(table({'kBF';'beta';'Ks'},{'h^-1';'nM^-1 h^-1';'nM'},lo',hi',...
  'VariableNames',{'parameter','unit','lower','upper'}),fullfile(out,'saturation_bounds.csv'));
 checks=[];
 for j=1:numel(cs)
  c=cs(j);t=unique([linspace(0,max(c.t),1001)';c.t;c.wash]);p=[cfg.oldFE,cfg.oldFB,kref,bref];
  lin=simulate(p,t,c,cfg,g);exact=simulateSaturation(p,Inf,t,c,cfg,g);
  large=simulateSaturation(p,1e10,t,c,cfg,g);
  zero=simulateSaturation([p(1:3),0],50,t,c,cfg,g);
  lg=logistic(g,t,c.y(1));
  e0=max(abs(zero(:,3)-lg))/max(lg);eInf=max(abs(exact(:,3)-lin(:,3)))/max(lin(:,3));
  eLarge=max(abs(large(:,3)-lin(:,3)))/max(lin(:,3));
  % Independent affine matrix exponential for the PK and washout boundary.
  A=[-(cfg.oldFE+cfg.oldFB),kref;cfg.oldFB,-kref];u=[cfg.kEF*cfg.VE/cfg.VI*c.dose;0];
  M=[A,u;0,0,0];yw=expm(M*c.wash)*[0;0;1];truth=zeros(numel(t),2);
  for i=1:numel(t)
   if t(i)<=c.wash,v=expm(M*t(i))*[0;0;1];else,v=[expm(A*(t(i)-c.wash))*yw(1:2);1];end
   truth(i,:)=v(1:2)';
  end
  pkerr=max(abs(exact(:,1:2)-truth),[],'all')/max(truth,[],'all');
  assert(e0<1e-4 && eInf<1e-8 && eLarge<1e-4 && pkerr<1e-5);
  assert(all(exact(:,3)>0) && min(exact(:,1:2),[],'all')>=-1e-8);
  checks=[checks;c.dose,e0,eInf,eLarge,pkerr];
 end
 writetable(array2table(checks,'VariableNames',{'dose_nM','zero_logistic_error','infinite_Ks_error','large_Ks_error','matrix_PK_error'}),fullfile(out,'saturation_checks.csv'));
 params=[];allstarts=[];profiles=[];pred=[];res=[];
 for stage=1:2
  folder=fullfile(out,sprintf('saturation_stage%d_checkpoints',stage));mkdir(folder);rng(271828+stage,'twister');
  if stage==1
   fn=@(x) saturationLoss(x,cs,cfg,g,kref);guess=[bref,50];lower=lo(2:3);upper=hi(2:3);near=[bref,1e7];
  else
   fn=@(x) saturationLoss(x(2:3),cs,cfg,g,x(1));guess=best1;lower=lo;upper=hi;near=[kref,bref,1e7];
  end
  [x,sse,runs]=multifit(fn,guess,lower,upper,n,sprintf('Saturation stage %d',stage),folder);
  [xe,se,re]=multifit(fn,near,lower,upper,1,sprintf('Stage %d linear start',stage),folder);runs=[runs;re];
  if se<sse,x=xe;sse=se;end
  if stage==2 && fn(best1)<sse,x=best1;sse=fn(best1);end
  if stage==1,p=[kref,x];rr=[repmat(kref,size(runs,1),1),runs];else,p=x;rr=runs;end
  % Full profiles re-optimize remaining parameters. Include current optimum.
  if strcmp(mode,'sat-full')
   grid=unique([logspace(-3,8,25),p(3)*[0.5,0.8,1,1.25,2],Inf]);grid=grid(grid>=lo(3) & (grid<=hi(3)|isinf(grid)));
   for z=1:numel(grid)
    ks=grid(z);
    if stage==1
     fun=@(b) saturationLoss([b,ks],cs,cfg,g,kref);
     [b,v,flag,ev]=profileBeta(fun,lo(2),hi(2),p(2));q=[kref,b,ks];
    else
     fun=@(a) saturationLoss([a(2),ks],cs,cfg,g,a(1));
     [a,v,rt]=multifit(fun,p(1:2),lo(1:2),hi(1:2),2,sprintf('Stage2 profile %02d',z),folder);
     [~,ix]=min(rt(:,3));flag=rt(ix,4);ev=sum(rt(:,5));q=[a,ks];
    end
    profiles=[profiles;stage,q,v,flag,ev];
    writetable(array2table(profiles,'VariableNames',{'stage','kBF','beta','Ks_nM','scaled_SSE','exitflag','evaluations'}),fullfile(out,'saturation_profile.csv'));
    if isfinite(ks) && v<sse,p=q;sse=v;end
   end
   % Refine if profile discovers a better region.
   if stage==1,seed=p(2:3);else,seed=p;end
   [xf,sf,rf]=multifit(fn,seed,lower,upper,1,sprintf('Stage %d profile refine',stage),folder);
   if stage==1,rr=[rr;kref,rf];else,rr=[rr;rf];end
   if sf<sse,sse=sf;if stage==1,p=[kref,xf];else,p=xf;end;end
  end
  if stage==1,best1=p;end
  params=[params;stage,p,p(2)*p(3),sse];allstarts=[allstarts;repmat(stage,size(rr,1),1),(1:size(rr,1))',rr];
  for j=1:numel(cs)
   c=cs(j);t=unique([linspace(0,12,241)';linspace(12,max(c.t),1501)';c.t;c.wash]);
   Y=simulateSaturation([cfg.oldFE,cfg.oldFB,p(1:2)],p(3),t,c,cfg,g);
   Yo=simulateSaturation([cfg.oldFE,cfg.oldFB,p(1:2)],p(3),c.t,c,cfg,g);
   kill=p(2)*Y(:,2)./(1+Y(:,2)/p(3));r=Yo(:,3)-c.y;
   assert(all(Y(:,3)>0) && all(kill>=-1e-10));
   pred=[pred;repmat([stage,c.dose],numel(t),1),t,Y,kill];
   res=[res;repmat([stage,c.dose],numel(c.t),1),c.t,c.y,Yo(:,3),r,r/max(c.y)];
  end
  assert(abs(sum(res(res(:,1)==stage,7).^2)-sse)<1e-8);
  save(fullfile(out,sprintf('saturation_stage%d.mat',stage)),'p','sse','rr','cfg','g','cs');
 end
 parameters=array2table(params,'VariableNames',{'stage','kBF','beta','Ks_nM','kmax_per_h','scaled_SSE'});
 starts=array2table(allstarts,'VariableNames',{'stage','start','kBF','beta','Ks_nM','scaled_SSE','exitflag','evaluations'});
 writetable(parameters,fullfile(out,'saturation_parameters.csv'));writetable(starts,fullfile(out,'saturation_starts.csv'));
 writetable(array2table(pred,'VariableNames',{'stage','dose_nM','time_h','CF','CB','N','killing_rate'}),fullfile(out,'saturation_predictions.csv'));
 writetable(array2table(res,'VariableNames',{'stage','dose_nM','time_h','observed_N','predicted_N','residual','scaled_residual'}),fullfile(out,'saturation_residuals.csv'));
 for file={'delay_summary.csv','delay_predictions.csv','delay_residuals.csv'}
  copyfile(fullfile(baseDir,file{1}),fullfile(out,['baseline_',file{1}]));
 end
 results=struct('parameters',parameters,'starts',starts,'checks',checks);
 save(fullfile(out,'saturation_results.mat'),'results','params','pred','res','profiles','cfg','g','cs');
 status.status='success';status.elapsed_seconds=toc(timer);status.scaled_SSE=params(:,6)';status.numerical_checks='passed';writeStatus(out,status);
 disp(parameters);fprintf('SATURATION CI PASSED %.2f seconds\n',toc(timer));
catch err
 status.status='failed';status.elapsed_seconds=toc(timer);status.error_report=getReport(err,'extended','hyperlinks','off');writeStatus(out,status);rethrow(err);
end
end
function writeStatus(out,status)
f=fopen(fullfile(out,'ci_status.json'),'w');assert(f~=-1);c=onCleanup(@() fclose(f));fprintf(f,'%s\n',jsonencode(status));
end

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

