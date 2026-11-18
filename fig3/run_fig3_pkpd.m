function run_fig3_pkpd(mode, out)
% Independent Fig. 3 fit, using base MATLAB and no legacy output directory.
if nargin<1, mode='fig3-quick'; end
mode=validatestring(mode,{'fig3-quick','fig3-full'});
root=fileparts(mfilename('fullpath'));
if nargin<2, out=fullfile(fileparts(root),'fig3_results'); end
assert(~isfolder(out) && ~isfile(out),'Choose a new Fig. 3 output directory.');
mkdir(out);
diary(fullfile(out,'run.log'));
logCleanup=onCleanup(@() diary('off')); %#ok<NASGU>
rng(42,'twister');
clock=tic;
status=struct('status','running','mode',mode,'matlab_version',version, ...
    'commit',getenv('GITHUB_SHA'),'run_id',getenv('GITHUB_RUN_ID'),'mesh_step_h',0.5);
writeStatus(out,status);
try
    data=readtable(fullfile(root,'observations.csv'),'TextType','string');
    assert(isequal(unique(data.panel),string(('a':'i')')));
    step=status.mesh_step_h;
    for j=1:9
        D=data(data.panel==char(96+j),:);
        assert(height(D)>=2 && all(diff(D.time_h)>0) && all(D.cells>0));
        assert(D.time_h(1)>=0 && numel(unique(D.dose_nM))==1 && numel(unique(D.wash_h))==1);
        c(j)=struct('panel',D.panel(1),'dose',D.dose_nM(1),'wash',D.wash_h(1), ...
            't',D.time_h,'y',D.cells,'mesh',fig3_mesh(D.time_h,D.wash_h(1),step), ...
            'platform',fig3_platform(D.dose_nM(1),D.wash_h(1))); %#ok<AGROW>
        window(j,:)=[j,height(D),D.time_h(1),D.cells(1),D.time_h(end),max(D.cells)]; %#ok<AGROW>
    end
    writematrixWithHeader(out,'observation_windows',window, ...
        {'panel_id','count','first_h','N0','last_h','residual_scale_cells'});
    checks=fig3_checks();
    writetable(checks,fullfile(out,'numerical_checks.csv'));
    quick=strcmp(mode,'fig3-quick');
    nStarts=3; nGrid=7;
    if ~quick, nStarts=10; nGrid=25; end
    kinds={'linear','hill1','hill'};
    pkModes={'fixed','zero','dynamic'};
    models=struct([]);
    for family=1:3
        for type=1:3
            m=3*(family-1)+type;
            [names,lo,hi,guess]=parameterBounds(kinds{type},pkModes{family});
            seeds=guess;
            if family>1
                inherited=models(m-3).p;
                if family==3, inherited=[inherited,0.0049373889]; end %#ok<AGROW>
                seeds=[seeds;inherited]; %#ok<AGROW>
            end
            if type==3
                previous=models(m-1).p;
                seeds=[seeds; [previous(1:3),1,previous(4:end)]]; %#ok<AGROW>
            end
            label=[pkModes{family},'_',kinds{type}];
            [best,sse,runs,theta]=fit(c,kinds{type},pkModes{family},lo,hi,seeds, ...
                nStarts,NaN,label,out);
            models(m).name=label;
            models(m).kind=kinds{type};
            models(m).pkMode=pkModes{family};
            models(m).p=best;
            models(m).sse=sse;
            models(m).theta=theta;
            models(m).lo=lo;
            models(m).hi=hi;
            models(m).names=names;
            if family==3
                grid=unique([0,logspace(log10(lo(end)),log10(hi(end)),nGrid),best(end), ...
                    max(lo(end),best(end)/3),min(hi(end),3*best(end))]);
                % Skip redundant re-fit at the exact joint optimum only as a seed.
                profile=zeros(numel(grid),3+numel(best)+9);
                warm=best(1:end-1);
                for g=1:numel(grid)
                    if grid(g)==0
                        z=models(m-3);
                        p=z.p; v=z.sse; th=z.theta; flag=NaN;
                    else
                        profLabel=sprintf('%s_profile_%02d',label,g);
                        [p,v,pr,th]=fit(c,kinds{type},'dynamic',lo(1:end-1), ...
                            hi(1:end-1),[best(1:end-1);warm],2,grid(g),profLabel,out);
                        [~,ix]=min(pr(:,end-2)); flag=pr(ix,end-1);
                        warm=p;
                    end
                    profile(g,:)=[grid(g),v,flag,p,grid(g),th];
                end
                models(m).profile=profile;
                % A profile can discover a better basin for the unrestricted fit.
                [~,ix]=min(profile(2:end,2)); ix=ix+1;
                if profile(ix,2)<sse-1e-8
                    [p,v,extra,th]=fit(c,kinds{type},'dynamic',lo,hi, ...
                        [profile(ix,4:3+numel(best));best],2,NaN,[label,'_refine'],out);
                    runs=[runs;extra]; %#ok<AGROW>
                    if v<sse, models(m).p=p; models(m).sse=v; models(m).theta=th; end
                    profile(end+1,:)=[models(m).p(end),models(m).sse,NaN,models(m).p,models(m).theta]; %#ok<AGROW>
                    models(m).profile=sortrows(profile,1);
                end
                assert(min(profile(2:end,2))>=models(m).sse-1e-4, ...
                    'Conditional fit outperformed unrestricted fit for %s.',label);
                profileNames=[{'kBF_fixed','SSE','exitflag'},names,thetaNames()];
                writematrixWithHeader(out,[label,'_kbf_profile'],models(m).profile,profileNames);
            end
            models(m).runs=runs;
        end
    end
    fig3_outputs(out,c,models,step);
    save(fullfile(out,'fig3_fit.mat'),'models','data','checks','window','mode');
    status.status='success';
    status.model_count=numel(models);
    status.point_count=height(data);
    status.elapsed_seconds=toc(clock);
    writeStatus(out,status);
catch failure
    status.status='failed';
    status.error=getReport(failure,'extended','hyperlinks','off');
    status.elapsed_seconds=toc(clock);
    writeStatus(out,status);
    rethrow(failure);
end
end

function [names,lo,hi,p]=parameterBounds(kind,pkMode)
names={'kp'};
lo=1e-4; hi=0.12; p=0.031094959444983555;
switch kind
    case 'linear'
        names=[names,{'beta'}]; lo=[lo,1e-9]; hi=[hi,0.1]; p=[p,5e-4];
    otherwise
        names=[names,{'Kmax','Khalf'}]; lo=[lo,1e-6,1e-3]; hi=[hi,1,1e7]; p=[p,0.07,100];
        if strcmp(kind,'hill')
            names=[names,{'n'}]; lo=[lo,0.2]; hi=[hi,8]; p=[p,1.5];
        end
end
if strcmp(pkMode,'dynamic')
    names=[names,{'kBF'}]; lo=[lo,1e-7]; hi=[hi,10]; p=[p,0.0049373889];
end
end

function [best,sse,runs,theta]=fit(c,kind,pkMode,lo,hi,seeds,n,kFixed,label,out)
% Project the nine independent theta parameters, then fit the shared rates.
% Reuse the repository's log-parameter, bounded fminsearch convention.
l=log(lo); u=log(hi);
options=optimset('Display','off','MaxIter',650,'MaxFunEvals',1800, ...
    'TolX',1e-5,'TolFun',1e-8);
n=max(n,size(seeds,1)); d=numel(lo);
runs=zeros(n,2*d+12); sse=Inf;
if strcmp(pkMode,'fixed') || strcmp(pkMode,'zero') || isfinite(kFixed)
    for j=1:9
        if strcmp(pkMode,'fixed')
            c(j).CB=repmat(c(j).platform,size(c(j).mesh.nodes));
        else
            k=0; if isfinite(kFixed), k=kFixed; end
            [~,c(j).CB]=fig3_pk(c(j).mesh.nodes,c(j).dose,c(j).wash,k);
        end
    end
end
for i=1:n
    if i<=size(seeds,1)
        z0=log(seeds(i,:));
    elseif mod(i,3)==0
        z0=l+rand(1,d).*(u-l);
    else
        z0=log(seeds(1,:))+1.8*randn(1,d);
        if strcmp(pkMode,'dynamic') && ~isfinite(kFixed)
            z0(end)=l(end)+rand*(u(end)-l(end));
        end
    end
    z0=min(u,max(l,z0));
    [z,~,flag,info]=fminsearch(@objective,z0,options);
    [value,th]=objective(z);
    runs(i,:)=[exp(z0),exp(z),th,value,flag,info.funcCount];
    if value<sse, best=exp(z); sse=value; theta=th; end
    fprintf('%s start %d/%d SSE %.9g exit %d evaluations %d\n',label,i,n,value,flag,info.funcCount);
    progress=runs(1:i,:); %#ok<NASGU>
    save(fullfile(out,[label,'_starts.mat']),'progress','lo','hi','sse','best','theta');
end
assert(isfinite(sse) && sse<1e8,'No valid fit for %s.',label);
headers=[compose("start_%d",1:d),compose("end_%d",1:d), ...
    string(thetaNames()),"SSE","exitflag","function_evaluations"];
writematrixWithHeader(out,[label,'_starts'],runs,cellstr(headers));

    function [v,th]=objective(z)
        th=nan(1,9);
        outside=max(l-z,0)+max(z-u,0);
        if any(outside>0), v=1e8+sum(outside.^2); return; end
        p=exp(z);
        v=0;
        for q=1:9
            if isfield(c,'CB')
                CB=c(q).CB;
            else
                [~,CB]=fig3_pk(c(q).mesh.nodes,c(q).dose,c(q).wash,p(end));
            end
            rate=fig3_kill(CB,kind,p(2:end));
            [base,density]=fig3_population(c(q).mesh,rate,p(1),c(q).y(1));
            loss=@(x) sum(((base./(1+density*exp(-x))-c(q).y)/max(c(q).y)).^2);
            [th(q),part]=fitTheta(loss);
            v=v+part;
        end
    end
end

function [theta,value]=fitTheta(loss)
persistent grid options
if isempty(grid)
    grid=linspace(log(100),log(1e8),9);
    options=optimset('Display','off','TolX',1e-5);
end
values=arrayfun(loss,grid);
[value,k]=min(values);
winner=grid(k);
basins=find(values(2:end-1)<=values(1:end-2) & values(2:end-1)<values(3:end))+1;
for j=unique([k,basins])
    [x,v]=fminbnd(loss,grid(max(1,j-1)),grid(min(numel(grid),j+1)),options);
    if v<value, winner=x; value=v; end
end
theta=exp(winner);
end

function names=thetaNames()
names=arrayfun(@(x) ['theta_',x],'a':'i','UniformOutput',false);
end

function writeStatus(out,status)
fid=fopen(fullfile(out,'ci_status.json'),'w');
assert(fid>=0);
guard=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(status,'PrettyPrint',true));
end

function writematrixWithHeader(out,name,values,names)
writetable(array2table(values,'VariableNames',names),fullfile(out,[name,'.csv']));
end
