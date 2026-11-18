function fig3_outputs(out,c,models,step)
labels=string({models.name})';
n=numel(models);
loss=zeros(n,9);
parameters=table();
bounds=table();
residuals=[];
curves=[];
kinetics=[];
turns=[];
grid=[0,1e-7,1e-6,1e-5,8.55e-5,1e-4,1e-3,1e-2];
for m=1:n
    entry=models(m);
    p=entry.p;
    k=0;
    if strcmp(entry.pkMode,'dynamic'), k=p(end); end
    param=[p,entry.theta];
    pname=[entry.names,arrayfun(@(x) ['theta_',x],'a':'i','UniformOutput',false)];
    lo=[entry.lo,repmat(100,1,9)]; hi=[entry.hi,repmat(1e8,1,9)];
    models(m).count=numel(param);
    params=table(repmat(labels(m),numel(param),1),string(pname(:)),param(:), ...
        'VariableNames',{'model','parameter','estimate'});
    if isempty(parameters), parameters=params; else, parameters=[parameters;params]; end %#ok<AGROW>
    atBound=min(log(param./lo),log(hi./param))<=0.01;
    rows=table(params.model,params.parameter,lo(:),hi(:),atBound(:), ...
        'VariableNames',{'model','parameter','lower','upper','near_bound'});
    if isempty(bounds), bounds=rows; else, bounds=[bounds;rows]; end %#ok<AGROW>
    for j=1:9
        g=c(j);
        tt=unique([0;g.t;g.wash;linspace(0,max(g.wash,g.t(end)),401)']);
        [cf,cb]=fig3_pk(tt,g.dose,g.wash,k);
        if strcmp(entry.pkMode,'fixed'), cb(:)=g.platform; cf(:)=NaN; end
        rate=fig3_kill(cb,entry.kind,p(2:end));
        select=tt>=g.t(1) & tt<=g.t(end);
        mesh=fig3_mesh(tt(select),g.wash,step);
        if strcmp(entry.pkMode,'fixed')
            input=repmat(g.platform,size(mesh.nodes));
        else
            [~,input]=fig3_pk(mesh.nodes,g.dose,g.wash,k);
        end
        [base,den,A]=fig3_population(mesh,fig3_kill(input,entry.kind,p(2:end)),p(1),g.y(1));
        N=nan(size(tt)); logN=N;
        logN(select)=log(g.y(1))+A-log1p(den/entry.theta(j));
        N(select)=exp(logN(select));
        [~,ix]=ismember(g.t,tt);
        e=(N(ix)-g.y)/max(g.y);
        loss(m,j)=sum(e.^2);
        residuals=[residuals; repmat([m,j],numel(ix),1),g.t,g.y,N(ix),e]; %#ok<AGROW>
        curves=[curves;repmat([m,j],numel(tt),1),tt,cf,cb,repmat(g.platform,numel(tt),1),rate,N,logN]; %#ok<AGROW>
        before=g.t<=48;
        [observedMin,obsI]=min(g.y);
        [predMin,predI]=min(N(select)); st=tt(select);
        turns=[turns; m,j, g.t(obsI), observedMin, g.y(end)-observedMin, ...
            st(predI),predMin,N(ix(end))-predMin,sum(e(before).^2), ...
            sum(e(~before).^2),nnz(before)]; %#ok<AGROW>
    end
    assert(abs(sum(loss(m,:))-entry.sse)<2e-6,'Prediction SSE disagrees with optimization for %s.',entry.name);
end
writetable(parameters,fullfile(out,'parameters.csv'));
writetable(bounds,fullfile(out,'bounds.csv'));
summary=table((1:n)',labels,[models.count]',sum(loss,2), ...
    'VariableNames',{'model_id','model','free_parameters','SSE'});
for j=1:9
    summary.(['SSE_',char(96+j)])=loss(:,j);
    g=c(j);
    for k=unique([grid,arrayfun(@(m) m.p(end),models(7:9))])
        [peakT,peak,t95]=pkMetrics(g,k);
        [~,last]=fig3_pk(g.t(end),g.dose,g.wash,k);
        [~,t600]=fig3_pk(600,g.dose,g.wash,k);
        [~,w600]=fig3_pk(g.wash+600,g.dose,g.wash,k);
        [~,atWash]=fig3_pk(g.wash,g.dose,g.wash,k);
        kinetics=[kinetics;j,k,g.platform/g.dose,g.platform,peakT,peak,t95, ...
            g.t(end),last,100*(last/peak-1),100*(t600/peak-1),100*(w600/peak-1),...
            100*(atWash/g.platform-1),100*(peak/g.platform-1)]; %#ok<AGROW>
    end
end
writetable(summary,fullfile(out,'summary.csv'));
output('residuals',residuals,{'model_id','panel_id','time_h','observed_N','predicted_N','scaled_residual'});
output('curves',curves,{'model_id','panel_id','time_h','CF_nM','CB_nM','CB_const_nM','kill_h_inv','N','log_N'});
output('regrowth',turns,{'model_id','panel_id','observed_min_t_h','observed_min_N', ...
    'observed_end_minus_min','predicted_min_t_h','predicted_min_N','predicted_end_minus_min', ...
    'SSE_through_48h','SSE_after_48h','count_through_48h'});
output('platform',kinetics,{'panel_id','kBF','alpha','irreversible_platform_nM', ...
    'peak_time_h','peak_nM','first_at_95pct_irreversible_h','last_obs_h','CB_at_last_obs_nM', ...
    'last_obs_change_from_peak_pct','change_from_peak_at_600h_pct', ...
    'change_from_peak_600h_after_wash_pct','wash_vs_irreversible_pct','peak_vs_irreversible_pct'});
for family=1:3
    selected=3*(family-1)+(1:3);
    plotGrid('fit',curves,8,selected);
    plotGrid('residual',residuals,6,selected);
    plotGrid('CB',curves,5,selected);
    plotGrid('killing',curves,7,selected);
end
f=figure('Visible','off','Color','w','Position',[40,40,1150,520]);
for j=1:3
    m=j+6; prof=models(m).profile; positive=prof(:,1)>0;
    subplot(1,3,j);
    semilogx(prof(positive,1),prof(positive,2),'o-'); hold on
    yline(prof(1,2),'k--','zero');
    xline(models(m).p(end),':');
    title(labels(m),'Interpreter','none'); xlabel('kBF (h^{-1})'); ylabel('Conditional SSE'); grid on;
end
saveas(f,fullfile(out,'kbf_profiles.png')); close(f);

    function output(name,values,headers)
        writetable(array2table(values,'VariableNames',headers),fullfile(out,[name,'.csv']));
    end
    function plotGrid(name,values,column,selected)
        f=figure('Visible','off','Color','w','Position',[30,30,1400,920]);
        for panel=1:9
            g=c(panel);
            subplot(3,3,panel); hold on;
            colors=lines(3);
            for i=1:3
                rows=values(values(:,1)==selected(i) & values(:,2)==panel,:);
                plot(rows(:,3),rows(:,column),'Color',colors(i,:),'DisplayName',labels(selected(i)));
            end
            if strcmp(name,'fit'), plot(g.t,g.y,'k.','DisplayName','Fig. 3 counts'); end
            if strcmp(name,'CB'), yline(g.platform,'k--','DisplayName','Irreversible plateau'); end
            if strcmp(name,'residual'), yline(0,'k:','HandleVisibility','off'); end
            xline(g.wash,':','HandleVisibility','off');
            title(sprintf('%s | %g nM / %g h',g.panel,g.dose,g.wash));
            xlabel('Hours from dosing'); ylabel(name); grid on;
            if panel==1, legend('Interpreter','none','Location','best'); end
        end
        saveas(f,fullfile(out,sprintf('%s_%d.png', name,selected(1)))); close(f);
    end
end

function [tPeak,peak,t95]=pkMetrics(g,k)
% Peak for positive k is after washout. The zero-k limit is reached at infinity.
if k==0
    tPeak=Inf; peak=g.platform;
else
    kFE=1.25; kFB=.21;
    gap=sqrt((kFE+kFB-k)^2+4*kFB*k);
    fast=(kFE+kFB+k+gap)/2;
    slow=kFE*k/fast;
    [cf,cb]=fig3_pk(g.wash,g.dose,g.wash,k);
    slope=kFB*cf-k*cb;
    tPeak=g.wash+log1p(gap*slope/(slow*(fast*cb+slope)))/gap;
    [~,peak]=fig3_pk(tPeak,g.dose,g.wash,k);
end
t95=NaN;
if peak>=.95*g.platform
    hi=min(tPeak,g.wash+100);
    t95=fzero(@crossing,[0,hi]);
end
    function value=crossing(t)
        [~,cb]=fig3_pk(t,g.dose,g.wash,k);
        value=cb-.95*g.platform;
    end
end
