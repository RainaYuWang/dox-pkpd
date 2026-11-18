function checks = fig3_checks()
% Compare the fitted solver with independent, split ode15s calculations.
names={}; errors=[]; tolerances=[];
kp=.031; theta=35000; N0=3000;
kFE=1.25; kFB=.21; input=3.5e-6*250/.004629;
options=odeset('RelTol',1e-9,'AbsTol',1e-11);
for wash=[6,12,24]
    for k=[0,1e-4,.005,10]
        dose=1250;
        t=unique([0,1,4,wash,wash+1,48,100,300,720])';
        initialTime=4;
        pk1=ode15s(@(~,y) pkRhs(y,dose,k),[0,wash],[0;0],options);
        pk2=ode15s(@(~,y) pkRhs(y,0,k),[wash,720],deval(pk1,wash),options);
        [cf,cb]=fig3_pk(t,dose,wash,k);
        first=t<=wash;
        actual=[deval(pk1,t(first)),deval(pk2,t(~first))]';
        add('analytic_PK_vs_ode15s',scaled([cf,cb],actual),2e-7);
        [leftCF,leftCB]=fig3_pk(wash-1e-8,dose,wash,k);
        [rightCF,rightCB]=fig3_pk(wash+1e-8,dose,wash,k);
        add('washout_CF_CB_continuity',max(abs([rightCF-leftCF,rightCB-leftCB]))/max([cf;cb]),1e-7);
        if k==0
            add('irreversible_platform',abs(cb(end)/fig3_platform(dose,wash)-1),1e-10);
        end
        delay=t(t>=initialTime);
        for kind={'linear','hill1','hill'}
            if strcmp(kind{1},'linear'), p=0.0005; else, p=[.08,100,1.5]; end
            m=fig3_mesh(delay,wash,.5);
            [~,bound]=fig3_pk(m.nodes,dose,wash,k);
            [base,den]=fig3_population(m,fig3_kill(bound,kind{1},p),kp,N0);
            pred=base./(1+den/theta);
            % Initialize the combined state with PK already evolved from t=0.
            init=[deval(pk1,initialTime); log(N0)];
            pre=ode15s(@(~,y) coupled(y,dose,k,kind{1},p),[initialTime,wash],init,options);
            atWash=deval(pre,wash);
            post=ode15s(@(~,y) coupled(y,0,k,kind{1},p),[wash,720],atWash,options);
            before=delay<=wash;
            reference=exp([deval(pre,delay(before),3),deval(post,delay(~before),3)]');
            add('population_vs_ode15s',scaled(pred,reference),2e-7);
            add('N_washout_continuity',abs(exp(deval(pre,wash-1e-8,3)-deval(post,wash+1e-8,3))-1),1e-7);
            add('nonnegative_states_and_rate',max(0,-min([cf;cb;pred;fig3_kill(cb,kind{1},p)])),1e-10);
        end
    end
end
t=[0;12;24;120;600];
m=fig3_mesh(t,12,.5);
for kill=[0,kp,.01,.08]
    [base,den]=fig3_population(m,kill+zeros(size(m.nodes)),kp,N0);
    pred=base./(1+den/theta);
    a=kp-kill;
    if a==0, integral=t; else, integral=expm1(a*t)/a; end
    exact=N0*exp(a*t)./(1+(kp*N0/theta)*integral);
    if kill==0
        add('zero_killing_logistic',scaled(pred,theta./(1+(theta/N0-1)*exp(-kp*t))),1e-10);
    end
    add('constant_input_analytic',scaled(pred,exact),1e-10);
end
cb=[0;logspace(-5,5,100)'];
beta=.0005; Ks=200;
add('Hill1_old_saturation',scaled(fig3_kill(cb,'hill1',[beta*Ks,Ks]), ...
    beta*cb./(1+cb/Ks)),1e-12);
assert(all(errors<=tolerances));
checks=table(string(names(:)),errors(:),tolerances(:),'VariableNames',{'check','error','tolerance'});

    function dy = pkRhs(y,ce,k)
        dy=[input*ce-(kFE+kFB)*y(1)+k*y(2);kFB*y(1)-k*y(2)];
    end
    function dy = coupled(y,ce,k,kind,p)
        dy=[pkRhs(y,ce,k);kp*(1-exp(y(3))/theta)-fig3_kill(y(2),kind,p)];
    end
    function add(name,errorValue,tolerance)
        names{end+1}=name; errors(end+1)=errorValue; tolerances(end+1)=tolerance;
        assert(isfinite(errorValue) && errorValue<=tolerance, ...
            '%s error %.9g exceeds tolerance %.9g.',name,errorValue,tolerance);
    end
end

function e=scaled(x,y)
e=max(abs(x(:)-y(:)))/max(abs(y(:)));
end
