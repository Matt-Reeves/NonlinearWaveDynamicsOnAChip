function Q = rk4(simgrids,simparams)
%xt = @(t) -pi+0.001*(1-cos(t));


        Q1=simparams.dt.*evolve(...
            simparams.t, ...
            simgrids.Q, ...
            simgrids.k,...
            simparams.N,...
            simparams.d,...
            simparams.kappa,...
            simparams.Gamma,...
            simparams.visc);

       % Q1(1) = xt(simparams.t);

        Q2=simparams.dt.*evolve(simparams.t+simparams.dt/2,simgrids.Q+Q1/2,simgrids.k,simparams.N,simparams.d,simparams.kappa,simparams.Gamma,simparams.visc);
       % Q2(1) = xt(simparams.t+simparams.dt/2);
        Q3=simparams.dt.*evolve(simparams.t+simparams.dt/2,simgrids.Q+Q2/2,simgrids.k,simparams.N,simparams.d,simparams.kappa,simparams.Gamma,simparams.visc);
        %Q3(1) = xt(simparams.t+simparams.dt/2);
        Q4=simparams.dt.*evolve(simparams.t+simparams.dt,simgrids.Q+Q3,simgrids.k,simparams.N,simparams.d,simparams.kappa,simparams.Gamma,simparams.visc);
       % Q4(1) = xt(simparams.t+simparams.dt);

        Q=simgrids.Q+(Q1+2*(Q2+Q3)+Q4)/6;
       % simparams.t = simparams.t+simparams.dt;

    return