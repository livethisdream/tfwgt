function [etsol,ezsol,mutsol,muzsol] = runSolver(Smeas,wval,etguess,ezguess,mutguess,muzguess)

global eps0 mu0 solveCase alg ket kmut;

% set the options for the solver
% options=optimset('Display','off','PlotFcns',...
%         {@optimplotfunccount @optimplotx @optimplotfval @optimplotstepsize});
if strcmp(alg,'LM')==1
    options=optimset('Display','iter','Algorithm',{'levenberg-marquardt',1e-6});
elseif strcmp(alg,'TRR')==1
    options=optimset('Display','iter');
    lb=-50;
    ub=50;
end



% eval the solveCase flag
switch solveCase
    case 1 % isotropic, dielectric, non-mag -> ez=et and muz=mut=mu0
        if strcmp(alg,'LM')==1
            Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess)],...
                wval,Smeas,[],[],options);
        elseif strcmp(alg,'TRR')==1 
            Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess)],...
                wval,Smeas,[0 lb],[ub 0],options);
        end
    etsol=Y(1)+1j*Y(2);
    ezsol=Y(1)+1j*Y(2);
    mutsol=1;
    muzsol=1;


    case 2 % isotropic, non-dielectric, mag -> ez=et=eps0 and muz=mut
        if strcmp(alg,'LM')==1
            Y=lsqcurvefit(@Sparams,[real(mutguess) imag(mutguess)],...
                wval,Smeas,[],[],options);
        elseif strcmp(alg,'TRR')==1 
            Y=lsqcurvefit(@Sparams,[real(mutguess) imag(mutguess)],...
                wval,Smeas,[0 lb],[ub 0],options);
        end
    etsol=1;
    ezsol=1;
    mutsol=Y(1)+1j*Y(2);
    muzsol=Y(1)+1j*Y(2);

    case 3 % isotropic, dielectric, mag -> ez=et and muz=mut
        if strcmp(alg,'LM')==1
            Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess) ...
                real(mutguess) imag(mutguess)],...
                wval,Smeas,[],[],options);
        elseif strcmp(alg,'TRR')==1
        Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess) ...
            real(mutguess) imag(mutguess)],...
            wval,Smeas,[0 lb 0 lb],[ub 0 ub 0],options);
        end
    etsol=Y(1)+1j*Y(2);
    ezsol=Y(1)+1j*Y(2);
    mutsol=Y(3)+1j*Y(4);
    muzsol=Y(3)+1j*Y(4);

    case 4 % uniaxial, dielectric, non-mag -> ez,et and muz=mut=mu0
        if strcmp(alg,'LM')==1
            Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess)...
                real(ezguess) imag(ezguess)],...
                wval,Smeas,[],[],options);
        elseif strcmp(alg,'TRR')==1
%             Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess)...
%                 real(ezguess) imag(ezguess)],...
%                 wval,Smeas,[lb lb lb lb],[ub ub ub ub],options);
              Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess)...
                real(ezguess) imag(ezguess)],...
                wval,Smeas,[0 lb 0 lb],[ub 0 ub 0],options);
        end
    etsol=Y(1)+1j*Y(2);
    ezsol=Y(3)+1j*Y(4);
    mutsol=1;
    muzsol=1;

    case 5 % uniaxial, non-dielectric, mag -> ez=et=eps0 and muz,mut
        if strcmp(alg,'LM')==1
            Y=lsqcurvefit(@Sparams,[real(mutguess) imag(mutguess)...
                real(muzguess) imag(muzguess)],...
                wval,Smeas,[],[],options);
        elseif strcmp(alg,'TRR')==1
            Y=lsqcurvefit(@Sparams,[real(mutguess) imag(mutguess)...
                real(muzguess) imag(muzguess)],...
                wval,Smeas,[0 lb 0 lb],[ub 0 ub 0],options);
        end
    etsol=1;
    ezsol=1;
    mutsol=Y(1)+1j*Y(2);
    muzsol=Y(3)+1j*Y(4);

    case 6 % uniaxial, dielectric, mag -> ez,et and muz,mut
        if strcmp(alg,'LM')==1
            Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess)...
                real(ezguess) imag(ezguess) real(mutguess) imag(mutguess)...
                real(muzguess) imag(muzguess)],...
                wval,Smeas,[],[],options);
        elseif strcmp(alg,'TRR')==1
            Y=lsqcurvefit(@Sparams,[real(etguess) imag(etguess)...
                real(ezguess) imag(ezguess) real(mutguess) imag(mutguess)...
                real(muzguess) imag(muzguess)],...
                wval,Smeas,[0 lb 0 lb 0 lb 0 lb],[ub 0 ub 0 ub 0 ub 0],options);
        end
    etsol=Y(1)+1j*Y(2);
    ezsol=Y(3)+1j*Y(4);
    mutsol=Y(5)+1j*Y(6);
    muzsol=Y(7)+1j*Y(8);
    
    case 7 % uniaxial, dielectric, mag with known et and mut
        if strcmp(alg,'LM')==1
            Y=lsqcurvefit(@Sparams,[real(ezguess) imag(ezguess) real(muzguess) imag(muzguess)],...
                wval,Smeas,[],[],options);
        elseif strcmp(alg,'TRR')==1
            Y=lsqcurvefit(@Sparams,[real(ezguess) imag(ezguess) real(muzguess) imag(muzguess)],...
                wval,Smeas,[0 lb 0 lb],[ub 0 ub 0],options);
        end
    etsol=ket;
    ezsol=Y(1)+1j*Y(2);
    mutsol=kmut;
    muzsol=Y(3)+1j*Y(4);

end 