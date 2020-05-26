function foptpath(numIterations,num,rep)

% This function finds the fastest base running path from home to second
% base on a baseball field based on my ability to run. It assumes that when
% I run, the magnitude of jerk is always zero.
%    Please see my research paper (link below) for a detailed explanation
%    of this function. For a brief explanation of how it goes about finding
%    the fastest path, see my video (link below).
%
%    Research paper: https://github.com/alensugimoto/optimum-path-finding-function/blob/master/ResearchPaper.pdf
%    Video: https://github.com/alensugimoto/optimum-path-finding-function/blob/master/ExplanationOfFunction.mp4

% WORKS CITED
% 
% Beltman, Floris. “Optimization of Ideal Racing Line.” PDF. Amsterdam: Floris Beltman,
%     2008. https://beta.vu.nl/nl/Images/werkstuk-beltman_tcm235-91313.pdf.
% Cardamone, L., Loiacono, D., Lanzi, P.L., and Bardelli, A.P. “Searching for the Optimal
%     Racing Line Using Genetic Algorithms.” ResearchGate. September 2010.
%     https://www.researchgate.net/profile/Daniele_Loiacono/publication/224180066_S
%     earching_for_the_Optimal_Racing_Line_Using_Genetic_Algorithms/links/00b7d
%     525e4a7499325000000/Searching-for-the-Optimal-Racing-Line-Using-GeneticAlgorithms.pdf.
% Carozza, D., Johnson, S. and Morgan, F. “Baserunner’s Optimal Path.” Math
%     Intelligencer 32, no. 4 (2010): 10-15. https://doi.org/10.1007/s00283-009-9106-2.
% O’Haver, Tom. “Fast smoothing function.” MATLAB Central File Exchange.
%     MathWorks. Updated February 14, 2017.
%     https://www.mathworks.com/matlabcentral/fileexchange/19998-fast-smoothingfunction.
% Tipping, Michael E. (Girton, GB), Hatton, Mark Andrew (Eye, GB), and Herbrich, Ralf
%     (Cambridge, GB). 2013. Racing line optimization. US Patent 8,393,944 B2, filed
%     January 13, 2011, and issued June 9, 2011.
% Xiong, Ying. Racing Line Optimization. Cambridge: Massachusetts Institute of
%     Technology, Computation for Design and Optimization Program, 2010.

% UNITS
% 
% feet and seconds

% INPUT ARGUMENTS
% 
% numIterations >>
% Compare [value of numIterations] paths before termination.
% 
% num >>
% ~1 == normal algorithm (recommended)
%  1 == precise but slow algorithm
% 
% rep >> 
% After every [value of rep] paths, display fastest path found so far.


switch nargin
    case 3
    case 2
        rep = 1;
    case 1
        num = 0;
        rep = 1;
    otherwise
        numIterations = 1;
        num = 0;
        rep = 1;
end

% The following values are obtained from a collection of data and are based
% on my running ability.
maxPerpDist = 30.5;
v_max = 268/12;
a_max = v_max/3.11;
a_min = -v_max/(9.63-8.22);
sqac = [3.73 1.75];

% the body of the function starts here
numPerpPoints = round(maxPerpDist+1);
numParaPoints = 6;
distBetwBases = 89;
pathPrcsn = 721;
desarcstep = 0.06;
paraStep = distBetwBases/(numParaPoints - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize graph

f = figure(1);
set(f,'WindowStyle','Normal')
f.WindowState = 'maximized';
set(f,'WindowStyle','Docked')

% bases
b = plot(0,0,'kd',distBetwBases,0,'k^',distBetwBases,distBetwBases,'ks');
b(1).MarkerSize = 10;
b(2).MarkerSize = 10;
b(3).MarkerSize = 10;
hold off

axis([-10 125 -35 100])
daspect([1 1 1])
xlabel('x (ft)')
xticks([0 89])
xticklabels({'0','89'})
ylabel('y (ft)')
yticks([0 89])
yticklabels({'0','89'})
title({'A graph of the fastest path from home to second base','produced after comparing 0 path(s)'})
set(gca,'FontSize',13)
legend([b(1) b(2) b(3)],{'home base','inner corner of first base','second base',...
    },'Location','NorthWest')
grid on
grid minor

pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%

pnt = linspace(0,-maxPerpDist,numPerpPoints);

for i = 1:numIterations
    nx1 = linspace(0,distBetwBases,numParaPoints);
    nx2 = nx1;

    % home to first
    y1 = zeros(1,numParaPoints);
    y1(2) = pnt(randi(end));

    for j = 2:numParaPoints-2
        k = 1;
        while round(pnt(k)-y1(j),4) > round(-y1(j)/(numParaPoints-j),4)
            k = k+1;
        end

        l = numPerpPoints;
        while round(pnt(l)-y1(j),4) < round(y1(j)-y1(j-1),4)
            l = l-1;
        end

        y1(j+1) = pnt(randi([k,l]));
    end

    mBeforeFirst = paraStep/y1(end-1);
    ny1 = y1;
    y1 = fastsmooth(y1,3,1,0);

    % first to second
    y2 = zeros(1,numParaPoints);

    for j = 1:numParaPoints-2
        k = 1;
        while round(pnt(k)-y2(j),4) > round(-y2(j)/(numParaPoints-j),4)
            k = k+1;
        end

        l = numPerpPoints;
        if j == 1 && ~isinf(mBeforeFirst)
            while round(pnt(l)/paraStep,4) < round(mBeforeFirst,4)
                l = l-1;
            end
        elseif j > 1
            while round(pnt(l)-y2(j),4) < round(y2(j)-y2(j-1),4)
                l = l-1;
            end
        end

        y2(j+1) = pnt(randi([k,l]));
    end

    ny2 = y2;
    y2 = fastsmooth(y2,3,1,0);

    x1 = nx1;
    x2 = nx2;
    
    % smoothing first
    theta = atan((-ny2(2)+(x1(end-1)-distBetwBases))/(x2(2)+ny1(end-1)));
    ycomp1 = (x1(end-1)-distBetwBases)*sin(theta)+ny1(end-1)*cos(theta);
    ycomp2 = -ny2(2)*sin(theta)+x2(2)*cos(theta);
    result = (ycomp1+ycomp2)/3;
    xcomp = -1*result*sin(-theta);
    ycomp = result*cos(-theta);

    deltax = -xcomp;
    deltay = -ycomp;
    
    for k = 2:5
        y1(k) = y1(k)+deltay*((k-1)/5);
        x1(k) = x1(k)+deltax*((k-1)/5);
    end

    deltax = -ycomp;
    deltay = xcomp;
    
    for k = 2:5
        y2(k) = y2(k)+deltay*((6-k)/5);
        x2(k) = x2(k)+deltax*((6-k)/5);
    end

    % make pchips
    pp1 = pchip(x1,y1);
    pp2 = pchip(x2,y2);

    % smooth path through first
    mA = -1/pp2.coefs(1,3);
    mB = -1/(3*pp1.coefs(end,1)*(x1(end)-x1(end-1))^2+2*pp1.coefs(end,2)*(x1(end)-...
        x1(end-1))+pp1.coefs(end,3));

    rise = 1/sqrt(1+mB^2) + 1/sqrt(1+1/mA^2);
    run = 1/sqrt(1+1/mB^2) + 1/sqrt(1+mA^2);
    mav = rise/run;

    f_p = pp1.coefs(end,3);
    A = [3*(x1(end)-x1(end-1))^2 2*(x1(end)-x1(end-1)); (x1(end)-x1(end-1))^3 (x1(end)...
        -x1(end-1))^2];
    B = [mav-f_p; -y1(end-1)-f_p*(x1(end)-x1(end-1))];
    x = linsolve(A,B);
    pp1.coefs(end,:) = [x(1,:) x(2,:) f_p y1(end-1)];

    f_p = pp2.coefs(2,3);
    c = -1/mav;
    A = [3*x2(2)^2 2*x2(2); x2(2)^3 x2(2)^2];
    B = [f_p-c; y2(2)-c*x2(2)];
    x = linsolve(A,B);
    pp2.coefs(1,:) = [x(1,:) x(2,:) c 0];

    coefs = [pp1.coefs; pp2.coefs];

    % find arcs for every node
    arclen = zeros(1,pp1.pieces+pp2.pieces);

    for k = 1:pp1.pieces+pp2.pieces
        if k <=pp1.pieces
            xq = linspace(0,x1(k+1)-x1(k));
        else
            xq = linspace(0,x2(k+1-pp1.pieces)-x2(k-pp1.pieces));
        end
        xq_i = xq(1:end-1);
        xq_f = xq(2:end);
        expr_i = sqrt(1+(3.*coefs(k,1).*xq_i.^2 + 2.*coefs(k,2).*xq_i + coefs(k,3)).^2);
        expr_f = sqrt(1+(3.*coefs(k,1).*xq_f.^2 + 2.*coefs(k,2).*xq_f + coefs(k,3)).^2);
        if k == 1
            arclen(k) = sum(0.5.*(xq_f-xq_i).*(expr_f+expr_i),2);
        else
            arclen(k) = arclen(k-1) + sum(0.5.*(xq_f-xq_i).*(expr_f+expr_i),2);
        end
    end

    % find speed vs distance graph
    a = [x1(1:numParaPoints-1) x2(1:numParaPoints-1)+distBetwBases];
    timePrcsn = round((arclen(end)/desarcstep)+1);
    arcstep = arclen(end)/(timePrcsn-1);
    s = linspace(0,arclen(end),timePrcsn);
    v = zeros(1,timePrcsn);
    v0 = v;
    prevPiece = 0;
    xprev = 0;
    yprev = sqrt(1+coefs(1,3)^2);

    for k = 1:timePrcsn
        pieces = find(sort([arclen s(k)])==s(k));
        piece = pieces(1);

        if prevPiece ~= piece
            c = coefs(piece,:);
        end

        if k ~= 1
            x = s(k)-s(k-1)+xprev;
            x0 = -1;
        else
            x = s(k);
        end

        if x ~= 0
            while round(x0,4) ~= round(x,4) 
                y = sqrt(1+(3*c(1)*(x-a(piece)).^2+2*c(2)*(x-a(piece))+c(3)).^2);
                if num == 1
                    expr = @(x) sqrt(1+(3*c(1)*(x-a(piece)).^2+2*c(2)*(x-a(piece))+c(3)).^2);
                    intfun = -arcstep + integral(expr,xprev,x);
                    intfun_p = y;
                    x0 = x;
                    x = x0 - (intfun/intfun_p);
                else
                    fun = 0.5*(x-xprev)*(y+yprev) - arcstep;
                    fun_p = 0.5*((y+yprev) + (x-xprev)*(3*c(1)*(x-a(piece))^2+2*c(2)*(x...
                        -a(piece))+c(3))*(6*c(1)*(x-a(piece))+2*c(2))/y);
                    x0 = x;
                    x = x0 - (fun/fun_p);
                end
            end
        end

        f_2p = 6*c(1)*(x-a(piece))+2*c(2);
        f_p = 3*c(1)*(x-a(piece))^2+2*c(2)*(x-a(piece))+c(3);

        r = (sqrt(1+f_p^2)^3)/(abs(f_2p));
        v1 = sqac(1)*sqrt(r)+sqac(2);

        v0(k) = min(v1,v_max);

        % correcting the graph
        if k > 1
            if k~=timePrcsn
                v(k) = min(v1,v_max);
            end
            for l = k-1:-1:1
                v_maxSlope = sqrt(v(l)^2+2*a_max*(s(l+1)-s(l)));
                v_minSlope = sqrt(v(l+1)^2-2*a_min*(s(l+1)-s(l)));
                if v(l+1) > v_maxSlope && l~=timePrcsn-1
                    v(l+1) = v_maxSlope;
                    break
                elseif v(l) > v_minSlope && l~=1
                    v(l) = v_minSlope;
                else
                    break
                end
            end
        end

        prevPiece = piece;
        xprev = x;
        yprev = sqrt(1+(3*c(1)*(x-a(piece))^2+2*c(2)*(x-a(piece))+c(3))^2);
    end

    t = 0;
    for k = 2:timePrcsn-2
        m = (v(k+1)-v(k))/arcstep;
        if m~=0
            t = t + log(v(k+1)/v(k))/m;
        else
            t = t + arcstep/v(k);
        end
    end

    if i==1 || t<=t_min
        X1 = x1;
        Y1 = y1;
        X2 = x2;
        Y2 = y2;
        PP1 = pp1;
        PP2 = pp2;
        if i==1 || t<t_min
            bests = 0;
        end
        t_min = t;
    end
        
    bests = bests+1;

    if (rem(i,rep)==0) || i==numIterations
        figure(1)
            
        % random paths
        xq = linspace(0,x1(end-1),pathPrcsn);
        xq2 = linspace(x2(2),distBetwBases,pathPrcsn);
        p1 = pchip(x1,y1,xq);
        p2 = pchip(x2,y2,xq2);

        rand1 = plot(xq,p1,'k-',-p2+distBetwBases,xq2,'k-');
        hold on

        xq = linspace(0,x1(end)-x1(end-1));
        xq2 = linspace(0,x2(2)-x2(1));
        rand2 = plot(xq+x1(end-1),polyval(pp1.coefs(end,:),xq),'k-',-polyval(...
            pp2.coefs(1,:),xq2)+distBetwBases,xq2,'k-');
        hold on
        
        rand1(1).Color = [0.7 0.7 0.7];
        rand1(2).Color = [0.7 0.7 0.7];
        rand2(1).Color = [0.7 0.7 0.7];
        rand2(2).Color = [0.7 0.7 0.7];
        rand1(1).LineWidth = 2;
        rand1(2).LineWidth = 2;
        rand2(1).LineWidth = 2;
        rand2(2).LineWidth = 2;
        
        % optimums
        xq = linspace(0,X1(end-1),pathPrcsn);
        xq2 = linspace(X2(2),distBetwBases,pathPrcsn);
        p1 = pchip(X1,Y1,xq);
        p2 = pchip(X2,Y2,xq2);
        
        opt1 = plot(xq,p1,'k-',-p2+distBetwBases,xq2,'k-');
        hold on

        xq = linspace(0,X1(end)-X1(end-1));
        xq2 = linspace(0,X2(2)-X2(1));
        opt2 = plot(xq+X1(end-1),polyval(PP1.coefs(end,:),xq),'k-',-polyval(PP2.coefs(1,:),...
            xq2)+distBetwBases,xq2,'k-');
        hold on
        
        opt1(1).LineWidth = 2;
        opt1(2).LineWidth = 2;
        opt2(1).LineWidth = 2;
        opt2(2).LineWidth = 2;

        % bases
        b = plot(0,0,'kd',distBetwBases,0,'k^',distBetwBases,distBetwBases,'ks');
        b(1).MarkerSize = 10;
        b(2).MarkerSize = 10;
        b(3).MarkerSize = 10;
        hold off
        
        axis([-10 125 -35 100])
        daspect([1 1 1])
        xlabel('x (ft)')
        xticks([0 89])
        xticklabels({'0','89'})
        ylabel('y (ft)')
        yticks([0 89])
        yticklabels({'0','89'})
        title({'A graph of the fastest path from home to second base',[ 'produced after comparing '...
            num2str(i) ' path(s)']})
        set(gca,'FontSize',13)
        legend([b(1) b(2) b(3) rand1(1) opt1(1)],{'home base','inner corner of first base',...
            'second base','randomly selected paths','fastest path found so far'},'Location',...
            'NorthWest')
        grid on
        grid minor
        
        text(30,36,{'Fastest path''s time: ',['\fontsize{20}' num2str(t_min) ' sec']})
    end
end

fprintf('Fastest path''s time: %4.2f sec (2DP)\n', t_min);

end
