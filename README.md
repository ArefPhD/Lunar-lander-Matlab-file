# Lunar-lander-Matlab-file

function lula
% LULA
% lunar lander
% Use arrow keys to operate the spaceship.
% tillmann.stuebler@akka.eu
%
    function callback(h,e)
        switch h
            case f
                pressed=strcmp(e.EventName,'KeyPress');
                switch e.Key
                    case 'leftarrow'
                        left=pressed;
                        if pressed
                            leftright=-1;
                        else
                            leftright=double(right);
                        end
                    case 'rightarrow'
                        right=pressed;
                        if pressed
                            leftright=1;
                        else
                            leftright=-left;
                        end
                    case 'uparrow'
                        up=pressed;
                        if pressed
                            updown=1;
                        else
                            updown=-down;
                        end
                    case 'downarrow'
                        down=pressed;
                        if pressed
                            updown=-1;
                        else
                            updown=double(up);
                        end
                end
        end
    end
    function dy=der(y)
        
        r1=y(1:2);
        v1=y(3:4);
        theta1=y(5);
        omega1=y(6);
        thrust1=y(7);
        m_fuel1=y(8);
        
        g_x=r1(1)+(-sin(theta1)*gx([3 10])+cos(theta1)*gy([3 10]));
        g_y=r1(2)+(cos(theta1)*gx([3 10])+sin(theta1)*gy([3 10]));
        vert_dist=g_y-fsurf(g_x);
        g_vx=v1(1)+omega1*(-cos(theta1)*gx([3 10])-sin(theta1)*gy([3 10]));
        g_vy=v1(2)+omega1*(-sin(theta1)*gx([3 10])+cos(theta1)*gy([3 10]));
        
        %[g_vx g_vy]
        
        surf_slope=fgrad(g_x);
        ey=surf_slope./(1+surf_slope.^2).^.5;
        ex=(1-ey.^2).^.5;
        v_lat=g_vx.*ex+g_vy.*ey;
        v_lon=g_vy.*ex-g_vx.*ey;
        K=-min(0,vert_dist);
        
        F_dissx=(vert_dist<0).*(-3e3*v_lat.*ex+1e3*v_lon.*ey);
        F_dissy=(vert_dist<0).*(-3e3*v_lat.*ey-1e3*v_lon.*ex);
        
        F_x=-1e4*K.*ey+F_dissx;
        F_y=1e4*K.*ex+F_dissy;
        
        domega=-.1*leftright+sum((g_x-r1(1)).*F_y-(g_y-r1(2)).*F_x)/((m+m_fuel1)*D^2);
        dthrust=max(2e3,thrust1)*updown;
        dtheta=omega1;
        dr=v1;
        dv=((m_fuel1>0)*thrust1*[cos(theta1);sin(theta1)]+[sum(F_x);sum(F_y)])./(m+m_fuel1)+[0;-g+v1(1)^2/R];
        dfuel=-thrust1/I_sp;
        ddiss=-.5*sum(g_vx.*F_dissx+g_vy.*F_dissy);
        
        dy=[dr;dv;dtheta;domega;dthrust;dfuel;ddiss];
    end
    function showhs
        % show high score list
        if exist(hslist,'file')
            hsfid=fopen(hslist,'r');
            hsc=textscan(hsfid,'%s %f %f %q');
            fclose(hsfid);
            [user,kine,remfuel,com]=deal(hsc{:});
            [uuser,~,iu]=unique(user);
            hsfig=figure('Name','LUNAR LANDER HIGH SCORES','MenuBar','none','NumberTitle','off','Color','k');
            hsa=subplot(1,1,1,'Parent',hsfig);
            cols=jet(numel(uuser));
            for k=1:numel(uuser)
                plot(hsa,kine(iu==k),remfuel(iu==k),'o','Color',cols(k,:),'MarkerFaceColor',cols(k,:),'DisplayName',uuser{k});
                hold(hsa,'on');
            end
            for k=find(~cellfun(@isempty,com))'
                text(hsa,kine(k),remfuel(k),sprintf(' \\leftarrow %s',com{k}),'HorizontalAlignment','left','VerticalAlignment','middle','Color','w');
            end
            hold(hsa,'off');
            xlabel(hsa,'residual kinetic energy, J');
            ylabel(hsa,'remaining fuel, kg');
            hsa.Color='k';
            hsa.XColor='w';
            hsa.YColor='w';
            hsa.GridColor='w';
            hsa.XScale='log';
            legend(hsa,'show','Location','eastoutside','Color','k','TextColor','w');
            grid(hsa,'on');
        end
    end
hslist='lulahs'; % high score list
up=false;
down=false;
left=false;
right=false;
leftright=0;
updown=0;
g=1.62; % moon gravity
m=1e3; % lander mass
D=2; % lander mass distribution radius
I_sp=5e2; %specific momentum
R=1738e3; % moon radius
nom_thrust=2e4; % nominal thrust
% initial conditions
r=[0;5e2]; % position
v=[-200;-30]; % velocity
theta=1; % angle
omega=0; % rotational speed
thrust=0;
m_fuel=1e3;
diss=0;
% lander polygon
fx=[61;84;84;70;77;72;52;30;15]*3e-2;
fx=[fx;-flip(fx)];
fy=[38;57;117;127;147;185;231;235;259]*3e-2-3;
fy=[fy;flip(fy)];
% landing gear
gx=[84;105;160;130;84;105]*3e-2;
gx=[gx;nan;-gx];
gy=[117;110;0;57;57;110]*3e-2-3;
gy=[gy;nan;gy];
% moon surface altitude profile
n=1e3;
xa=[0;cumsum(exp(randn(n-1,1)))];
ya=cumsum(sinh(randn(n,1)));
ya=movmedian(ya,5);
xg=linspace(0,1,1e4)';
ya=movmean(interp1(xa./xa(end),ya,xg),10);
xg=1e3*xg;
yg=ya-min(ya)+10;
f=figure('Name','LUNAR LANDER','MenuBar','none','Color','k','NumberTitle','off','KeyPressFcn',@callback,'KeyReleaseFcn',@callback);
a0=axes(f,'Position',[0 0 1 1]);
a1=axes(f,'Position',[0 0 1 1]);
a2=axes(f,'Position',[.75 .65 .2 .3]);
xg=[-2*xg(end)+xg;flip(-xg);xg;2*xg(end)-flip(xg)];
yg=[yg;flip(yg);yg;flip(yg)];
[xg,j]=unique(xg);
yg=yg(j);
nstars=40;
scatter(a0,2*rand(nstars,1)-1,rand(nstars,1),5*rand(nstars,1),'w','filled');
hold(a0,'on');
phig=linspace(0,2*pi);
xe=.1;
ye=.5;
re=.1;
patch(a0,'XData',xe+re*cos(phig),'YData',ye+re*sin(phig),'EdgeColor','none','FaceColor',[.2 .3 .7]);
africa=[35 -6;21 -17;14 -17;4 -8;5 7;-11 14;-18 11;-35 19;-33 27;-24 35;-20 34;-15 41;-4 40;12 51;11 44;28 34;31 32;33 20;31 20;33 11;37 11;35 -2]*pi./180;
eurasia=[28 35;13 44;18 56;22 60;26 56;24 52;29 48;29 90;40 90;50 90;60 90;74 90;64 36;66 33;66 40;68 40;71 26;62 5;58 8;52 16;66 21;65 25;60 22;55 20;54 6;48 -5;44 -2;44 -9;37 -9;36 -6;37 -1;43 3;44 9;40 16;38 16;38 12;37 15;38 15;39 17;40 16;40 19;45 12;46 14;42 20;37 22;38 24;40 23;40 26;36 28;37 36]*pi./180;
southamerica=[-6 -35;-1 -49;5 -52;10 -63;12 -71;8 -78;-4 -81;-18 -71;-50 -76;-55 -66;-49 -69;-21 -41;-13 -40]*pi./180;
patch(a0,'XData',xe+re.*cos(africa(:,1)).*sin(africa(:,2)),'YData',ye+re.*sin(africa(:,1)),'FaceColor',[.7 .7 .5],'EdgeColor','none');
patch(a0,'XData',xe+re.*cos(eurasia(:,1)).*sin(eurasia(:,2)),'YData',ye+re.*sin(eurasia(:,1)),'FaceColor',[.7 .7 .5],'EdgeColor','none');
patch(a0,'XData',xe+re.*cos(southamerica(:,1)).*sin(southamerica(:,2)),'YData',ye+re.*sin(southamerica(:,1)),'FaceColor',[.7 .7 .5],'EdgeColor','none');
hold(a0,'off');
axis(a0,'off');
axis(a0,'equal');
ylim(a0,[0 1]);
h_dust=plot(a1,nan,nan,'.','Color',[.7 .7 .7],'MarkerSize',2);
hold(a1,'on');
area(a1,xg,yg,'FaceColor',[.4 .5 .6],'BaseValue',-inf);
hs=scatter(a1,1:3,1:3,20,[0 0 0],'filled');
hp=patch(a1,'XData',1:3,'YData',[2 3 1],'FaceColor',[.5 .5 .5],'EdgeColor','none');
hl=plot(a1,1:3,1:3,'LineWidth',2,'Color',[.7 .7 .7]);
hold(a1,'off');
axis(a1,'equal');
axis(a1,'off');
hte=uicontrol(f,'Style','text','Units','normalized','Position',[0 0 1 .1],'FontName','FixedWidth','ForegroundColor','w','BackgroundColor','k');
hs2=scatter(a2,1:3,1:3,20,[0 0 0],'filled');
hold(a2,'on');
hp2=patch(a2,'XData',1:3,'YData',[2 3 1],'FaceColor',[.5 .5 .5],'EdgeColor','none');
hl2=plot(a2,1:3,1:3,'LineWidth',2,'Color',[.7 .7 .7]);
hold(a2,'off');
axis(a2,'equal');
a2.Color=[.2 .2 .2];
a2.XTick=[];
a2.YTick=[];
fsurf=griddedInterpolant(xg,yg);
fgrad=griddedInterpolant(xg,gradient(yg,xg));
complete=0;
colormap(f,'autumn');
t0=tic;
t2=toc(t0);
while ishandle(a1)
    
    t1=t2;
    t2=toc(t0);
    dt=t2-t1;
    
    u=[r;v;theta;omega;thrust;m_fuel;diss];
    
    % maximum step size 5ms
    N=ceil(dt/5e-3);
    for j=1:N
        % dormand-prince method
        du1=der(u);
        du2=der(u+.2*dt/N*du1);
        du3=der(u+dt/N*(3/40*du1+9/40*du2));
        du4=der(u+dt/N*(44/45*du1-56/15*du2+32/9*du3));
        du5=der(u+dt/N*(19372/6561*du1-25360/2187*du2+64448/6551*du3-212/729*du4));
        du6=der(u+dt/N*(9017/3168*du1-355/33*du2+46732/5247*du3+49/176*du4-5103/18656*du5));
        du=35/384*du1+500/1113*du3+125/192*du4-2187/6784*du5+11/84*du6;
        u=u+du*dt/N;
    end
    
    r=u(1:2);
    v=u(3:4);
    theta=u(5);
    omega=u(6);
    m_fuel=max(0,u(8));
    diss=u(9);
    thrust=min(nom_thrust,max(0,u(7)))*(m_fuel>0);
    r(1)=mod(r(1)+1e3,2e3)-1e3;
    
    alt=r(2)-fsurf(r(1));
    
    x_proj=r(1)-cot(theta)*alt;
    y_proj=fsurf(x_proj);
    
    S=max(20,abs(r(2)));
    
    M=[-sin(theta) cos(theta);cos(theta) sin(theta)];
    f_rot=(M*[fx fy]')';
    g_rot=(M*[gx gy]')';
    g_rot([3 10],2)=max(r(2)+g_rot([3 10],2),fsurf(r(1)+g_rot([3 10],1)))-r(2);
    hp.XData=r(1)+f_rot(:,1);
    hp.YData=r(2)+f_rot(:,2);
    hl.XData=r(1)+g_rot(:,1);
    hl.YData=r(2)+g_rot(:,2);
    
    d=(-2.5+1.5*[-2;-1;0]*thrust/(2e3+thrust));
    hs.XData=r(1)+cos(theta)*d+[.1;.1;0].*randn(3,1);
    hs.YData=r(2)+sin(theta)*d+[.1;.1;0].*randn(3,1);
    hs.SizeData=max(1,thrust)*[.06;.03;.01].*20./S;
    hs.CData=[0;.5;1];
    hs2.XData=cos(theta)*d+[.1;.1;0].*randn(3,1);
    hs2.YData=sin(theta)*d+[.1;.1;0].*randn(3,1);
    hs2.SizeData=max(1,thrust)*[.06;.03;.01];
    hs2.CData=[0;.5;1];
    
    if (r(1)-x_proj)*cos(theta)+(r(2)-y_proj)*sin(theta)>0
        n_dust=floor(.3*thrust/(1+alt.^2));
        x_dust=x_proj+.3*randn(n_dust,1)*max(4,min(40,alt));
        phi_dust=2*atan(fgrad(x_dust))-theta+pi+.3*randn(n_dust,1);
        r_dust=1e-3*thrust*randn(n_dust,1).^2;
        h_dust.XData=x_dust+cos(phi_dust).*r_dust;
        h_dust.YData=fsurf(x_dust)+sin(phi_dust).*r_dust;
    else
        h_dust.XData=nan;
        h_dust.YData=nan;
    end
    
    if any(r(2)+f_rot(:,2)<=fsurf(r(1)+f_rot(:,1)))
        hold(a1,'on');
        n=10;
        scatter(a1,r(1)+randn(n,1),r(2)+randn(n,1),exp(5+randn(n,1)),rand(n,1),'filled','MarkerFaceAlpha',.5);
        hold(a1,'off');
        errordlg('MISSION FAILED !','LUNAR LANDER');
        break
    end
    
    if thrust==0 && sum(v.^2)<1e-2
        complete=complete+dt;
    else
        complete=0;
    end
    
    if complete>2
        comm=inputdlg(sprintf('MISSION COMPLETE !\nResidual kinetic energy was %1.0f J.\nRemaining fuel was %1.0f kg.\n\nPress OK to sign in to the high score list (with optional comment).\nPress CANCEL to skip entering your mission to the high score list.\n',diss,m_fuel),'LUNAR LANDER');
        if ~isempty(comm)
            fid=fopen(hslist,'a');
            if fid>0
                fprintf(fid,'%s %1.0f %1.0f "%s"\n',getenv('username'),diss,m_fuel,comm{:});
                fclose(fid);
            end
            showhs;
        end
        break
    end
    
    if r(2)>20
        hp2.XData=f_rot(:,1);
        hp2.YData=f_rot(:,2);
        hl2.XData=g_rot(:,1);
        hl2.YData=g_rot(:,2);
        axis(a2,'on');
    else
        hp2.XData=nan;
        hp2.YData=nan;
        hl2.XData=nan;
        hl2.YData=nan;
        axis(a2,'off');
    end
    
    hte.String=sprintf('FUEL % 5.0fkg   THRUST % 5.0fN   ALTITUDE % 5.0fm\nHORIZONTAL SPEED %+ 4.1fm/s   VERTICAL SPEED %+ 4.1fm/s   ROTATIONAL SPEED %+ 2.3f/s',m_fuel,thrust,alt,v(1),v(2),omega);
    
    xlim(a1,r(1)+[-S S]);
    ylim(a1,[0 1.1*S+5]);
    xlim(a2,[-10 10]);
    ylim(a2,[-10 10]);
    
    drawnow;
    
end
end
