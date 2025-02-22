clear all; clc; clf;

%% SET INITIAL DATA

% Set times for trajectory
t1 = 0;                % initial time for trajectories, from Jan.1,2017, in days
t2 = 2000;             % final time, in days
dt = 1/(24*60);        % integration time interval, in days

% Set radius in au
au = 1.495978707e11;        % one astronomical unit, in m
Rt=6371000/au;              % Earth's radius
Rm=1737000/au;              % Moon's radius
Rs=695700000/au;            % Sun's radius
R=[Rs Rt Rm];

%% FIND TRAJECTORIES AND DISTANCES

% Find trajectories, we call function solarSystem
[r,body] = solarSystem( t1, t2, dt ); 
[nd,nb,nt] = size(r);
rt=r(:,2,:);   %Earth's trayectory
rm=r(:,3,:);   %Moon's trayectory
rs=r(:,1,:);   %Sun's trayectory

% Find distances between bodies
rst=rt-rs; 
rsm=rm-rs;
rmt=rt-rm;
rsta=zeros(nt,1);
rsma=rsta;
rmta=rsta;
for i=1:nt
    rsta(i,:)=norm(rst(:,:,i));  %Sun-Earth, in au
    rsma(i,:)=norm(rsm(:,:,i));  %Sun-Moon, in au
    rmta(i,:)=norm(rmt(:,:,i));  %Moon-Earth, in au
end

%% FIND RADIUS OF THE SHADOW
%Angulo entre rst y rsm
thetas=zeros(nt,1);
xes=zeros(nt,1);
Ros=zeros(nt,1);
thetal=thetas;
xel=xes;
Rol=Ros;

for i=1:nt 
    thetas(i,:)=acos((dot(rst(:,:,i),rsm(:,:,i)))/(rsta(i,:)*rsma(i,:)));  %Angle between Sun-Earth and Sun-Moon
    xes(i,:)=sin(thetas(i,:))*rsta(i,:);                                   %Distance Shadow-Earth
    Ros(i,:)=Rm-((rmta(i,:)*(Rs))/rsta(i,:));                              %Radius of the shadow
    Ros(i,:)=abs(Ros(i,:));
  
    thetal(i,:)=acos((dot(rst(:,:,i),rsm(:,:,i)))/(rsta(i,:)*rsma(i,:)));  %Angle between Sun-Earth and Sun-Moon
    xel(i,:)=sin(thetal(i,:))*rsma(i,:);                                   %Distance Shadow-Moon
    Rol(i,:)=Rt-((rmta(i,:)*(Rs))/rsma(i,:));                              %Radius of the shadow
    Rol(i,:)=abs(Rol(i,:));
end

%% FIND WHEN SOLAR AND LUNAR ECLIPSES OCCUR

se=[0]; 
le=[0];

for i=1:nt
    if xes(i,:)<Ros(i,:)+Rt && rsta(i,:)>rsma(i,:)        % Condition for a solar eclipse
    se(i,:)=i;                                            % Times when a solar eclipse occur
    elseif xel(i,:)<Rol(i,:)+Rm && rsta(i,:)<rsma(i,:)    % Condition for a lunar eclipse
    le(i,:)=i;                                            % Times when a lunar eclipse occur
    end      
end

% Deploy error if no eclipse occur
if numel(se)==0 && numel(le)==0
    error('No solar or lunar eclipse during the time given')
end

% Find when a eclipse begin and end
se=se(se~=0);
A=circshift(se,[1,0])+1;
A1=find(se~=A);            % What times solar eclipses begin

le=le(le~=0);
B=circshift(le,[1,0])+1;
B1=find(le~=B);            % What times lunar eclipses begin

%% DISPLAY SOLAR ECLIPSE INFORMATION 

% Find what hour and date is the initial time given t1
t0=t1*24*60;
year0=fix(t0/(60*24*365.25))+2017;
rest0=mod(t0,(60*24*365.25));
month0=fix(rest0/(60*24*30.4375))+1;
rest02=mod(rest0,(60*24*30.4375));
day0=fix(rest02/(60*24))+1;
rest03=mod(rest02,(60*24));
hour0=fix(rest03/60);
minute0=mod(rest03,60);

% Find exact time and date of the solar eclipse
if numel(se)>0
fprintf(['There will be solar eclipses on the dates:\n'])
elseif numel(se)==0
fprintf(['There will be no solar eclipses\n'])
end

for j=1:numel(A1)
    if j<numel(A1)
    time1=se(A1(j,:),:);
    time2=se(A1(j+1,:)-1,:);
    
    % Time and date of the beginning of the eclipse
    year1=fix(time1/(60*24*365.25));
    rest1=mod(time1,(60*24*365.25));
    month1=fix(rest1/(60*24*30.4375));
    rest2=mod(rest1,(60*24*30.4375));
    day1=fix(rest2/(60*24));
    rest3=mod(rest2,(60*24));
    hour1=fix(rest3/60);
    minute1=mod(rest3,60);

    % We set the date and hour depending on t1 given
    YEAR1=year1+year0;
    month1=month1+month0;
    day1=day1+day0;
    hour1=hour1+hour0;
    minute1=minute1+minute0;

    %We modify date and time if they aren't correct
    if minute1>59
    minute1=minute1-60;
    hour1=hour1+1;
    end

    if hour1>23
    hour1=hour1-24;
    day1=day1+1;
    end

    if day1>30
    day1=day1-30;
    month1=month1+1;
    end

    if month1>12
    month1= month1-12;
    YEAR1=YEAR1+1;
    end

    % Time and date of the end of the eclipse
    year2=fix(time2/(60*24*365.25));
    rest12=mod(time2,(60*24*365.25));
    month2=fix(rest12/(60*24*30.4375));
    rest22=mod(rest12,(60*24*30.4375));
    day2=fix(rest22/(60*24));
    rest32=mod(rest22,(60*24));
    hour2=fix(rest32/60);
    minute2=mod(rest32,60);

    % We set the date and hour depending on t1 given
    YEAR2=year2+year0;
    month2=month2+month0;
    day2=day2+day0;
    hour2=hour2+hour0;
    minute2=minute2+minute0;

    %We modify date and time if they aren't correct
    if minute2>59
    minute2=minute2-60;
    hour2=hour2+1;
    end

    if hour2>23
    hour2=hour2-24;
    day2=day2+1;
    end

    if day2>30
    day2=day2-30;
    month2=month2+1;
    end

    if month2>12
    month2= month2-12;
    YEAR2=YEAR2+1;
    end

    % Display the information 
    
    % We add zeros to the hour if it doesn't have two digits
    if minute1<10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    end
    
    % We modify it to work for the last eclipse
    elseif j==numel(A1)  
    time1=se(A1(j,:));
    time2=se(end);
          
    year1=fix(time1/(60*24*365.25));
    rest1=mod(time1,(60*24*365.25));
    month1=fix(rest1/(60*24*30.4375));
    rest2=mod(rest1,(60*24*30.4375));
    day1=fix(rest2/(60*24));
    rest3=mod(rest2,(60*24));
    hour1=fix(rest3/60);
    minute1=mod(rest3,60);

    % We set the date and hour depending on t1 given
    YEAR1=year1+year0;
    month1=month1+month0;
    day1=day1+day0;
    hour1=hour1+hour0;
    minute1=minute1+minute0;

    %We modify date and time if they aren't correct
    if minute1>59
    minute1=minute1-60;
    hour1=hour1+1;
    end

    if hour1>23
    hour1=hour1-24;
    day1=day1+1;
    end

    if day1>30
    day1=day1-30;
    month1=month1+1;
    end

    if month1>12
    month1= month1-12;
    YEAR1=YEAR1+1;
    end

    year2=fix(time2/(60*24*365.25));
    rest12=mod(time2,(60*24*365.25));
    month2=fix(rest12/(60*24*30.4375));
    rest22=mod(rest12,(60*24*30.4375));
    day2=fix(rest22/(60*24));
    rest32=mod(rest22,(60*24));
    hour2=fix(rest32/60);
    minute2=mod(rest32,60);

    % We set the date and hour depending on t1 given
    YEAR2=year2+year0;
    month2=month2+month0;
    day2=day2+day0;
    hour2=hour2+hour0;
    minute2=minute2+minute0;

    %We modify date and time if they aren't correct
    if minute2>59
    minute2=minute2-60;
    hour2=hour2+1;
    end

    if hour2>23
    hour2=hour2-24;
    day2=day2+1;
    end

    if day2>30
    day2=day2-30;
    month2=month2+1;
    end

    if month2>12
    month2= month2-12;
    YEAR2=YEAR2+1;
    end

    % We add zeros to the hour if it doesn't have two digits
    if minute1<10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    end
    
    end
end

%% DISPLAY LUNAR ECLIPSE INFORMATION 

% Find exact time and date of the lunar eclipse
if numel(le)>0
fprintf(['\nThere will be lunar eclipses on the dates:\n'])
elseif numel(le)==0
fprintf(['\nThere will be no lunar eclipses.\n'])
end

for j=1:numel(B1)
    if j<numel(B1)
    time1=le(B1(j,:));
    time2=le(B1(j+1,:)-1);
    
    % Time and date of the beginning of the eclipse
    year1=fix(time1/(60*24*365.25));
    rest1=mod(time1,(60*24*365.25));
    month1=fix(rest1/(60*24*30.4375));
    rest2=mod(rest1,(60*24*30.4375));
    day1=fix(rest2/(60*24));
    rest3=mod(rest2,(60*24));
    hour1=fix(rest3/60);
    minute1=mod(rest3,60);

    % We set the date and hour depending on t1 given
    YEAR1=year1+year0;
    month1=month1+month0;
    day1=day1+day0;
    hour1=hour1+hour0;
    minute1=minute1+minute0;

    %We modify date and time if they aren't correct
    if minute1>59
    minute1=minute1-60;
    hour1=hour1+1;
    end

    if hour1>23
    hour1=hour1-24;
    day1=day1+1;
    end

    if day1>30
    day1=day1-30;
    month1=month1+1;
    end

    if month1>12
    month1= month1-12;
    YEAR1=YEAR1+1;
    end
    
    % Time and date of the end of the eclipse
    year2=fix(time2/(60*24*365.25));
    rest12=mod(time2,(60*24*365.25));
    month2=fix(rest12/(60*24*30.4375));
    rest22=mod(rest12,(60*24*30.4375));
    day2=fix(rest22/(60*24));
    rest32=mod(rest22,(60*24));
    hour2=fix(rest32/60);
    minute2=mod(rest32,60);

    % We set the date and hour depending on t1 given
    YEAR2=year2+year0;
    month2=month2+month0;
    day2=day2+day0;
    hour2=hour2+hour0;
    minute2=minute2+minute0;

    %We modify date and time if they aren't correct
    if minute2>59
    minute2=minute2-60;
    hour2=hour2+1;
    end

    if hour2>23
    hour2=hour2-24;
    day2=day2+1;
    end

    if day2>30
    day2=day2-30;
    month2=month2+1;
    end

    if month2>12
    month2= month2-12;
    YEAR2=YEAR2+1;
    end
    
    % Display the information
    % We add zeros to the hour if it doesn't have two digits
    if minute1<10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    end
    
    % We modify it to work for the last eclipse
    elseif j==numel(B1) 
    time1=le(B1(j,:));
    time2=le(end);
          
    year1=fix(time1/(60*24*365.25));
    rest1=mod(time1,(60*24*365.25));
    month1=fix(rest1/(60*24*30.4375));
    rest2=mod(rest1,(60*24*30.4375));
    day1=fix(rest2/(60*24));
    rest3=mod(rest2,(60*24));
    hour1=fix(rest3/60);
    minute1=mod(rest3,60);

    % We set the date and hour depending on t1 given
    YEAR1=year1+year0;
    month1=month1+month0;
    day1=day1+day0;
    hour1=hour1+hour0;
    minute1=minute1+minute0;

    %We modify date and time if they aren't correct
    if minute1>59
    minute1=minute1-60;
    hour1=hour1+1;
    end

    if hour1>23
    hour1=hour1-24;
    day1=day1+1;
    end

    if day1>30
    day1=day1-30;
    month1=month1+1;
    end

    if month1>12
    month1= month1-12;
    YEAR1=YEAR1+1;
    end

    year2=fix(time2/(60*24*365.25));
    rest12=mod(time2,(60*24*365.25));
    month2=fix(rest12/(60*24*30.4375));
    rest22=mod(rest12,(60*24*30.4375));
    day2=fix(rest22/(60*24));
    rest32=mod(rest22,(60*24));
    hour2=fix(rest32/60);
    minute2=mod(rest32,60);

    % We set the date and hour depending on t1 given
    YEAR2=year2+year0;
    month2=month2+month0;
    day2=day2+day0;
    hour2=hour2+hour0;
    minute2=minute2+minute0;

    %We modify date and time if they aren't correct
    if minute2>59
    minute2=minute2-60;
    hour2=hour2+1;
    end

    if hour2>23
    hour2=hour2-24;
    day2=day2+1;
    end

    if day2>30
    day2=day2-30;
    month2=month2+1;
    end

    if month2>12
    month2= month2-12;
    YEAR2=YEAR2+1;
    end

    % Display the information
    % We add zeros to the hour if it doesn't have two digits
    if minute1<10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2<10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2>=10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , %2.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1<10 && hour1>=10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , %2.0f:0%1.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    elseif minute1>=10 && hour1<10 && minute2>=10 && hour2<10
    fprintf(['%4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f - %4.0f/%2.0f/%2.0f , 0%1.0f:%2.0f\n'], ...
     day1,month1,YEAR1,hour1,minute1,day2,month2,YEAR2,hour2,minute2)
    end
    end
end



%% CREATE SUN, EARTH AND MOON

% Sun's sphere and trajectory 
xs(:)=r(1,1,:);
ys(:)=r(2,1,:);
zs(:)=r(3,1,:);

[X,Y,Z]=sphere(50);
XS=X*R(1,1);
YS=Y*R(1,1);
ZS=Z*R(1,1);
    
% Moon's sphere and trajectory
xl(:)=r(1,3,:);
yl(:)=r(2,3,:);
zl(:)=r(3,3,:);

[X,Y,Z]=sphere(25);
XL=X*R(1,3);
YL=Y*R(1,3);
ZL=Z*R(1,3);
    
% Earth's sphere and trajectory 
xt(:)=r(1,2,:);
yt(:)=r(2,2,:);
zt(:)=r(3,2,:);
    
[xx,yy,zz]=sphere(25);
XT2=xx*R(1,2);
YT2=yy*R(1,2);
ZT2=zz*R(1,2);

%% PLOT THE MOVIE OF THE ECLIPSE 

% Find a good position to plot the eclipse
medx=((xt+xl)/2);
medy=((yt+yl)/2);
medz=((zt+zl)/2);

% Time interval to plot the eclipse
if numel(se)==0
timeeclipse1=(le(end)-340);                            % Initial time of the plot
timeeclipse2=(le(B1(end))+340);                        % Final time of the plot
elseif numel(se)>0    
timeeclipse1=(se(end)-340);                            % Initial time of the plot
timeeclipse2=(se(A1(end))+340);                        % Final time of the plot
end
medtime=round((timeeclipse1+timeeclipse2)/2);          % Time to center the plot
timeeclipse= (timeeclipse1-1000):(timeeclipse2+1000);  % Time interval to plot the trajectory

% Create a movie of the eclipse
for it = timeeclipse1:timeeclipse2
figure(1)
set(figure(1),'Name','Eclipse','NumberTitle','off')
% Movie in 3D
subplot(1,2,1)

% Plot trajectories
plot3(xt(timeeclipse),yt(timeeclipse),zt(timeeclipse),'Color',[0.4,0.6,0.2])
hold on 
plot3(xl(timeeclipse),yl(timeeclipse),zl(timeeclipse),'k')

grid on
axis equal

% Plot sunlight and spheres
plot3([xs(it)  xl(it)],[ys(it)  yl(it)],[zs(it)  zl(it)],'Color',[1,0.7,0]);
plot3([xs(it)  xt(it)],[ys(it)  yt(it)],[zs(it)  zt(it)],'r');
plot3(XT2+xt(it),YT2+yt(it),ZT2+zt(it),'Color',[0.4,0.5,0.8])
plot3(XL+xl(it),YL+yl(it),ZL+zl(it),'k');

legend('Earth´s trajectory','Moon´s trajectory', ...
       'Sunlight of the Moon','Sunlight of the Earth','Location','northwest')
set(gca,'fontsize',10)

% Set limits
xlim([(medx(medtime)-70*4.26*10^-5)    (medx(medtime)+70*4.26*10^-5)])
ylim([(medy(medtime)-70*4.26*10^-5)    (medy(medtime)+70*4.26*10^-5)])
zlim([(medz(medtime)-40*4.26*10^-5)    (medz(medtime)+40*4.26*10^-5)])

xlabel('x (au)'), ylabel('y (au)'), zlabel('z (au)')
title('3D plot of the eclipse')
set(gca,'Color',[0.7,0.9,1])
hold off

% Movie in 2D
subplot(1,2,2)

%Plot trajectories
plot3(xt(timeeclipse),yt(timeeclipse),zt(timeeclipse),'Color',[0.4,0.6,0.2])
hold on 
plot3(xl(timeeclipse),yl(timeeclipse),zl(timeeclipse),'k')
view(2)
grid on
axis equal

%Plot sunlight
plot3([xs(it)  xl(it)],[ys(it)  yl(it)],[zs(it)  zl(it)],'Color',[1,0.7,0])
plot3([xs(it)  xt(it)],[ys(it)  yt(it)],[zs(it)  zt(it)],'r')
plot3(XT2+xt(it),YT2+yt(it),ZT2+zt(it),'Color',[0.4,0.5,0.8])
plot3(XL+xl(it),YL+yl(it),ZL+zl(it),'k')

%Set limits
xlim([(medx(medtime)-70*4.26*10^-5)    (medx(medtime)+70*4.26*10^-5)])
ylim([(medy(medtime)-70*4.26*10^-5)    (medy(medtime)+70*4.26*10^-5)])
zlim([(medz(medtime)-40*4.26*10^-5)    (medz(medtime)+40*4.26*10^-5)])

xlabel('x (au)'), ylabel('y (au)'), zlabel('z (au)')
title('2D plot of the eclipse')

%We inlcude a watch
year=fix(it/(60*24*365.25));
REST1=mod(it,(60*24*365.25));
month=fix(REST1/(60*24*30.4375));
REST2=mod(REST1,(60*24*30.4375));
day=fix(REST2/(60*24));
REST3=mod(REST2,(60*24));
hour=fix(REST3/60);
minute=mod(REST3,60);

% We set the date and hour depending on t1 given
year=year+year0;
month=month+month0;
day=day+day0;
hour=hour+hour0;
minute=minute+minute0;

%We modify date and time if they aren't correct
if minute>59
minute=minute-60;
hour=hour+1;
end

if hour>23
hour=hour-24;
day=day+1;
end

if day>30
day=day-30;
month=month+1;
end

if month>12
month= month-12;
year=year+1;
end

%We display the watch as a title, adding zeros if the hour is not correct
if numel(se)>0
    if hour<10 && minute<10
    txt1 = ['Solar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = ['0' num2str(hour) ': 0' num2str(minute)]; 
    elseif hour>=10 && minute<10
    txt1 = ['Solar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = [num2str(hour) ': 0' num2str(minute)];
    elseif hour<10 && minute>=10
    txt1 = ['Solar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = ['0' num2str(hour) ':' num2str(minute)];
    elseif hour>=10 && minute>=10
    txt1 = ['Solar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = [num2str(hour) ':' num2str(minute)];
    end
elseif numel(se)==0
    if hour<10 && minute<10
    txt1 = ['Lunar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = [ '0' num2str(hour) ': 0' num2str(minute)];
    elseif hour>=10 && minute<10
    txt1 = ['Lunar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = [num2str(hour) ': 0' num2str(minute)];
    elseif hour<10 && minute>=10
    txt1 = ['Solar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = ['0' num2str(hour) ':' num2str(minute)];
    elseif hour>=10 && minute>=10
    txt1 = ['Lunar eclipse on the day ' num2str(day) '/' num2str(month) '/' num2str(year) ];
    txt2 = [num2str(hour) ':' num2str(minute)];
    end
end

txt= [txt1 ' , ' txt2];
sgtitle(txt)
set(gca,'Color',[0.7,0.9,1])

hold off

pause(0.000000000000000001)
end

%% PLOT THE SOLAR SYSTEM

figure(2)
plot3(xs,ys,zs,'Color',[0.3,0.4,0.5])
hold on
plot3(xt,yt,zt,'Color','g')
plot3(xl,yl,zl,'k')
plot3(XS+xs(se(1)),YS+ys(se(1)),ZS+zs(se(1)),'Color',[0.3,0.4,0.5])

xlabel('x (au)'), ylabel('y (au)'), zlabel('z (au)')
grid on
axis equal
legend('Sun´s trajectory', 'Earth´s trajectory','Moon´s trajectory')
title('Solar system in scale')
set(gca,'Color',[0.7,0.9,1])
set(figure(2),'Name','Solar system','NumberTitle','off')

hold off


%% ADDITIONAL INFORMATION
nmin=find(xes==min(xes));
figure(3)
plot([0   0],[0  Rt],[1  1],[0  Ros(nmin)],[2  2],[0  xes(nmin)],'linewidth',2)
axis([-0.5  2.5  0  7*10^-4])
ylabel('Distance (au)'), grid on;
legend('Radius of Earth','Radius of Shadow','Distance Earth-Shadow')


figure(4)
t=(t1:dt:t2);
plot(t,xes,t,(rsma*(2.5*10^-3)));
xlabel('Time (days)'), ylabel('Distance (au)'), grid on;
legend('Distance Earth-shadow','Distance Sun-Moon (2.5*10^-3)')
set(figure(4),'Name','Additional data','NumberTitle','off')



