clear all
global nx ny nz k_a k_b aerfa a b c lambda theta_num front_constant gama beta1;
%set parameters
prompt={'nx','ny','nz','lambda(in A)'};
dlg_title='Set Parameters for PXRD Simulation';
num_lines=4;
def={'100','100','1','1.54'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer=inputdlg(prompt,dlg_title,num_lines,def,options);
nx=str2num(answer{1});
ny=str2num(answer{2});
nz=str2num(answer{3});
lambda=str2num(answer{4});
prompt={'interval of theta(in ¡ã)','maximum of theta(in ¡ã)','minimum of theta(in ¡ã)'};
dlg_title='Set Parameters for PXRD Simulation';
num_lines=4;
def={'0.05','5','1.5'};
answer=inputdlg(prompt,dlg_title,num_lines,def);
deta_theta=str2num(answer{1}).*pi/180;
theta_max=str2num(answer{2}).*pi/180;
theta_min=str2num(answer{3}).*pi/180;
theta_length=0.02/180*pi;
theta_fine_num=11;
prompt={'k_a','k_b'};
dlg_title='Set Parameters for PXRD Simulation';
num_lines=2;
def={'0','0'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answer=inputdlg(prompt,dlg_title,num_lines,def,options);
k_a=str2num(answer{1});
k_b=str2num(answer{2});

%read in cif structure
addpath(pwd);
[filename, pathname]=uigetfile('*.cif','select the cif file');
samplename=filename(1:end-4);
filename=strcat(pathname,filename);
cd(pathname);
fid=fopen(filename,'r');

i = 0;
targetline1 = '_symmetry_equiv_pos_as_xyz';
tag = 1;
%find the line of symmetry operation
while tag
   i=i+1;
   readline = fgetl(fid);
   if strcmp(readline,targetline1)==1
       tag = 0;
   end
end

 readline = fgetl(fid);
 count = 1;
 %save symmetry operation
while isempty(strfind(readline,'_cell_length_a'))
    [str1, str2, str3] = strread(readline,'%s %s %s','delimiter',',');
    value(count,:) = [str1 str2 str3];
    count = count + 1;
     readline = fgetl(fid);
end
clear i str1 str2 str3 count;

%read in cell parameters: a, b and c
i =1;
while i <= 3
    [str_length,value_length] = strread(readline,'%s %f','delimiter',' ');
    para_len(1,i) = value_length/10;
    i = i + 1;
    readline = fgetl(fid);
end
clear i str_length value_length;

%read in cell angles£ºAlpha,Beta,Gamma
i =1;
while i <= 3
    [str_angle,value_angle] = strread(readline,'%s %f','delimiter',' ')
    para_anle(1,i) = 180-value_angle;
    i = i + 1;
    readline = fgetl(fid);
end
clear i str_angle value_angle;

%find atom coordinations
targetline2 = '_atom_site_occupancy';
tag = 1;
while tag
   readline = fgetl(fid);
   if strcmp(readline,targetline2)==1
       tag = 0;
   end
end
clear tag;
%read in elements
flag = fopen('Atom.txt');
line = fgetl(flag);
i =1;
while i< 99
    
    Atom{i,1} = strread(line,'%s');
    line = fgetl(flag);
    i = i + 1;
end
    
%read in atom coordinations
numcount = 1;
readline = fgetl(fid);
i=1;
while ~strcmp(readline,'loop_')
 [at1,at2,X,Y,Z,uiso,adp,occup] = strread(readline,'%s%s%f%f%f%f%s%f','delimiter',' ');
 %find element
 if at2{1,1} == 'H'
     coord(numcount,1) = 1;
 elseif at2{1,1} == 'C'
      coord(numcount,1) = 6;
 elseif at2{1,1} == 'N'
      coord(numcount,1) = 7;
 elseif at2{1,1} == 'O'  
      coord(numcount,1) = 8;
 elseif at2{1,1} == 'Al'  
      coord(numcount,1) = 13; 
 elseif at2{1,1} == 'Cu'  
      coord(numcount,1) = 29;
 elseif at2{1,1} == 'Zn'
      coord(numcount,1) = 30;
 elseif at2{1,1} == 'Zr'
      coord(numcount,1) = 40;
 elseif at2{1,1} == 'Hf'
      coord(numcount,1) = 72;
 %find elements in perodic table
 else
     while i < 99
         if Atom{i,1}{1,1} == at2{1,1}
             coord(numcount,1) = i;
             break;
         end
       i = i+1;
     end
 end                
 coord(numcount,2:4) = [X Y Z];
 numcount = numcount + 1;
 readline = fgetl(fid);
end
clear at1 at2 X Y Z uiso adp occup;
fclose(fid);
clear fid targetline1 targetline2 readline;

%Coordinations are saved in Coord
count = 1;
syms x y z;
i = 1;
while i <= size(coord,1)
   j = 1;
   while j <= size(value,1)
    Coord(count,1) = coord(i,1);
    Coord(count,2:4) = subs(value(j,:),[x y z],[coord(i,2) coord(i,3) coord(i,4)]);
    Coord(count,2:4) = Coord(count,2:4) - floor(Coord(count,2:4));
    j = j + 1;
    count = count + 1;
   end
   i = i + 1;
end
 Coord = unique(Coord,'rows');
clear i j x y z;
%set core numbers
%matlabpool local 6

aerfa=para_anle(1,3)/180*pi;%angle between axis a and b
gama=para_anle(1,1)/180*pi; %angle between axis b and c
beta1=para_anle(1,2)/180*pi;%angle between axis a and c
a = para_len(1,1)*10;
b = para_len(1,2)*10;
c = para_len(1,3)*10;

%calculation of axis c_2, which is vertical to the plane a and b
A=0;
B=pi/2;
SIA=(sqrt(1-A.^2));
c_2=c*(A.*sin(B).^2.*cos(beta1).*sin(aerfa).^2 + SIA.*sin(B).*sin(aerfa).*(A.^4.*cos(B).^2.*sin(aerfa).^2 - SIA.^2.*cos(gama).^2 - A.^2.*cos(gama).^2 + A.^2.*sin(B).^2.*sin(aerfa).^2 + cos(B).^2.*SIA.^4.*sin(aerfa).^2 - A.^2.*cos(aerfa).^2.*cos(beta1).^2 + SIA.^2.*sin(B).^2.*sin(aerfa).^2 - SIA.^2.*cos(aerfa).^2.*cos(beta1).^2 - sin(B).^2.*cos(beta1).^2.*sin(aerfa).^2 + 2.*A.^2.*cos(B).^2.*SIA.^2.*sin(aerfa).^2 + 2.*A.^2.*cos(aerfa).*cos(beta1).*cos(gama) - A.^2.*cos(B).^2.*cos(beta1).^2.*sin(aerfa).^2 + 2.*SIA.^2.*cos(aerfa).*cos(beta1).*cos(gama) - cos(B).^2.*SIA.^2.*cos(beta1).^2.*sin(aerfa).^2).^(1./2) + A.^3.*cos(B).^2.*cos(beta1).*sin(aerfa).^2 - cos(B).*SIA.^3.*cos(gama).*sin(aerfa) + A.*cos(B).^2.*SIA.^2.*cos(beta1).*sin(aerfa).^2 - A.^2.*cos(B).*SIA.*cos(gama).*sin(aerfa) + cos(B).*SIA.^3.*cos(aerfa).*cos(beta1).*sin(aerfa) + A.^2.*cos(B).*SIA.*cos(aerfa).*cos(beta1).*sin(aerfa))./(A.^4.*cos(B).^2.*sin(aerfa).^2 + A.^2.*sin(B).^2.*sin(aerfa).^2 + cos(B).^2.*SIA.^4.*sin(aerfa).^2 + SIA.^2.*sin(B).^2.*sin(aerfa).^2 + 2.*A.^2.*cos(B).^2.*SIA.^2.*sin(aerfa).^2);
clearvars A B SIA;

%find possible 3D peak positions 
front_constant=4*pi*1i/lambda;
theta_extra=theta_min:deta_theta:theta_max;
hmax=fix(2*sin(theta_max)*a/lambda);
kmax=fix(2*sin(theta_max)*b/lambda);
lmax=fix(2*sin(theta_max)*c/lambda);
h=-hmax:hmax;
k=-kmax:kmax;
l=-lmax:lmax;
all_theta=[];
hkl=[];
 for h=-hmax:hmax
     for k=-kmax:kmax
         for l=-lmax:lmax
             matrix_1=(h/a)*det([h/a,cos(aerfa),cos(beta1);k/b,1,cos(gama);l/c,cos(gama),1]);
             matrix_2=(k/b)*det([1,h/a,cos(beta1);cos(aerfa),k/b,cos(gama);cos(beta1),l/c,1]);
             matrix_3=(l/c)*det([1,cos(aerfa),h/a;cos(aerfa),1,k/b;cos(beta1),cos(gama),l/c]);
             matrix_4=det([1,cos(aerfa),cos(beta1);cos(aerfa),1,cos(gama);cos(beta1),cos(gama),1]);
             d=((matrix_1+matrix_2+matrix_3)/matrix_4)^(-1/2);
             theta_2=asin(lambda/2/d);
             all_theta=[all_theta theta_2];
             hkl=[hkl;h k l];
         end
     end
 end
for h=-hmax:hmax
     for k=-kmax:kmax
        theta_2=asin(lambda/2/sin(aerfa)*sqrt((h/a)^2+(k/b)^2-2*h*k*cos(aerfa)/a/b));
        all_theta=[all_theta theta_2];
     end
end
[theta_0,loc]=unique(all_theta(find(all_theta<theta_max&all_theta>theta_min)));
theta_num_1=length(theta_0);
theta=[];
for(t=1:theta_num_1)
    theta_1=linspace(theta_0(t)-theta_length,theta_0(t)+theta_length,theta_fine_num);
    theta=[theta theta_1];
end
theta=sort([theta theta_extra]);
theta=unique(theta);
theta_num=length(theta);
front_constant=4*pi*1i/lambda;
flag=zeros(theta_num,1);
root=cell(theta_num,1);

%find possible orientations contributing to the intensity
h = waitbar(0,'Finding the Roots ...');

for q=1:theta_num
    A_max=0;
    A_min=fix(2.*pi.*sin(theta(q)).*a.*cos(pi)./lambda./pi-0.0000001);
    y_A=A_min:A_max;
    root_1=[];
    for r=1:size(y_A,2)
        root1_1=y_A(r)/2/sin(theta(q))/a*lambda;
        root1=root1_1(imag(root1_1)==0);
        root_1=[root_1 root1'];
    end
   if(any((-1-0.000000001<=root_1)&(-1+0.000001>root_1))&&(sin(2*pi*sin(theta(q))*b*(cos(0)*cos(aerfa)))<10^-6))
        flag(q)=1;
        root_1=root_1(root_1>=-1+0.000001);
    end
    num_root=size(root_1,2);
    root_4=[];
    for w=1:num_root
        B_max=fix(2.*sin(theta(q)).*b.*(root_1(w).*cos(aerfa)-sin(aerfa).*sqrt(1-(root_1(w))^2).*cos(pi)+0.000000000001)./lambda);
        B_min=fix(2.*sin(theta(q)).*b.*(root_1(w).*cos(aerfa)-sin(aerfa).*sqrt(1-(root_1(w))^2).*cos(0)-0.000000000001)./lambda);
        y_B=B_min:B_max;
        root_2=[];
        for u=1:size(y_B,2)
            if(((root_1(w).*cos(aerfa)-y_B(u)/2/sin(theta(q))/b*lambda)/sin(aerfa)/sqrt(1-(root_1(w))^2)>1)||((root_1(w).*cos(aerfa)-y_B(u)/2/sin(theta(q))/b*lambda)/sin(aerfa)/sqrt(1-(root_1(w))^2)<-1))
                angle=fix((root_1(w).*cos(aerfa)-y_B(u)/2/sin(theta(q))/b*lambda)/sin(aerfa)/sqrt(1-(root_1(w))^2));
            else
                angle=(root_1(w).*cos(aerfa)-y_B(u)/2/sin(theta(q))/b*lambda)/sin(aerfa)/sqrt(1-(root_1(w))^2);
            end
            root2_1=acos(angle);
            root2_2=[root2_1 (2*pi-root2_1)'];
            root2=root2_2(imag(root2_2)==0);
            root_2=[root_2 root2];
        end
        B_num=size(root_2,2);
        root_3=[repmat(root_1(w),[1 B_num]);root_2];
        root_4=[root_4 root_3];
    end
    root_last=root_4;
    root{q}=root_last;
        waitbar(q/theta_num);
end
close(h);

%read in the coordination and calculate the disperison factors
para_1 =xlsread('Atomic scattering factor for X-ray.xlsx');
para=[zeros(size(para_1,1),1) para_1 ];
locat=unique(Coord(:,1));
coef=XAtomScattering(theta,lambda,locat,para);
times_1=size(locat,1);
num_1=size(Coord,1);
num_2=size(theta,2);
coef_2=zeros(num_1,num_2);
for k=1:times_1
    x=find(ismember(Coord(:,1),locat(k))==1);
    coef_2(x,:)=repmat(coef(k,:),[size(x,1) 1]);
end

%prepare for the integration
x=Coord(:,2);
y=Coord(:,3);
z=Coord(:,4);
k=pi;
f=zeros(theta_num,1);
front_constant=4*pi*1i/lambda;
edge_fun=@(A,B,point_num2_candidate)(fif_fun_1( point_num2_candidate,front_constant,A,a,nx ).*fif_fun_2( point_num2_candidate,front_constant,A,B,b,aerfa,ny).*fif_fun_4( point_num2_candidate,front_constant,A,B,a,b,nz,aerfa,c_2,k_a,k_b ));
% fun=@(A,B,point_num2_candidate)(fif_fun_1( point_num2_candidate,front_constant,A,a,nx ).*fif_fun_2( point_num2_candidate,front_constant,A,B,b,aerfa,ny).*fif_fun_4( point_num2_candidate,front_constant,A,B,a,b,nz,aerfa,c_2,k_a,k_b ));

load('matlab64.mat') % parameters for Gauss Integration Method

%integrating intensity for every theta
h=waitbar(0,'Integrating ...');
for point_num2=1:theta_num
    point_num3=size(root{point_num2},2);
    V=zeros(point_num3,1);
    
    %find divergent points that need to be treated seperately
    if(flag(point_num2)==1);
        edge_length=k*lambda/(nx)/2/pi/sin(theta(point_num2))/a;
        V_1= quad2dg(@(x,y)edge_fun(x,y,theta(point_num2)),-1,-1+edge_length,0,2*pi).*fj_2( theta(point_num2),coef_2,point_num2,front_constant,-1,0,a,b,x,y,z,aerfa,k_a,k_b,c_2);
    else
        V_1=0;
        edge_length=0;
    end
   
    %for other points
    for point_num4=1:point_num3
        %set the integration range for each point
        co_angle_A=root{point_num2}(1,point_num4);
        angle_B=root{point_num2}(2,point_num4);
        x_length=k*lambda/(nx)/2/pi/sin(theta(point_num2))/a*20;
        x_edge_1=max(-1+edge_length,co_angle_A-x_length);
        x_edge_2=min(0,co_angle_A+x_length);
        if angle_B>=0&&angle_B<=pi;
        y_edge_1=acos(min(1,cos(angle_B)+10*lambda/2/sin(theta(point_num2))/b/(ny)/sin(aerfa)/sqrt(1-co_angle_A^2)));
        y_edge_2=acos(max(-1,cos(angle_B)-10*lambda/2/sin(theta(point_num2))/b/(ny)/sin(aerfa)/sqrt(1-co_angle_A^2)));
        else
        y_edge_1=2*pi-acos(max(-1,cos(angle_B)-10*lambda/2/sin(theta(point_num2))/b/(ny)/sin(aerfa)/sqrt(1-co_angle_A^2)));
        y_edge_2=2*pi-acos(min(1,cos(angle_B)+10*lambda/2/sin(theta(point_num2))/b/(ny)/sin(aerfa)/sqrt(1-co_angle_A^2)));
        end
        %calculate the prefactor
        fj_psra=fj_2( theta(point_num2),coef_2,point_num2,front_constant,co_angle_A,angle_B,a,b,x,y,z,aerfa,k_a,k_b,c_2);
        %perform the integration, where Co_1 and Co_2 are weights for Gauss
        %Integration, fif_fun_1,2,4 are the intensity functions for a b c directions 
        func=Co_1.*Co_2.*fif_fun_1( theta(point_num2),front_constant,(x_edge_1-x_edge_2)./2.*X+(x_edge_1+x_edge_2)./2,a,nx ).*fif_fun_2( theta(point_num2),front_constant,(x_edge_1-x_edge_2)./2.*X+(x_edge_1+x_edge_2)./2,(y_edge_1-y_edge_2)./2.*Y+(y_edge_1+y_edge_2)./2,b,aerfa,ny).*fif_fun_4( theta(point_num2),front_constant,(x_edge_1-x_edge_2)./2.*X+(x_edge_1+x_edge_2)./2,(y_edge_1-y_edge_2)./2.*Y+(y_edge_1+y_edge_2)./2,a,b,nz,aerfa,c_2,k_a,k_b );
        V(point_num4)=sum(sum(func.*fj_psra)).*(x_edge_1-x_edge_2).*(y_edge_1-y_edge_2).*0.25;
    end
     f(point_num2)=2.*sum(V)+2*V_1;
     waitbar(point_num2/theta_num,h,sprintf('Integrating... %.2f %%',100*point_num2/theta_num));

end
close(h);
%plot the simulation result
t=2*theta.*360./(2.*pi);
if length(theta) == length(unique(theta)) && length(t) == length(unique(t))
    [fitresult, gof] = createFit(t',f);
    t_1=linspace(min(t)+0.001,max(t)-0.001,10000);
    f_1=feval(fitresult,t_1);
    [t_e,loc_e]=sort([t t_1]);
    f_e=[f;f_1];
    f_e=f_e(loc_e);
    plot(t_e,f_e);
    t=t';
    title(strcat('Simulated PXRD pattern for   ',samplename));
    xlabel('2-theta');
    ylabel('Intensity');
    twotheta=t_e;
    intensity=f_e;
    save(samplename,'twotheta','intensity');
else
    plot(t,f);
    t=t';
    title(strcat('Simulated PXRD pattern for   ',samplename));
    xlabel('2-theta');
    ylabel('Intensity');
    twotheta=t;
    intensity=f;
    save(samplename,'twotheta','intensity');
end

% matlabpool close;
