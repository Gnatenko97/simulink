data=readtable('25PDB.csv');
new_table=data(3:height(data),1:4);
names=table2array(new_table(:,1));
sequence=table2array(new_table(:,2));
coded=table2array(new_table(:,3));
result=table2array(new_table(:,4));

ind=500; string=''; cod='';

%string=cell2mat(sequence(ind,1));
%cod=cell2mat(coded(ind,1));
for ind=1:1:10
string=[string,cell2mat(sequence(ind,1))];
cod=[cod,cell2mat(coded(ind,1))];
end

char_array=char(string);
[Helic_S, Coilic_S, Sheets_S]=SingleCorelation(char_array,cod);

arr=''; coil=1;
Helic=''; Coilic=''; Sheets='';
for i=2:1:length(char_array)-1
if cod(i)=='H' && cod(i+1)=='H'
    Helic=[Helic,char_array(1,i)];
elseif cod(i)=='C' && cod(i+1)=='C'
    Coilic=[Coilic,char_array(i)];
elseif cod(i)=='E' && cod(i+1)=='E'
    Sheets=[Sheets,char_array(i)];
end 
end

[List_of_Helic]=Double(Helic,1);
[List_of_Sheets]=Double(Sheets,1);

Ploting(List_of_Helic,'Helical');
Ploting(List_of_Sheets,'Sheets');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total=length(Helic_S)+length(Sheets_S)+length(Coilic_S);
p_h=length(Helic_S)/total*100;
p_e=length(Sheets_S)/total*100;
p_c=length(Coilic_S)/total*100;

v1=['Amino Acids forming Alpha-Helical structures: ',num2str(length(Helic_S)),' (',num2str(p_h),'%)'];
v2=['Amino Acids forming Beta-Sheets structures: ',num2str(length(Sheets_S)),' (',num2str(p_e),'%)'];
v3=['Amino Acids forming T-turns S-bend, or unassigned structures: ',num2str(length(Coilic_S)),' (',num2str(p_c),'%)'];
nm1=['Protein id: ',cell2mat(names(ind))];
nm2=['Number of AA: ',num2str(length(cod))];

disp(nm1);
disp(nm2);
disp(v1);
disp(v2);
disp(v3);
disp(char_array);
disp(cod);


function[]=Ploting(List_of_Helic,a)
Y=List_of_Helic(:,2);
Yy=cell2mat(Y);
X=List_of_Helic(:,1);
X=X';
order_table=cell2table(X);
order_matrix=table2array(order_table);
Xx = categorical(X,order_matrix);
figure;
t1=['Bargraph of ',a,' Amino Acid Freqency'];
bar(Xx,Yy); title(t1);
xlabel('Amino Acids'); ylabel('Freqency/occurance');
end
function [out]=Double(Helic,spac)
cell={};
for i=1:1:length(Helic)-spac
h2=Helic(i:i+spac);
[a,b]=size(cell);
plI=a+1;
for j=1:1:a
if char(cell(j,1))== h2
    plI=j;
end
end

if plI==a+1
cell(plI,1)={h2}; cell(plI,2)={1};
else
    cell(plI,2)={cell2mat(cell(plI,2))+1};
end
end

out = sortrows(cell,2,'descend');
end
function[Helic,Coilic,Sheets]=SingleCorelation(char_array,cod)
Helic=''; Coilic=''; Sheets='';
for i=1:1:length(char_array)
if cod(i)=='H'
    Helic=[Helic,char_array(1,i)];
elseif cod(i)=='C'
    Coilic=[Coilic,char_array(i)];
elseif cod(i)=='E'
    Sheets=[Sheets,char_array(i)];
end 
end
%C-ether a T turn, S-bend, or unassigned  
%H-alpha helix
%E-beta Bridge
end