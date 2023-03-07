function A = read_transportation_network(dataset_name)

% addpath('Datasets/TransportationNetworks/')
% dataset_name = 'Barcelona';
% filename_node = sprintf('Datasets/TransportationNetworks/%s/%s_node.tntp',dataset_name,dataset_name);
filename_net = sprintf('Datasets/TransportationNetworks/%s/%s_net.tntp',dataset_name,dataset_name);


% a=dlmread(filename_node,'\t',1,0);
% z=a(:,1);
% x=a(:,2);
% y=a(:,3);

a=textread(filename_net, '%[^\n]\n');

strn = regexp(a{2},'\d*','match');
n = str2double(strn{1});

j=0;
for i=1:length(a)
    if(length(a{i}>0))
        c=a{i}(1);
        if(c~='~' & c~='<' & c~=' ' & c~='\t' & c~='\n')            
            b=sscanf(a{i},'%f');
            j=j+1;
            Tail(j,1)=b(1);
            Head(j,1)=b(2);
        end
    end
end

%Edges = [Tail, Head];
A = sparse(Tail,Head,1,n,n);

end

% save 'ChicagoRegionalBaseNet.mat' AllLinks x y


% fid  = fopen(filename, 'r') ;
% nodeinfo = cell(1e6, 4) ;                    % Prealloc.
% rCnt = 0 ;                               % Row counter.
% 
% while ~feof(fid)
%     rCnt = rCnt + 1 ;
%     nodeinfo{rCnt,1} = fscanf(fid, '%d', 1) ;
%     nodeinfo{rCnt,2} = fscanf(fid, '%s', 1) ;
%     nodeinfo{rCnt,3} = fscanf(fid, '%s', 1) ;
%     nodeinfo{rCnt,4} = fscanf(fid, '%s', 1) ;
% end
% 
%  fclose(fid) ;
%  nodeinfo = nodeinfo(1:rCnt,:) ;                  % Truncate.
% 
