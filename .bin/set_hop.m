%%  usage:  H_xyz = set_hop(amp,hi,hj,vector_list,H_xyz,mode)
function H_xyz = set_hop(amp,hi,hj,vector_list,H_xyz,mode)
    if nargin <6
        mode = 'set';
    end
    [NRPTS,~]=size(H_xyz);
    WAN_NUM=length(H_xyz(1).Hnum);
    WAN_NUM_half = WAN_NUM/2;
    [n_vector,~] = size(vector_list);
    %
    V = [H_xyz.vector];
    V =reshape(V,3,length(V)/3)';
    % disp(vector);
    % disp(V);
    for i = 1:n_vector
        vector = vector_list(i,:);
        [~,seq]=ismember(vector,V,'rows');

        %% test amp
        length_amp = length(amp);
        if strcmp(mode,'mat')
            if seq == 0
                seq = NRPTS +1;
                H_xyz(seq,:).seq = seq;
                H_xyz(seq,:).vector = vector;
                H_xyz(seq,:).Degen = 1;
                H_xyz(seq,:).Hnum = zeros(WAN_NUM);
                %H_xyz(seq,:).Hsym = sym(zeros(WAN_NUM));
                H_xyz(seq,:).Hcoe = sym(zeros(WAN_NUM));
            end
            H_xyz(seq,:).seq = seq;
            H_xyz(seq,:).vector = vector;
            H_xyz(seq,:).Degen = 1;
            H_xyz(seq,:).Hnum = amp  ;
        elseif strcmp(mode,'symmat')  
            if seq == 0
                seq = NRPTS +1;
                H_xyz(seq,:).seq = seq;
                H_xyz(seq,:).vector = vector;
                H_xyz(seq,:).Degen = 1;
                H_xyz(seq,:).Hnum = zeros(WAN_NUM);
                %H_xyz(seq,:).Hsym = sym(zeros(WAN_NUM));
                H_xyz(seq,:).Hcoe = sym(zeros(WAN_NUM));
            end
            H_xyz(seq,:).seq = seq;
            H_xyz(seq,:).vector = vector;
            H_xyz(seq,:).Degen = 1;
            H_xyz(seq,:).Hcoe = amp  ;
       
        elseif length_amp == 1
            % single mode
            H_xyz = set_hop_single(amp,hi,hj,vector,H_xyz,mode);
        elseif length_amp == 2
            % spin mode
            H_xyz = set_hop_single(amp(1,1),hi,hj,vector,H_xyz,mode);
            H_xyz = set_hop_single(amp(1,2),hi,hj+WAN_NUM_half,vector,H_xyz,mode);
            H_xyz = set_hop_single(amp(2,1),hi+WAN_NUM_half,hj,vector,H_xyz,mode);
            H_xyz = set_hop_single(amp(2,2),hi+WAN_NUM_half,hj+WAN_NUM_half,vector,H_xyz,mode);
        elseif length_amp == 4
            % four-band mode
            disp('if needed');
        elseif length_amp == 8
            % eight-band mode
            disp('if needed');
        else
        end
    end

end

function H_xyz = set_hop_single(amp,hi,hj,vector,H_xyz,mode)
    [NRPTS,~]=size(H_xyz);
    WAN_NUM=length(H_xyz(1).Hnum);
    WAN_NUM_half = WAN_NUM/2;
    %
    V = [H_xyz.vector];
    V =reshape(V,3,length(V)/3)';

    [~,seq]=ismember(vector,V,'rows');
    if strcmp(mode,'mat')
        if seq == 0
            seq = NRPTS +1;
            H_xyz(seq,:).seq = seq;
            H_xyz(seq,:).vector = vector;
            H_xyz(seq,:).Degen = 1;
            H_xyz(seq,:).Hnum = zeros(WAN_NUM);
            %H_xyz(seq,:).Hsym = sym(zeros(WAN_NUM));
            H_xyz(seq,:).Hcoe = sym(zeros(WAN_NUM));
        end
        H_xyz(seq,:).seq = seq;
        H_xyz(seq,:).vector = vector;
        H_xyz(seq,:).Degen = 1;
        H_xyz(seq,:).Hnum = amp  ;
        return
    else
        % for new block
        if seq == 0
            seq = NRPTS +1;
            H_xyz(seq,:).seq = seq;
            H_xyz(seq,:).vector = vector;
            H_xyz(seq,:).Degen = 1;
            H_xyz(seq,:).Hnum = zeros(WAN_NUM);
            %H_xyz(seq,:).Hsym = sym(zeros(WAN_NUM));
            H_xyz(seq,:).Hcoe = sym(zeros(WAN_NUM));
        end

        if strcmp(mode,'set')
            if H_xyz(seq).Hnum(hi,hj) ~= 0
                warning('You should use add mode on this NRPT and hi hj');
                fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3),hi,hj);
            end
            H_xyz(seq).Hnum(hi,hj) = amp ;
        elseif strcmp(mode,'add')
            H_xyz(seq).Hnum(hi,hj) = H_xyz(seq).Hnum(hi,hj) + amp  ;
        elseif strcmp(mode,'sym')
            H_xyz(seq).Hcoe(hi,hj) =  amp  ;
        elseif strcmp(mode,'symadd')
            H_xyz(seq).Hcoe(hi,hj) = H_xyz(seq).Hcoe(hi,hj) + amp  ;
        end
    end
end
