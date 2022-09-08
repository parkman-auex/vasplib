classdef CmdLineProgressBar < handle
% class for command-line progress-bar notification.
% Use example:
%   pb = CmdLineProgressBar('Doing stuff...');
%   for k = 1 : 10
%       pb.print(k,10)
%       % do stuff
%   end
%
% Author: Itamar Katz, itakatz@gmail.com
% Modified: parkman

    properties
        last_msg_len = 0;
    end
    methods
        %--- ctor
        function obj = CmdLineProgressBar(msg)
            fprintf('%s', msg)
        end
        %--- print method
        function print(obj, n, tot,msg_tail)
            if nargin < 4
                msg_tail = '';
            end
            if numel(n)>1 && numel(tot) == numel(n)
                %mulity
                fprintf('%s', char(8*ones(1, obj.last_msg_len))) % delete last info_str
                info_str = [''];
                for i = 1: numel(n)
                    info_str = [info_str,num2str(n(i)),char(tot{i})];
                end
                info_str = [info_str,msg_tail];
                fprintf('%s', info_str);
                %--- assume user counts monotonically
                if n(1) == n(2)
                    fprintf('\n')
                end
                obj.last_msg_len = length(info_str);
            else
                fprintf('%s', char(8*ones(1, obj.last_msg_len))) % delete last info_str
                if length(n) == 1 && length(tot) == 1
                    info_str = [sprintf('%d/%d', n, tot),msg_tail];
                else
                    info_str = [num2str(n),'/',num2str(tot),msg_tail];
                end
                fprintf('%s', info_str);
                %--- assume user counts monotonically
                if n == tot
                    fprintf('\n')
                end
                obj.last_msg_len = length(info_str);
            end
        end
        %--- dtor
        function delete(obj)
            %fprintf('%s', char(8*ones(1, obj.last_msg_len))) % delete last info_str
            delete(obj);
            %fprintf('\n');
        end
    end
end