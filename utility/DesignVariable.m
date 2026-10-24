classdef DesignVariable

    properties (SetAccess = private)
        value   double
        indx    double {mustBePositive, mustBeInteger}
        rootLen (1,1) double {mustBeNonnegative, mustBeInteger}
        Jac     (:,:) double
    end

    methods

        function obj = DesignVariable(value, indx, rootLen, Jac)
            if nargin < 2
                indx = reshape(1:numel(value), size(value));
            elseif ~isequal(size(value), size(indx)) && ~isempty(indx)
                try
                    indx = reshape(indx, size(value));
                catch exception
                    error('value and indx must have the same number of elements.')
                end
            end

            if nargin < 3 || isempty(rootLen)
                rootLen = numel(value);
            end

            if any(indx > rootLen)
                error('indx cannot exceed rootLen.')
            end

            if nargin < 4 || isempty(Jac)
                Jac = sparse(1:numel(value),indx(:),1,numel(value),rootLen);
            elseif ~isequal(size(Jac), [numel(value),rootLen])
                error('Jac must have size [numel(value), rootLen].')
            end

            obj.value = value;
            obj.indx = indx;
            obj.rootLen = rootLen;
            obj.Jac = Jac;
        end
        function out = subsref(obj, S)
            switch S(1).type
                case '()'
                    % S(1).subs contains the indices
                    out = obj.at(S(1).subs{:});
    
                case '.'
                    out = builtin('subsref', obj, S);
    
                otherwise
                    error("Unsupported indexing type.")
            end

            if numel(S) > 1 && strcmp(S(1).type, '()')
                out = subsref(out, S(2:end));
            end
        end

        % Overloaded operators
        function out = at(obj, indx)
            out = DesignVariable(obj.value(indx), obj.indx(indx), obj.rootLen, obj.Jac(indx,:));
        end

        function obj = reshape(obj, varargin)
            obj.value = reshape(obj.value, varargin{:});
            obj.indx  = reshape(obj.indx, varargin{:});
        end

        function out = exp(obj)
            expx = exp(obj.value);
            out = DesignVariable(expx, [], obj.rootLen, expx*obj.Jac);
        end

        function out = size(obj)
            out = size(obj.value);
        end

        % NOTE: cannot override numel() as this causes indexing errors.
        %       use prod(size(obj)) or numel(obj.value) or size(obj.Jac, 1)

        function obj = vertcat(varargin)
            N = length(varargin);
            values = cell(N,1);
            indxes = cell(N,1);
            Jacs   = cell(N,1);
            for i = 1:N
                if i == 1
                    rtLen = varargin{i}.rootLen;
                elseif varargin{i}.rootLen ~= rtLen
                    error('All DesignVariables must have the same rootLen.')
                end
                values{i} = varargin{i}.value;
                indxes{i} = varargin{i}.indx;
                Jacs{i}   = varargin{i}.Jac;
            end
            obj = DesignVariable(vertcat(values{:}), vertcat(indxes{:}), rtLen, vertcat(Jacs{:}));
        end
    end
end