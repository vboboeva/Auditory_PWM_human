% [varargout] = extract_from_struct(s, varargin)
%
% Does a deal.m operation on structs, extracting numeric vectors out of the
% struct. For example, if a structure A contains a fielname 'this', such
% that A(1).this = 10; A(2).this = 30; and so on, then
%
%     this = extract_from_struct(A, 'this')
%
% will return a vector with [10 ; 30 ; etc.] as its elements. The length of
% the returned vector will be the same as the length of A.
%
% You can extract mutliple things at once. For example, if A also contains
% a fieldname 'that', then you can call
%
%     [this, that] - extract_from_struct(A, 'this', 'that')
%

function [varargout] = extract_from_struct(s, varargin)

for i=1:numel(varargin),
   x = cell(numel(s),1); 
   [x{:}] = deal(s.(varargin{i})); 
   varargout{i} = cell2mat(x);
end;
