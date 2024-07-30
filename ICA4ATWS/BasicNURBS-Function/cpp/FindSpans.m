function [spans,isend] = FindSpans(knots,p,U,offset)
%% Instruction of programs
%
% FINDSPANS:
%	 Find the spans of B-spline knot vector at the parameter points.
%
% Calling Sequence:
%    [spans,isend] = FindSpans(knots,p,U,offset)
%
% INPUT:
%    knots  : Non-decreasing knot sequence
%    p      : Degree of B-spline
%    U      : Parameter point
%    offset : Offset value
%
% OUTPUT:
%    spans  : Span index corresponding to the parameter point
%    isend  : Whether the parameter point is at the end of the interval
%
%% Body of programs
if min(U) < knots(p+1) || max(U) > knots(end-p)
    error('Some value is outside the knot-vector span.')
end

nu = length(U);
spans = zeros(nu,1);
isend = zeros(nu,1);

% Compute span index
if nargin == 3
    for i = 1:nu
        if U(i) == knots(end-p)
            spans(i) = find(knots < U(i),1,'last');
            isend(i) = 1;
        else
            spans(i) = find(knots <= U(i),1,'last');
            if knots(spans(i)) == U(i)
                isend(i) = 1;
            end
        end
    end
else % Add offset
    for i = 1:nu
        if U(i) == knots(end-p)
            spans(i) = find(knots < U(i),1,'last')+offset(i);
            isend(i) = 1;
        else
            spans(i) = find(knots <= U(i),1,'last')-offset(i);
            if knots(spans(i)) == U(i)
                isend(i) = 1;
            end
        end
    end
end

end
%% Demo of programs
% knots = [0 0 0 0.2 0.4 0.6 0.8 1 1 1];
% p = 2;
% U = [0.1 0.3 0.5 0.7 0.9];
% spans = FindSpans(knots,p,U);
%