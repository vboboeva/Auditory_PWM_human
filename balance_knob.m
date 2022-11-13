% [snd] = balance_knob(balance, snd)
%
% Takes a stereo sound in the 2-by-n matrix snd, and passes it through a
% balance knob: -1 means all left, 0 means equal, and +1 means all right.
% The multiplication is done so that the sum of the square of the factors
% multiplying left and right channels is kept constant, independent of
% balance.
%

function [snd] = balance_knob(balance, snd)

if balance<-1, balance = -1; end;
if balance>1,  balance = 1;  end;

if balance < -1 + 1e-10,
   L2 = 2; R2 = 0;
else
   alpha = (1 - balance) / (1 + balance);
   L2 = 2*alpha/(1+alpha);
   R2 = 2/(1+alpha);
end;

snd(1,:) = sqrt(L2)*snd(1,:);
snd(2,:) = sqrt(R2)*snd(2,:);
