wm = input('Enter a w_m value:');

wn = input('Enter a w_n value:');

if wm == 0 && wn == 0 
    disp('Case 1')
elseif wm ~=0 && wn==0
    disp('Case 2')
elseif wm==0 && wn~=0
    disp('Case 3')
elseif wm~=0 && wn~=0 && wm==wn
    disp('Case 4')
elseif wm~=0 && wn~=0 && wm~=wn
    disp('Case 5')
end
