function Dice = DiceFunction(Left_Atrium , participants)

Dice = zeros(1 , participants); % Final dice vector

for i = 1:participants

    k = [1 2 3];
    k(i) = [];
    % k = [2 3] --> [1 3] --> [1 2]
    % [Ben Eric] --> [Jude Eric] --> [Jude Ben]

    overlap = sum( Left_Atrium{k(1)} .* Left_Atrium{k(2)} , "all" ); 
    % product finds overlap (1*1=1 , 1*0=0 , 0*0=0)

    total = sum( Left_Atrium{k(1)} + Left_Atrium{k(2)} , "all" );
    % sum finds total

    Dice(i) = (2 * overlap) / total; % dice equation 
    % 1 = perfect overlaap , 0 = no overlap

end