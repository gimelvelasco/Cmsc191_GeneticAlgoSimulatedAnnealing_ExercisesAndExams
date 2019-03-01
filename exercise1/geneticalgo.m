function geneticalgo
%%
%Visualize a 10x8x40 box
%10 Generations with
%40 Chromosomes per generation consisting of
%8 bits (1 sign bit and 7 value bits) per chromosome
%
%The 3-Dimensional Array Variable named "TS" w/c stands for "Time-Space"
%will hold the entire chromosomes in each generation
%%
%The Fittest Chromosome in a given generation is the one who has the value
%closest to 0.
%The Fitness Function is defines as:
%           f(x) = x^2
%on the interval [-100,100]
%%
%Initializing TS box to zero
%i - generation counter
%j - chromosome counter
%k - bit counter
%Format: TS(i,k,j)
for i=1:10  %across 10 generations
    for j=1:40  %with a chromosome population of 40 in each generation
       %TS(i,:,j) = zeros(1,8);   %with 8 bits per chromosome
    end
end
%%
%INITIALIZATION
for j=1:40
    TS(1,:,j) = round(rand(1,8));   %first generation
    while DECODE(TS(1,:,j)) > 100 || DECODE(TS(1,:,j)) < -100%this will prevent values exceeding +-100
        TS(1,:,j) = round(rand(1,8));
    end
end
%%
for i=1:10  %for a maximum of 10 generations
%SELECTION
TS = SELECTION(TS,i);
%User Display
PRINT(TS,i);
%CROSSOVER
TStmp = CROSSOVER(TS,i);
TStmp = SELECTION(TStmp,i);%Ranking
%INSERTION
TS = INSERTION(TS,TStmp,i);
TS = SELECTION(TS,i);%Ranking
%Moves on to the next generation
TS(i+1,:,:) = TS(i,:,:);
%STOPPING CRITERIA %assuming we dont know the answer.
if i > 1    %prevents segmentation fault
    %Mean Deviation Stopping Criteria is implemented
    if GETAVEG(TS,i) < 30 %compares the mean fitness of current and previous generations
        PRINT(TS,i+1);
        fprintf('Termination by Population Mean Deviation.\n');
        break; %if equal, break.
    end
end
end
if i == 10
    fprintf('Termination by Generations Completed.\n'); 
end
%%
%Displays the Solution to the problem.
fprintf('1st Generation Fittest: %d\tGenerations Completed: %d',GETFITNESS(TS(1,:,1)),i);
fprintf('\nThe Most Fit Solution is %d\n',GETFITNESS(TS(i,:,1)));
end

function TSret = SELECTION(TSrec,i)
    %Uses Selection by Fitness
    %bubblesort
    %for sorting, an ascending order is used w/c means that the fittest
    %would be found in the first index of the array
    for m=1:39
       for n=1:39
           if GETFITNESS(TSrec(i,:,n)) > GETFITNESS(TSrec(i,:,n+1))
               %swap
               temp = TSrec(i,:,n);
               TSrec(i,:,n) = TSrec(i,:,n+1);
               TSrec(i,:,n+1) = temp;
           end
       end
    end
    TSret = TSrec;
end

function TSret = CROSSOVER(TSrec,i)
    %Crossover first 15 chromosomes to the next 15 chromosomes
    %Single Point Crossover. 0000|0000
    for m=1:15
        for n=1:4
           temp(n,m) = TSrec(i,n,m);
           TSrec(i,n,m) = TSrec(i,n+4,m+15);
           TSrec(i,n+4,m+15) = temp(n,m);
        end
    end
    TSret = TSrec;
end

function TSret = INSERTION(TSrec,TStmprec,i)
    %replace the last 15 with the 15 new offsprings
    %Elitist Strategy used as the technique for insertion
    for j=1:15
       TSrec(i,:,j+25) = TStmprec(i,:,j);
    end
    TSret = TSrec;
end

function PRINT(TSrec,i)
    fprintf('Generation %d Fitness Ranking: (w/ Aveg Fitness of %.2f and w/ %d as the Fittest)\n',i,GETAVEG(TSrec,i),GETFITNESS(TSrec(i,:,1)));
    for j=1:40
       fprintf('%d ',GETFITNESS(TSrec(i,:,j))); 
    end
    fprintf('\n\n');
end

function AVEGFITret = GETAVEG(TSrec,i)
    AVEGFITret = 0;
    for j=1:40
        AVEGFITret = AVEGFITret + GETFITNESS(TSrec(i,:,j));
    end
    AVEGFITret = AVEGFITret/40;
end

function VALUEret = DECODE(chrmsm)
    %%
    %DECODING. Converting binary to integer
    VALUEret = 0;
    for k=2:8   %computes the value of the chromosome
       VALUEret = VALUEret + (2^(8-k))*chrmsm(k);
    end
    if chrmsm(1) == 1   %if sign bit is 1, value is negative. otherwise positive.
        VALUEret = (-1)*VALUEret;
    end
end

function FITNESSret = GETFITNESS(chrmsm)    %FITNESS FUNCTION.
    VALUEret = DECODE(chrmsm);
    %%
    %Actual Fitness Function f(x) = x^2
    FITNESSret = VALUEret^2;
end