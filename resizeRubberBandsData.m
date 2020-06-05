function [newData] = resizeRubberBandsData(PmatData)
newData = PmatData;
%Defining number of materials
N = length( fieldnames(PmatData) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(PmatData);
%Useful constants
test1 = 1;
test2 = 1;
test3 = 1;
%Thickness is in millimeters
Nat100RThickness = [0, 0, 0, 0, 0, 0, ...
    0.33, 0.42,...  %Yellow
    0.48, 0.52,...  %Red
    0.64, 0.75,...  %Blue
    0.97, 1.03,...  %Green
    1.17, 1.31,...  %Black
    1.32, 1.49];    %Orange

for i=7:N    
    M = length( fieldnames(PmatData. ( fields.N{i} ) ) );
    fields.M = fieldnames(PmatData.( fields.N{i} ) );
    newData = rmfield(newData,(fields.N{i}));
    for j=1:M
        P = length( PmatData. ( fields.N{i} ).( fields.M{j}) );
        for k=1:P
            switch (fields.M{j}) %disR50 disR250 disR500
                case {"disR50","L5min"}
                    testCount = test1;
                    test1 = test1 + 1;
                case {"disR250","L15min"}
                    testCount = test2;
                    test2 = test2 + 1;
                case {"disR500","L180min"}
                    testCount = test3;
                    test3 = test3 + 1;
            end            
            newData.Nat100R.(fields.M{j})(testCount).paths = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).paths;
            newData.Nat100R.(fields.M{j})(testCount).rload = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).rload;
            newData.Nat100R.(fields.M{j})(testCount).rdis = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).rdis;
            newData.Nat100R.(fields.M{j})(testCount).rtime = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).rtime;
            newData.Nat100R.(fields.M{j})(testCount).pload = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).pload;
            newData.Nat100R.(fields.M{j})(testCount).pdis = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).pdis;
            newData.Nat100R.(fields.M{j})(testCount).ptime = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).ptime;
            newData.Nat100R.(fields.M{j})(testCount).sload = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).sload;
            newData.Nat100R.(fields.M{j})(testCount).stats = ...
                PmatData.(fields.N{i}).(fields.M{j})(k).stats;
            newData.Nat100R.(fields.M{j})(testCount).stats.thickness = Nat100RThickness(i);
            newData.Nat100R.(fields.M{j})(testCount).name = ...
                fields.N{i};
        end
    end
end
%Finding Smallest Test
N = length( fieldnames(newData) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(newData);
for i=7
    M = length( fieldnames(newData. ( fields.N{i} ) ) );
    fields.M = fieldnames(newData.( fields.N{i} ) );
    for j=1:M
        P = length( newData. ( fields.N{i} ).( fields.M{j}) );
        for k=1:P
            u(k) = size(newData.(fields.N{i}).(fields.M{j})(k).pload,2);
        end
        smallest = min(u);
        for k=1:P
            newData.(fields.N{i}).(fields.M{j})(k).stats.smallest = smallest;
        end
    end
end
end