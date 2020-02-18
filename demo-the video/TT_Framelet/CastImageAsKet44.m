function [X] = CastImageAsKet44(M,Nway,I1,J1)
    i=1;j=1; % Nway = [16,16,16,16,3]
    N = 1:numel(Nway); % N = [1 2 3 4 5]
    X = zeros(Nway);

                    for i4 = 1: Nway(4)
                        [i,j] = index_ij(i,j,i4,N(4));
                        for i3 = 1: Nway(3)
                            [i,j] = index_ij(i,j,i3,N(3));
                            for i2 = 1: Nway(2)
                                [i,j] = index_ij(i,j,i2,N(2));
                                Mtemp = M(i:i+I1-1,j:j+J1-1,:);
                                X(i4,i3,i2,:) = Mtemp(:);
                            end
                        end
                    end
end

function [iN,jN] = index_ij(i0,j0,k_N,N)
    switch k_N
        case 1
            iN = i0;
            jN = j0;
        case 2
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 3
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 4
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;            
        case 5
            [K] = sumi3(N);
            iN = i0-K;
            jN = j0+4;
        case 6
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 7
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 8
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 9
            [K] = sumi3(N);
            iN = i0-K;
            jN = j0+4;
        case 10
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 11
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 12
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 13
            [K] = sumi3(N);
            iN = i0-K;
            jN = j0+4;
        case 14
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 15
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
        case 16
            iN = i0+4;
            K = sumi2(N);
            jN = j0+1 - K;
    end
end

function [K] = sumi2(N)
    K=4^(N-1)-3;
%     for n = 2:N
%         K = K+4^(n-2);
%     end
end

function [K] = sumi3(N)
    K=0;
    for n = 2:N
        K = K+3*4^(n-1);
    end
end

% function [K] = sumi4(N)
%     K=0;
%     for n = 2:N
%         K = K+3*4^(n-1);
%     end
% end
% 
% function [K] = sumi5(N)
%     K=0;
%     for n = 2:N
%         K = K+4*4^(n-1);
%     end
% end