function [h, r, epsilon] = make_input(points, connection_probability)
    % points: nx2 행렬
    % connection_probability: 0과 1 사이의 값
    
    % 1. points를 flat하게 만들어 h 생성
    h = reshape(points.',[], 1)

    % 2. 점들 사이의 거리를 계산하여 r 생성
    n = size(points, 1);
    edge_idx = 1;

    for i = 1:n
        for j = 1:n
            if i == j
                continue
            end

            if rand <= connection_probability
                % 점 i와 j 사이의 거리를 계산하여 r에 저장
                distance = sqrt(sum((points(i, :) - points(j, :) ).^2));
                r(edge_idx, 1) = distance;
                
                epsilon(edge_idx, 1:2) = [i, j];
                edge_idx = edge_idx + 1;
            end
        end
    end

    h = h + normrnd(0, 5, size(h));
    r = r + normrnd(0, 0.001, size(r));
end