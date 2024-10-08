%% 更新粒子群的位置和速度
function [new_particle,velocity] = update_particle_position(particle, pbest, gbest, w, c1, c2, velocity,vmax,max_num_cluster)
%     num_buildings = numel(particle); % 计算传入的粒子总数
    % 更新速度
    velocity = w * velocity + c1.* (pbest - particle) + c2.* (gbest - particle);
    
    % 更新位置
    new_particle = round(particle + velocity);
    
    % 边界处理
    velocity=min(velocity,vmax);
    velocity=max(velocity,-vmax);
    new_particle(new_particle < 1) = 1;
    new_particle(new_particle > max_num_cluster) = max_num_cluster;
end

