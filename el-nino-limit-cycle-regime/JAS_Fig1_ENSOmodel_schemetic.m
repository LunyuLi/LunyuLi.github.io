function enso_diagram_only()
    figure('Color', 'w', 'Position', [100, 100, 1100, 600]);
    hold on; axis equal; axis off;
    font_name = 'Times New Roman';
    font_title = 28;
    xlim([-9, 12]);
    ylim([-6, 8]);
    draw_real_enso_diagram();
    draw_vertical_temp_profile();
    set(findall(gcf, '-property', 'FontName'), 'FontName', font_name);
    drawnow;

    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end
    outBase = fullfile(scriptDir, 'JAS_Fig1_ENSOmodel_schemetic');
    exportgraphics(gcf, [outBase, '.png'], 'Resolution', 600);
    print(gcf, [outBase, '.pdf'], '-dpdf', '-r600', '-bestfit');
    savefig(gcf, [outBase, '.fig']);
    disp(['Saved figure: ', outBase, '.png']);

    figure('Color', 'w', 'Position', [100, 100, 1100, 600]);
    hold on; axis equal; axis off;
    xlim([-9, 12]);
    ylim([-6, 8]);
    draw_real_enso_diagram();
    draw_vertical_temp_profile('formula');
    set(findall(gcf, '-property', 'FontName'), 'FontName', font_name);
    drawnow;

    outBaseFormula = fullfile(scriptDir, 'JAS_Fig1_ENSOmodel_schemetic_chapter2_formula');
    exportgraphics(gcf, [outBaseFormula, '.png'], 'Resolution', 600);
    print(gcf, [outBaseFormula, '.pdf'], '-dpdf', '-r600', '-bestfit');
    savefig(gcf, [outBaseFormula, '.fig']);
    disp(['Saved figure: ', outBaseFormula, '.png']);
end
function draw_real_enso_diagram()
    col_sun     = [0.8, 0.1, 0.1];   
    col_atmos   = [0.6, 0.6, 0.6];   
    col_t1      = [1.0, 0.4, 0.4];   
    col_t2      = [1.0, 0.6, 0.6];  
    col_sub     = [0.75, 0.85, 0.95];
    col_deep    = [0.0, 0.45, 0.75]; 
    col_arrow   = 'k'; 
    lw_border = 1.5;
    sun_x = 0; sun_y = 6.0; sun_r = 0.7;
    for ang = 0:45:315
        rad = ang * pi/180;
        r_in = sun_r + 0.2;
        r_out = sun_r + 0.6;
        plot([sun_x + r_in*cos(rad), sun_x + r_out*cos(rad)], ...
             [sun_y + r_in*sin(rad), sun_y + r_out*sin(rad)], ...
             'k-', 'LineWidth', 2);
    end
    viscircles([sun_x, sun_y], sun_r, 'Color', 'r', 'LineWidth', 1);
    fill(sun_x + sun_r*cos(linspace(0,2*pi,50)), ...
         sun_y + sun_r*sin(linspace(0,2*pi,50)), col_sun, 'EdgeColor', 'r');
    rect_x = -5.5; rect_y = 2.5; rect_w = 11; rect_h = 2.0;
    rectangle('Position', [rect_x, rect_y, rect_w, rect_h], 'Curvature', 0.2, ...
              'FaceColor', col_atmos, 'EdgeColor', 'none');
    text(0, rect_y + rect_h/2, 'Atmosphere', 'Horiz', 'center', 'FontSize', 28);
    x = linspace(-5.5, 5.5, 200);
    y_surf = 0 * cos(x * 0.57-4.8) + 1.2;
    y_therm = 0 * cos(x * 0.57-4.8) + 0.167;
    y_deep_bound = -2.42 + 1.55 * tanh(0.4 * (x - 0.85)); 
    y_bot = -5 * ones(size(x));
    patch([x, fliplr(x)], [y_deep_bound, fliplr(y_bot)], col_deep, 'EdgeColor', 'k', 'LineWidth', lw_border);
    patch([x, fliplr(x)], [y_therm, fliplr(y_deep_bound)], col_sub, 'EdgeColor', 'k', 'LineWidth', lw_border);
    plot(x, y_deep_bound, 'g', 'LineWidth', 2);
    y_200 = -2.93;
    plot([-5.5, 5.5], [y_200, y_200], 'g--', 'LineWidth', 2);
    y_th_0 = interp1(x, y_therm, 0); 
    y_db_0 = interp1(x, y_deep_bound, 0);  
    plot([0, 0], [y_db_0, y_th_0], 'k--', 'LineWidth', 1.5);
    idx_mid = find(x >= 0, 1);
    x_left = x(1:idx_mid);
    y_surf_left = y_surf(1:idx_mid);
    y_therm_left = y_therm(1:idx_mid);
    patch([x_left, 0, fliplr(x_left), -5.5], ...
          [y_surf_left, y_therm_left(end), fliplr(y_therm_left), y_surf_left(1)], ...
          col_t1, 'EdgeColor', 'k', 'LineWidth', lw_border);
    x_right = x(idx_mid:end);
    y_surf_right = y_surf(idx_mid:end);
    y_therm_right = y_therm(idx_mid:end);
    patch([x_right, 5.5, fliplr(x_right), 0], ...
          [y_surf_right, y_therm_right(end), fliplr(y_therm_right), y_surf_right(1)], ...
          col_t2, 'EdgeColor', 'k', 'LineWidth', lw_border);
    lw_h_arrow = 1.5; 
    x_pos_h1 = -3.5;
    x_pos_h2 = 3.5;
    y_func = @(xi) -2.42 + 1.55 * tanh(0.4 * (xi - 0.85));
    draw_double_arrow(x_pos_h1, y_200, y_func(x_pos_h1), 'g', lw_h_arrow);
    draw_double_arrow(x_pos_h2, y_200, y_func(x_pos_h2), 'g', lw_h_arrow);
    text(-3.2, 0.6, 'T1', 'Horiz', 'center', 'FontSize', 24, 'FontWeight', 'bold');
    text(3.2, 0.6, 'T2', 'Horiz', 'center', 'FontSize', 24, 'FontWeight', 'bold');
    
    text(-2.5, -3.3, "h1'", 'Horiz', 'right', 'FontSize', 24, 'FontAngle', 'italic', 'FontWeight', 'bold');
    text(2.6, -2, "h2'", 'Horiz', 'left', 'FontSize', 24, 'FontAngle', 'italic', 'FontWeight', 'bold');

    p_tl = [-2.5, 0.6]; % Top Left
    p_tr = [ 2.5, 0.6]; % Top Right
    p_bl = [-2.5, -1.8]; % Bottom Left (near h1)
    p_br = [ 2.5, -0.9]; % Bottom Right (near h2)
    
    lw_arr = 4;
    
    % 箭头1: 表面 T2 -> T1 (向左)
    draw_straight_arrow(p_tr(1), p_tr(2), p_tl(1), p_tl(2), 'k', lw_arr);
    
    % 箭头2: 下沉 T1 -> h1 (向下)
    draw_straight_arrow(p_tl(1), p_tl(2), p_bl(1), p_bl(2), 'k', lw_arr);
    
    % 箭头3: 温跃层 h1 -> h2 (向右斜上)
    draw_straight_arrow(p_bl(1), p_bl(2), p_br(1), p_br(2), 'k', lw_arr);
    
    % 箭头4: 上升 h2 -> T2 (向上)
    draw_straight_arrow(p_br(1), p_br(2), p_tr(1), p_tr(2), 'k', lw_arr);
    
    % 6.2 风向虚线箭头 (Wind Stress)
    y_wind = 1.8;
    plot([3.5, -3.5], [y_wind, y_wind], 'k--', 'LineWidth', 2); % 虚线杆
    % 虚线箭头的头 (实线)
    patch([-3.5, -3.0, -3.0], [y_wind, y_wind+0.15, y_wind-0.15], 'k', 'EdgeColor', 'none');
    % =========================
% 左侧 Depth 坐标轴
% =========================

x_axis = -5.5;  % 稍微放在图外一点
y_top = 1.2;
y_bot = -5;

% 主轴线
plot([x_axis, x_axis], [y_bot, y_top], 'k', 'LineWidth', 1.5);

% 刻度（对应深度）
depth_vals = [0,50, 100, 150, 200, 250,300];
y_ticks = linspace(y_top, y_bot, length(depth_vals));

for i = 1:length(depth_vals)
    % 刻度线
    plot([x_axis-0.2, x_axis], [y_ticks(i), y_ticks(i)], 'k', 'LineWidth', 1.5);
    
    % 数字
    text(x_axis-0.4, y_ticks(i), num2str(depth_vals(i)), ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 20);
end

% 轴标签
text(x_axis-1.5, (y_top+y_bot)/2, 'Depth [m]', ...
    'Rotation', 90, ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 20);
end

% =========================================================================
%  通用绘图函数 (箭头等)
% =========================================================================
function draw_straight_arrow(x1, y1, x2, y2, color, lw)
    % 画箭杆
    plot([x1, x2], [y1, y2], 'Color', color, 'LineWidth', lw);
    
    % 画箭头
    angle = atan2(y2-y1, x2-x1);
    head_len = 0.6;
    head_ang = 25 * pi/180;
    
    p1x = x2 - head_len * cos(angle - head_ang);
    p1y = y2 - head_len * sin(angle - head_ang);
    p2x = x2 - head_len * cos(angle + head_ang);
    p2y = y2 - head_len * sin(angle + head_ang);
    
    patch([x2, p1x, p2x], [y2, p1y, p2y], color, 'EdgeColor', 'none');
end
function draw_vertical_temp_profile(legend_mode)
    if nargin < 1
        legend_mode = 'plain';
    end
    % --- 位置定义 ---
    x_box_left = 5.7;
    x_box_right = 9.2;
    y_surface = 1.2;
    y_mixed_base = 0.167;
    y_top = y_surface;
    y_bot = -5;
    
    % --- 绘制坐标框 ---
    rectangle('Position', [x_box_left, y_bot, (x_box_right-x_box_left), (y_top-y_bot)], ...
              'EdgeColor', 'k', 'LineWidth', 1.2);
    grid_y = linspace(y_top, y_bot, 6); % 对应 50, ..., 300m
    grid_y = linspace(y_top, y_bot, 7); % 0, 50, ..., 300 m
    for gy = grid_y(2:end-1)
        plot([x_box_left, x_box_right], [gy, gy], ':', 'Color', [0.8 0.8 0.8]);
    end
    plot([x_box_left, x_box_right], [y_mixed_base, y_mixed_base], ...
        '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0);

    % --- 温度廓线计算 (公式 4-5) ---
    % 参数设置 (根据图3虚线效果拟合)
    Te = 29.5; 
    Tb = 17.5; 
    z0 = 75; % 跃层中心深度
    H_star = 65; % 跃层强度/平滑度
    
    % 定义深度 z (-50 到 -300)
    z_abs = linspace(0, -250, 100);
    % Phi(z) 公式
    phi_z = Te - ((Te - Tb)/2) * (1 - tanh((z_abs-50 + z0) / H_star));
    target_mixed_temp = 27;
    phi_z = phi_z + (target_mixed_temp - phi_z(1));
    
    % --- 映射到绘图坐标系 ---
    % Temp 12~28 degC maps to x_box_left ~ x_box_right.
    temp_min = 12; temp_max = 28;
    plot_x = x_box_left + (phi_z - temp_min) / (temp_max - temp_min) * (x_box_right - x_box_left);
    % Depth -50~-300 映射到 y_top ~ y_bot
    plot_y = y_mixed_base + (z_abs / 250) * (y_mixed_base - y_bot);
    
    % 绘制非线性廓线 (nonlinear)
    plot(plot_x, plot_y, 'b-', 'LineWidth', 2.5);

    % Linearized profile.  It uses the same reference temperature at h2 = 0
    % and the linear model lapse rate (Te - Tb)/(2*H_star).
    phi_ref = phi_z(1);
    z_lin_end = max(min(z_abs), (temp_min - phi_ref) * 2 * H_star / (Te - Tb));
    z_lin = linspace(0, z_lin_end, 50);
    phi_lin = phi_ref + z_lin * (Te - Tb) / (2 * H_star);
    lin_x = x_box_left + (phi_lin - temp_min) / (temp_max - temp_min) * (x_box_right - x_box_left);
    lin_y = y_mixed_base + (z_lin / 250) * (y_mixed_base - y_bot);
    plot(lin_x, lin_y, 'r-', 'LineWidth', 2.0);

    % Extend both profiles upward through the mixed layer. The mixed layer
    % is vertically isothermal, so the nonlinear and linear profiles share
    % the same temperature there. Blue is drawn wider and red narrower so
    % the overlap remains visible without implying two different values.
    mixed_temp_x = x_box_left + (phi_ref - temp_min) / (temp_max - temp_min) * (x_box_right - x_box_left);
    mixed_y = [y_surface, y_mixed_base];
    plot([mixed_temp_x, mixed_temp_x], mixed_y, 'b-', 'LineWidth', 4.0);
    plot([mixed_temp_x, mixed_temp_x], mixed_y, 'r-', 'LineWidth', 1.8);

    % Manual legend placed above the profile panel.
    legend_x0 = x_box_left + 0.12;
    legend_x1 = legend_x0 + 0.38;
    legend_tx = legend_x1 + 0.12;
    if strcmp(legend_mode, 'formula')
        label_font = 11;
        formula_font = 8.6;
        legend_y_nl = y_top + 2.78;
        legend_y_nl_eq1 = legend_y_nl - 0.42;
        legend_y_nl_eq2 = legend_y_nl_eq1 - 0.33;
        legend_y_lin = legend_y_nl_eq2 - 0.54;
        legend_y_lin_eq = legend_y_lin - 0.42;
        legend_y_lin_eq2 = legend_y_lin_eq - 0.33;

        plot([legend_x0, legend_x1], [legend_y_nl, legend_y_nl], 'b-', 'LineWidth', 2.5);
        text(legend_tx, legend_y_nl, 'nonlinear', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', label_font);
        text(legend_tx, legend_y_nl_eq1, 'Phi(z) = Te - (Te - Tb)/2', ...
            'Interpreter', 'none', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', 'FontSize', formula_font);
        text(legend_tx, legend_y_nl_eq2, '[1 - tanh((z - H1 + z0)/H*)]', ...
            'Interpreter', 'none', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', 'FontSize', formula_font);

        plot([legend_x0, legend_x1], [legend_y_lin, legend_y_lin], 'r-', 'LineWidth', 2.0);
        text(legend_tx, legend_y_lin, 'linear', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', label_font);
        text(legend_tx, legend_y_lin_eq, 'Phi_L(z) = Phi(0) +', ...
            'Interpreter', 'none', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', 'FontSize', formula_font);
        text(legend_tx, legend_y_lin_eq2, '(Te - Tb)z/(2H*)', ...
            'Interpreter', 'none', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', 'FontSize', formula_font);
    else
        legend_font = 12;
        legend_y_nl = y_top + 2.05;
        legend_y_nl2 = legend_y_nl - 0.43;
        legend_y_lin = legend_y_nl - 1.12;
        legend_y_lin2 = legend_y_lin - 0.43;

        plot([legend_x0, legend_x1], [legend_y_nl, legend_y_nl], 'b-', 'LineWidth', 2.5);
        text(legend_tx, legend_y_nl, 'nonlinear profile', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', legend_font);
        text(legend_tx, legend_y_nl2, 'used in original model', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', legend_font);

        plot([legend_x0, legend_x1], [legend_y_lin, legend_y_lin], 'r-', 'LineWidth', 2.0);
        text(legend_tx, legend_y_lin, 'linear profile', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', legend_font);
        text(legend_tx, legend_y_lin2, 'in nondimensional model', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', legend_font);
    end
    
    % --- 坐标轴标签 ---
    text((x_box_left + x_box_right)/2, y_bot - 1.2, 'Tsub [\circC]', 'Horiz', 'center', 'FontSize', 20);
    
    % 温度刻度 [12, 16, 20, 24, 28]
    temp_labels = [12, 16, 20, 24, 28];
    for tl = temp_labels
        lx = x_box_left + (tl - temp_min) / (temp_max - temp_min) * (x_box_right - x_box_left);
        plot([lx, lx], [y_bot, y_bot-0.15], 'k', 'LineWidth', 1.2);
        text(lx, y_bot - 0.5, num2str(tl), 'Horiz', 'center', 'FontSize', 18);
    end
    

end
% --- 辅助函数：绘制双向箭头 ---
function draw_double_arrow(x, y1, y2, color, lw)
    % 画主干线
    plot([x, x], [y1, y2], 'Color', color, 'LineWidth', lw);
    
    head_len = 0.25;
    head_ang = 30 * pi/180;
    
    % 上箭头
    patch([x, x-head_len*sin(head_ang), x+head_len*sin(head_ang)], ...
          [max(y1,y2), max(y1,y2)-head_len*cos(head_ang), max(y1,y2)-head_len*cos(head_ang)], ...
          color, 'EdgeColor', 'none');
    % 下箭头
    patch([x, x-head_len*sin(head_ang), x+head_len*sin(head_ang)], ...
          [min(y1,y2), min(y1,y2)+head_len*cos(head_ang), min(y1,y2)+head_len*cos(head_ang)], ...
          color, 'EdgeColor', 'none');
end
