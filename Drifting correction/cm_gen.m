function color_map_list = cm_gen(num_plot)
cm = jet;
idx = 1:(num_plot-1);
idx = floor((idx./num_plot)*256);
idx = [1 idx];
color_map_list = cm(idx,:);
end
