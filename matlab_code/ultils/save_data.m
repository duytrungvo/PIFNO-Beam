function save_data(folder_name, file_name, BC, nx, seed_number, input, output, x, alpha, beta, q0, E, I0, h0, L)

filename = [file_name, '_', BC,'_', num2str(nx),'_', 'rng', num2str(seed_number),'.mat'];
filepath = fullfile(folder_name, filename);
save(filepath,...
    "input", "output", "x", "alpha", "beta", "q0", "E", "I0", "h0", "L")

end