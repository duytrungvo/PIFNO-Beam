try:
    import wandb
except ImportError:
    wandb = None

import matplotlib.pyplot as plt

def plot_pred(data_config, test_x, test_y, pred_y, error_index):
    for i in range(3):
        # key = np.random.randint(0, data_config['n_sample'])
        key = error_index[-1-i]
        x_plot = test_x[key]
        y_true_plot = test_y[key]
        y_pred_plot = pred_y[key]

        fig = plt.figure(figsize=(10, 12))
        if data_config['out_dim'] == 1:
            plt.subplot(3, 1, 1)
            plt.plot(x_plot[:, -1], x_plot[:, 0] / data_config['I0'])
            plt.xlabel('$x$')
            plt.ylabel('$I/I_0$')
            plt.title(f'Input $I(x)$')
            plt.xlim([0, 1])
            plt.ylim([0, 1])

            plt.subplot(3, 1, 2)
            plt.plot(x_plot[:, -1], x_plot[:, 1] / data_config['q0'])
            plt.xlabel('$x$')
            plt.ylabel('$q/q_0$')
            plt.title(f'Input $q(x)$')
            plt.xlim([0, 1])
            plt.ylim([-0.1, 1.1])

            plt.subplot(3, 1, 3)
            plt.plot(x_plot[:, -1], y_pred_plot[:, 0], 'r', label='predict sol')
            plt.plot(x_plot[:, -1], y_true_plot[:, 0], 'b', label='exact sol')
            plt.xlabel('$x$')
            plt.ylabel(r'$w$')
            # plt.ylim([0, 1])
            plt.legend()
            plt.grid(visible=True)
            plt.title(f'Predict and exact $w(x)$')
            plt.tight_layout()
        if data_config['out_dim'] == 2:
           plt.subplot(4, 1, 1)
           plt.plot(x_plot[:, -1], x_plot[:, 0] / data_config['I0'])
           plt.xlabel('$x$')
           plt.ylabel('$I/I_0$')
           plt.title(f'Input $I(x)$')
           plt.xlim([0, 1])
           # plt.ylim([0, 2])

           plt.subplot(4, 1, 2)
           plt.plot(x_plot[:, -1], x_plot[:, 1] / data_config['q0'])
           plt.xlabel('$x$')
           plt.ylabel('$q/q_0$')
           plt.title(f'Input $q(x)$')
           plt.xlim([0, 1])
           plt.ylim([-0.1, 1.1])

           plt.subplot(4, 1, 3)
           plt.plot(x_plot[:, -1], y_pred_plot[:, 0], 'r', label='predict sol')
           plt.plot(x_plot[:, -1], y_true_plot[:, 0], 'b', label='exact sol')
           plt.xlabel('$x$')
           plt.ylabel(r'$w$')
           # plt.ylim([0, 1])
           plt.legend()
           plt.grid(visible=True)
           plt.title(f'Predict and exact $w(x)$')

           plt.subplot(4, 1, 4)
           plt.plot(x_plot[:, -1], y_pred_plot[:, 1], 'r', label='predict sol')
           plt.plot(x_plot[:, -1], y_true_plot[:, 1], 'b', label='exact sol')
           plt.xlabel('$x$')
           plt.ylabel(r'$M$')
           # plt.ylim([0, 1])
           plt.legend()
           plt.grid(visible=True)
           plt.title(f'Predict and exact $M(x)$')
           plt.tight_layout()

        plt.show()

