from argparse import ArgumentParser
import yaml
import torch
from models import FNO1d
from train_utils import Adam
from train_utils.datasets import Loader_1D
from train_utils.train_1d import train_EB_Beam
from train_utils.losses import LpLoss, FDM_Order4_Euler_Bernoulli_Beam
from train_utils.plot_test import plot_pred
import numpy as np
import matplotlib.pyplot as plt

torch.manual_seed(0)
np.random.seed(0)
torch.cuda.manual_seed(0)
torch.cuda.manual_seed_all(0)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

def run(config, device):

    data_config = config['data']

    dataset = Loader_1D(data_config['datapath'],
                        nx=data_config['nx'],
                        sub=data_config['sub'],
                        in_dim=data_config['in_dim'],
                        out_dim=data_config['out_dim'])

    train_loader = dataset.make_loader(n_sample=data_config['n_sample'],
                                       batch_size=config['train']['batchsize'],
                                       start=data_config['offset'])

    # define model
    model_config = config['model']
    if model_config['name'] == 'fno':
        model = FNO1d(modes=model_config['modes'],
                      fc_dim=model_config['fc_dim'],
                      layers=model_config['layers'],
                      in_dim=data_config['in_dim'],
                      out_dim=data_config['out_dim'],
                      act=model_config['act'])

    # train
    train_config = config['train']
    optimizer = Adam(model.parameters(), betas=(0.9, 0.999),
                     lr=train_config['base_lr'])
    scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer,
                                                     milestones=train_config['milestones'],
                                                     gamma=train_config['scheduler_gamma'])

    pino_loss = FDM_Order4_Euler_Bernoulli_Beam
    train_EB_Beam(model,
             train_loader,
             optimizer,
             scheduler,
             device,
             config,
             pino_loss=pino_loss,
             log=False,
             use_tqdm=True)

    path = config['train']['save_dir']
    ckpt_dir = 'checkpoints/%s/' % path
    loss_history = np.loadtxt(ckpt_dir+'train_loss_history.txt')
    plt.semilogy(range(len(loss_history)), loss_history[:, :5])
    plt.xlabel('Epochs')
    plt.ylabel('$Loss$')
    plt.tight_layout()
    plt.legend(['train_loss', 'pde_mse', 'bc_l_mse', 'bc_r_mse', 'data_l2'])
    plt.show()

def test(config, device):

    data_config = config['data']
    dataset = Loader_1D(data_config['datapath'],
                        nx=data_config['nx'],
                        sub=data_config['sub'],
                        in_dim=data_config['in_dim'],
                        out_dim=data_config['out_dim'])
    data_loader = dataset.make_loader(n_sample=data_config['n_sample'],
                                      batch_size=config['test']['batchsize'],
                                      start=data_config['offset'])

    model_config = config['model']
    if model_config['name'] == 'fno':
        model = FNO1d(modes=model_config['modes'],
                      fc_dim=model_config['fc_dim'],
                      layers=model_config['layers'],
                      in_dim=data_config['in_dim'],
                      out_dim=data_config['out_dim'],
                      act=model_config['act']).to(device)

    # Load from checkpoint
    if 'ckpt' in config['test']:
        ckpt_path = config['test']['ckpt']
        ckpt = torch.load(ckpt_path)
        model.load_state_dict(ckpt['model'])
        print('Weights loaded from %s' % ckpt_path)

    myloss = LpLoss(size_average=True)
    model.eval()
    nsample = data_config['n_sample']
    in_dim = data_config['in_dim']
    out_dim = data_config['out_dim']
    s = int(np.ceil(data_config['nx'] / data_config['sub']))
    test_x = np.zeros((nsample, s, in_dim))
    preds_y = np.zeros((nsample, s, out_dim))
    test_y = np.zeros((nsample, s, out_dim))
    test_err = []
    test_err_w = []
    test_err_m = []
    with torch.no_grad():
        for i, data in enumerate(data_loader):
            data_x, data_y = data
            data_x, data_y = data_x.to(device), data_y.to(device)
            pred_y = model(data_x).reshape(data_y.shape)
            data_loss = myloss(pred_y, data_y).to(device)
            data_loss_w = myloss(pred_y[:, :, 0], data_y[:, :, 0])
            data_loss_m = np.nan
            test_err.append(data_loss.item())
            test_err_w.append(data_loss_w.item())
            test_err_m.append(data_loss_m)
            test_x[i] = data_x.cpu().numpy()
            test_y[i] = data_y.cpu().numpy()
            preds_y[i] = pred_y.cpu().numpy()

    mean_err = np.mean(test_err)
    std_err = np.std(test_err, ddof=1) / np.sqrt(len(test_err))
    mean_err_w = np.mean(test_err_w)
    std_err_w = np.std(test_err_w, ddof=1) / np.sqrt(len(test_err))
    mean_err_m = np.mean(test_err_m)
    std_err_m = np.std(test_err_m, ddof=1) / np.sqrt(len(test_err))

    print(f'==Test for {nsample} samples==')
    print(f'==Averaged relative L2 error mean(w & M): {mean_err}, std error: {std_err}==')
    print(f'==Averaged relative L2 error mean(w): {mean_err_w}, std error: {std_err_w}==')
    print(f'==Averaged relative L2 error mean(M): {mean_err_m}, std error: {std_err_m}==')

    # err_idx = np.argsort(test_err)
    err_idx = np.arange(0, nsample)
    np.random.shuffle(err_idx)
    plot_pred(data_config, test_x, test_y, preds_y, err_idx)

if __name__ == '__main__':

    parser = ArgumentParser(description='Basic paser')
    parser.add_argument('--config_path', type=str, help='Path to the configuration file')
    parser.add_argument('--log', action='store_true', help='Turn on the wandb')
    parser.add_argument('--mode', type=str, help='train or test')
    args = parser.parse_args()

    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

    # read data
    config_file = args.config_path
    with open(config_file, 'r') as stream:
        config = yaml.load(stream, yaml.FullLoader)
    if args.mode == 'train':
        run(config, device)
    else:
        test(config, device)

    print('Done')