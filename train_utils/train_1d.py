import torch
from .losses import LpLoss
from tqdm import tqdm
from .utils import save_checkpoint, save_loss, shape_function

try:
    import wandb
except ImportError:
    wandb = None

def train_EB_Beam(model,
                 train_loader,
                 optimizer,
                 scheduler,
                 device,
                 config,
                 pino_loss,
                 log=False,
                 use_tqdm=True):

    model.to(device)
    data_weight = config['train']['xy_loss']
    f_weight = config['train']['f_loss']
    bc_weight_l = config['train']['bc_loss_l']
    bc_weight_r = config['train']['bc_loss_r']
    myloss = LpLoss(size_average=True)
    pbar = range(config['train']['epochs'])
    if use_tqdm:
        pbar = tqdm(pbar, dynamic_ncols=True, smoothing=0.1)
    train_loss_epoch = torch.zeros(config['train']['epochs'], 6, device=device)

    for e in pbar:
        model.train()
        physic_mse = 0.0
        physic1_mse = 0.0
        bc_l_mse = 0.0
        bc_r_mse = 0.0
        data_l2 = 0.0
        train_loss = 0.0
        for x, y in train_loader:
            x, y = x.to(device), y.to(device)
            out = model(x).reshape(y.shape)

            # loss
            data_loss = myloss(out, y).to(device)
            f_loss1, f_loss2, bc_loss_l, bc_loss_r = pino_loss(config['data'], x, out)
            f_loss = f_loss1 + f_loss2
            total_loss = f_weight * f_loss + bc_weight_l * bc_loss_l + bc_weight_r * bc_loss_r + data_weight * data_loss

            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()

            physic_mse += f_loss.item()
            physic1_mse += f_loss1.item()
            bc_l_mse += bc_loss_l.item()
            bc_r_mse += bc_loss_r.item()
            data_l2 += data_loss.item()
            train_loss += total_loss.item()

        scheduler.step()
        physic_mse /= len(train_loader)
        physic1_mse /= len(train_loader)
        bc_l_mse /= len(train_loader)
        bc_r_mse /= len(train_loader)
        data_l2 /= len(train_loader)
        train_loss /= len(train_loader)

        if use_tqdm:
            pbar.set_description(
                (
                    f'Epoch {e}, train loss: {train_loss:.5E}; '
                    f'train f error: {physic_mse:.5E}; '
                    f'train f1 error: {physic1_mse:.5E}; '
                    f'train bc left error: {bc_l_mse:.5E}; '
                    f'train bc right error: {bc_r_mse:.5E}; '
                    f'data l2 error: {data_l2:.5E}'
                )
            )
        if wandb and log:
            wandb.log(
                {
                    'Train f error': physic_mse,
                    'Train bc left error': bc_l_mse,
                    'Train bc right error': bc_r_mse,
                    'Train L2 error': data_l2,
                    'Train loss': train_loss,
                }
            )

        if e % 100 == 0:
            save_checkpoint(config['train']['save_dir'],
                            config['train']['save_name'].replace('.pt', f'_{e}.pt'),
                            model, optimizer)

        train_loss_epoch[e, :] = torch.tensor([train_loss, physic_mse, bc_l_mse, bc_r_mse, data_l2, physic1_mse],
                                              device=device)

    save_loss(config['train']['save_dir'],
              config['train']['loss_save_name'],
              train_loss_epoch)

    save_checkpoint(config['train']['save_dir'],
                    config['train']['save_name'],
                    model, optimizer)

# train Euler Bernoulli beam using boundary shape function
def train_EB_Beam_bsf(model,
                 train_loader,
                 optimizer,
                 scheduler,
                 device,
                 config,
                 pino_loss,
                 log=False,
                 use_tqdm=True):

    model.to(device)
    data_weight = config['train']['xy_loss']
    f_weight = config['train']['f_loss']
    bc_weight_l = config['train']['bc_loss_l']
    bc_weight_r = config['train']['bc_loss_r']

    myloss = LpLoss(size_average=True)
    pbar = range(config['train']['epochs'])
    if use_tqdm:
        pbar = tqdm(pbar, dynamic_ncols=True, smoothing=0.1)
    train_loss_epoch = torch.zeros(config['train']['epochs'], 6, device=device)

    batchsize = train_loader.batch_size
    nx = train_loader.dataset[0][0].size(0)
    x_loss = train_loader.dataset[0][0][:, -1].repeat(batchsize).reshape(batchsize, nx)
    bc = shape_function(config['data']['BC'], x_loss, config['data']['L']).to(device)
    for e in pbar:
        model.train()
        physic_mse = 0.0
        physic1_mse = 0.0
        bc_l_mse = 0.0
        bc_r_mse = 0.0
        data_l2 = 0.0
        train_loss = 0.0
        for x, y in train_loader:
            x, y = x.to(device), y.to(device)
            out = model(x).reshape(y.shape)

            # loss
            data_loss = myloss(out, y).to(device)
            f_loss1, f_loss2, bc_loss_l, bc_loss_r = pino_loss(config['data'], x, out, bc)
            f_loss = f_loss1 + f_loss2
            total_loss = f_weight * f_loss + bc_weight_l * bc_loss_l + bc_weight_r * bc_loss_r + data_weight * data_loss

            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()

            physic_mse += f_loss.item()
            physic1_mse += f_loss1.item()
            bc_l_mse += bc_loss_l.item()
            bc_r_mse += bc_loss_r.item()
            data_l2 += data_loss.item()
            train_loss += total_loss.item()

        scheduler.step()
        physic_mse /= len(train_loader)
        physic1_mse /= len(train_loader)
        bc_l_mse /= len(train_loader)
        bc_r_mse /= len(train_loader)
        data_l2 /= len(train_loader)
        train_loss /= len(train_loader)

        if use_tqdm:
            pbar.set_description(
                (
                    f'Epoch {e}, train loss: {train_loss:.5E}; '
                    f'train f error: {physic_mse:.5E}; '
                    f'train f1 error: {physic1_mse:.5E}; '
                    f'train bc left error: {bc_l_mse:.5E}; '
                    f'train bc right error: {bc_r_mse:.5E}; '
                    f'data l2 error: {data_l2:.5E}'
                )
            )
        if wandb and log:
            wandb.log(
                {
                    'Train f error': physic_mse,
                    'Train bc left error': bc_l_mse,
                    'Train bc right error': bc_r_mse,
                    'Train L2 error': data_l2,
                    'Train loss': train_loss,
                }
            )

        if e % 100 == 0:
            save_checkpoint(config['train']['save_dir'],
                            config['train']['save_name'].replace('.pt', f'_{e}.pt'),
                            model, optimizer)

        train_loss_epoch[e, :] = torch.tensor([train_loss, physic_mse, bc_l_mse, bc_r_mse, data_l2, physic1_mse],
                                              device=device)

    save_loss(config['train']['save_dir'],
              config['train']['loss_save_name'],
              train_loss_epoch)

    save_checkpoint(config['train']['save_dir'],
                    config['train']['save_name'],
                    model, optimizer)
