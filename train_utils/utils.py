import os
import numpy as np
import torch

def zero_grad(params):
    '''
    set grad field to 0
    '''
    if isinstance(params, torch.Tensor):
        if params.grad is not None:
            params.grad.zero_()
    else:
        for p in params:
            if p.grad is not None:
                p.grad.zero_()

def save_checkpoint(path, name, model, optimizer=None):
    ckpt_dir = 'checkpoints/%s/' % path
    if not os.path.exists(ckpt_dir):
        os.makedirs(ckpt_dir)
    try:
        model_state_dict = model.module.state_dict()
    except AttributeError:
        model_state_dict = model.state_dict()

    if optimizer is not None:
        optim_dict = optimizer.state_dict()
    else:
        optim_dict = 0.0

    torch.save({
        'model': model_state_dict,
        'optim': optim_dict
    }, ckpt_dir + name)
    print('Checkpoint is saved at %s' % ckpt_dir + name)

def save_loss(path, name, loss):
    ckpt_dir = f'checkpoints/{path}/'
    os.makedirs(ckpt_dir, exist_ok=True)

    if isinstance(loss, torch.Tensor):
        loss = loss.detach().cpu().numpy()

    np.savetxt(ckpt_dir + name, loss)
    print('Training lost is saved at %s' % ckpt_dir + name)

def test_func_disp(BC):
    if BC == 'CC':
        return lambda x, l: x * (l - x)                                                 #psi1
        # return lambda x, l: - x ** 2 / l + x ** 3 / l ** 2                            #psi2
        # return lambda x, l: x ** 2 / l ** 2 - 2 * x ** 3 / l ** 3 + x ** 4 / l ** 4   #psi3
        # return lambda x, l: torch.cos(2 * torch.pi * x / l) - 1                       #psi4
    if BC == 'CS':
        return lambda x, l: x * (l - x)
        # return lambda x, l: - x ** 2 / l + x ** 3 / l ** 2
    if BC == 'CF':
        return lambda x, l: x * (l - x)
    if BC == 'SS':
        return lambda x, l: x * (l - x)

def test_func_moment(BC):
    if BC == 'CC':
        return lambda x, l: 1
    if BC == 'CS':
        return lambda x, l: l - x
        # return lambda x, l: l - x
    if BC == 'CF':
        return lambda x, l: x * (l - x)
    if BC == 'SS':
        return lambda x, l: x * (l - x)

def shape_function(BC, x, L):
    bc = torch.zeros((x.size(0), x.size(1), 16))
    if BC == 'CF':
        bc[:, :, 0] = 1
        bc[:, :, 1] = x
        bc[:, :, 6] = 1
        bc[:, :, 7] = x - L
    if BC == 'SS':
        bc[:, :, 0] = 1 - x / L
        bc[:, :, 2] = x / L
        bc[:, :, 4] = 1 - x / L
        bc[:, :, 6] = x / L
    if BC == 'CS':
        bc[:, :, 0] = 1 - x ** 2 / L ** 2
        bc[:, :, 1] = x - x ** 2 / L
        bc[:, :, 2] = x ** 2 / L ** 2
        bc[:, :, 6] = 1
        bc[:, :, 8] = - 2 / L ** 2
        bc[:, :, 9] = - 2 / L
        bc[:, :, 10] = 2 / L ** 2
    if BC == 'CC':
        bc[:, :, 0] = 1 - 3 * x ** 2 / L ** 2 + 2 * x ** 3 / L ** 3
        bc[:, :, 1] = x - 2 * x ** 2 / L + x ** 3 / L ** 2
        bc[:, :, 2] = 3 * x ** 2 / L ** 2 - 2 * x ** 3 / L ** 3
        bc[:, :, 3] = - x ** 2 / L + x ** 3 / L ** 2
        bc[:, :, 8] = - 6 / L ** 2 + 12 * x / L ** 3
        bc[:, :, 9] = - 4 / L + 6 * x / L ** 2
        bc[:, :, 10] = 6 / L ** 2 - 12 * x / L ** 3
        bc[:, :, 11] = - 2 / L + 6 * x / L ** 2

    return bc

def boundary_function(w, m, bc, nx, dx, BC):
    batchsize = w.size(0)

    b1, c1, d2b1dx2, d2c1dx2 = bc[:, :, 0], bc[:, :, 4], bc[:, :, 8], bc[:, :, 12]
    b2, c2, d2b2dx2, d2c2dx2 = bc[:, :, 1], bc[:, :, 5], bc[:, :, 9], bc[:, :, 13]
    b3, c3, d2b3dx2, d2c3dx2 = bc[:, :, 2], bc[:, :, 6], bc[:, :, 10], bc[:, :, 14]
    b4, c4, d2b4dx2, d2c4dx2 = bc[:, :, 3], bc[:, :, 7], bc[:, :, 11], bc[:, :, 15]

    w0 = torch.repeat_interleave(w[:, 0], nx, dim=0).reshape((batchsize, nx))
    wL = torch.repeat_interleave(w[:, -1], nx, dim=0).reshape((batchsize, nx))

    m0 = torch.repeat_interleave(m[:, 0], nx, dim=0).reshape((batchsize, nx))
    mL = torch.repeat_interleave(m[:, -1], nx, dim=0).reshape((batchsize, nx))

    dw0 = (-1.5 * w[:, 0] + 2 * w[:, 1] - 0.5 * w[:, 2]) / dx
    dwdx0 = torch.repeat_interleave(dw0, nx, dim=0).reshape((batchsize, nx))

    dwL = (0.5 * w[:, -3] - 2 * w[:, -2] + 1.5 * w[:, -1]) / dx
    dwdxL = torch.repeat_interleave(dwL, nx, dim=0).reshape((batchsize, nx))

    dm0 = (-1.5 * m[:, 0] + 2 * m[:, 1] - 0.5 * m[:, 2]) / dx
    dmdx0 = torch.repeat_interleave(dm0, nx, dim=0).reshape((batchsize, nx))

    dmL = (0.5 * m[:, -3] - 2 * m[:, -2] + 1.5 * m[:, -1]) / dx
    dmdxL = torch.repeat_interleave(dmL, nx, dim=0).reshape((batchsize, nx))

    G1 = b1 * w0 + b2 * dwdx0 + b3 * wL + b4 * dwdxL
    G2 = c1 * m0 + c2 * dmdx0 + c3 * mL + c4 * dmdxL
    d2G1dx2 = d2b1dx2 * w0 + d2b2dx2 * dwdx0 + d2b3dx2 * wL + d2b4dx2 * dwdxL
    d2G2dx2 = d2c1dx2 * m0 + d2c2dx2 * dmdx0 + d2c3dx2 * mL + d2c4dx2 * dmdxL

    if BC == 'CF':
        boundary_l = torch.stack((w0, dwdx0), 1)
        boundary_r = torch.stack((mL, dmdxL), 1)
    if BC == 'CS':
        boundary_l = torch.stack((w0, dwdx0), 1)
        boundary_r = torch.stack((wL, mL), 1)
    if BC == 'CC':
        boundary_l = torch.stack((w0, dwdx0), 1)
        boundary_r = torch.stack((wL, dwdxL), 1)
    if BC == 'SS':
        boundary_l = torch.stack((w0, m0), 1)
        boundary_r = torch.stack((wL, mL), 1)

    return G1, G2, d2G1dx2, d2G2dx2, boundary_l, boundary_r