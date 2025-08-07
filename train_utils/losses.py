import torch
import torch.nn.functional as F
from train_utils.utils import boundary_function

def pino_loss_reduced_order2_1d(func):
    def inner(*args, **kwargs):
        Du1, Du2, boundary_l, boundary_r = func(*args, **kwargs)
        f1 = torch.zeros(Du1.shape, device=args[2].device)
        f_loss1 = F.mse_loss(Du1, f1)
        f2 = torch.zeros(Du2.shape, device=args[2].device)
        f_loss2 = F.mse_loss(Du2, f2)

        loss_boundary_l = F.mse_loss(boundary_l, torch.zeros(boundary_l.shape, device=args[2].device))
        loss_boundary_r = F.mse_loss(boundary_r, torch.zeros(boundary_r.shape, device=args[2].device))

        return f_loss1, f_loss2, loss_boundary_l, loss_boundary_r
    return inner

class LpLoss(object):
    '''
    loss function with rel/abs Lp loss
    '''
    def __init__(self, d=2, p=2, size_average=True, reduction=True):
        super(LpLoss, self).__init__()

        #Dimension and Lp-norm type are postive
        assert d > 0 and p > 0

        self.d = d
        self.p = p
        self.reduction = reduction
        self.size_average = size_average

    def abs(self, x, y):
        num_examples = x.size()[0]

        #Assume uniform mesh
        h = 1.0 / (x.size()[1] - 1.0)

        all_norms = (h**(self.d/self.p))*torch.norm(x.view(num_examples,-1) - y.view(num_examples,-1), self.p, 1)

        if self.reduction:
            if self.size_average:
                return torch.mean(all_norms)
            else:
                return torch.sum(all_norms)

        return all_norms

    def rel(self, x, y):
        num_examples = x.size()[0]

        diff_norms = torch.norm(x.reshape(num_examples,-1) - y.reshape(num_examples,-1), self.p, 1)
        y_norms = torch.norm(y.reshape(num_examples,-1), self.p, 1)

        if self.reduction:
            if self.size_average:
                return torch.mean(diff_norms/y_norms)
            else:
                return torch.sum(diff_norms/y_norms)

        return diff_norms/y_norms

    def __call__(self, x, y):
        return self.rel(x, y)

@pino_loss_reduced_order2_1d
def FDM_Order4_Euler_Bernoulli_Beam(config_data, a, u):
    batchsize = u.size(0)
    nx = u.size(1)
    dx = 1 / (nx - 1)
    BC = config_data['BC']
    E = config_data['E']
    out_dim = config_data['out_dim']
    u = u.reshape(batchsize, nx, out_dim)
    q = a[:, 2:-2, 1]

    uxx = (u[:, 1:-3, 0] - 2 * u[:, 2:-2, 0] + u[:, 3:-1, 0]) / dx ** 2
    uxxx = (-0.5 * u[:, :-4, 0] + u[:, 1:-3, 0] - u[:, 3:-1, 0] + 0.5 * u[:, 4:, 0]) / dx ** 3
    uxxxx = (u[:, :-4, 0] - 4 * u[:, 1:-3, 0] + 6 * u[:, 2:-2, 0] - 4 * u[:, 3:-1, 0] + u[:, 4:, 0]) / dx ** 4

    I = a[:, 2:-2, 0]
    Ix = (-0.5 * a[:, 1:-3, 0] + 0.5 * a[:, 3:-1, 0]) / dx
    Ixx = (a[:, 1:-3, 0] - 2 * a[:, 2:-2, 0] + a[:, 3:-1, 0]) / dx ** 2

    Du1 = I * uxxxx + 2 * Ix * uxxx + Ixx * uxx + q / E

    w = u[:, :, 0]
    w0 = torch.repeat_interleave(w[:, 0], nx, dim=0).reshape((batchsize, nx))
    wL = torch.repeat_interleave(w[:, -1], nx, dim=0).reshape((batchsize, nx))

    dw0 = (-1.5 * w[:, 0] + 2 * w[:, 1] - 0.5 * w[:, 2]) / dx
    dwdx0 = torch.repeat_interleave(dw0, nx, dim=0).reshape((batchsize, nx))

    dwL = (0.5 * w[:, -3] - 2 * w[:, -2] + 1.5 * w[:, -1]) / dx
    dwdxL = torch.repeat_interleave(dwL, nx, dim=0).reshape((batchsize, nx))

    d2w0 = (2 * w[:, 0] - 5 * w[:, 1] + 4 * w[:, 2] - w[:, 3]) / dx ** 2
    d2wdx0 = torch.repeat_interleave(d2w0, nx, dim=0).reshape((batchsize, nx))

    d2wL = (-w[:, -4] + 4 * w[:, -3] - 5 * w[:, -2] + 2 * w[:, -1]) / dx ** 2
    d2wdxL = torch.repeat_interleave(d2wL, nx, dim=0).reshape((batchsize, nx))

    d3wL = (1.5 * w[:, -5] - 7 * w[:, -4] + 12 * w[:, -3] - 9 * w[:, -2] + 2.5 * w[:, -1]) / dx ** 3
    d3wdxL = torch.repeat_interleave(d3wL, nx, dim=0).reshape((batchsize, nx))

    if BC == 'CF':
        boundary_l = torch.stack((w0, dwdx0), 1)
        boundary_r = torch.stack((d2wdxL, d3wdxL), 1)
    if BC == 'CS':
        boundary_l = torch.stack((w0, dwdx0), 1)
        boundary_r = torch.stack((wL, d2wdxL), 1)
    if BC == 'CC':
        boundary_l = torch.stack((w0, dwdx0), 1)
        boundary_r = torch.stack((wL, dwdxL), 1)
    if BC == 'SS':
        boundary_l = torch.stack((w0, d2wdx0), 1)
        boundary_r = torch.stack((wL, d2wdxL), 1)

    Du2 = torch.zeros_like(Du1)

    return Du1, Du2, boundary_l, boundary_r

@pino_loss_reduced_order2_1d
def FDM_ReducedOrder2_Euler_Bernoulli_Beam(config_data, a, u):
    batchsize = u.size(0)
    nx = u.size(1)
    dx = 1 / (nx - 1)
    BC = config_data['BC']
    E = config_data['E']
    out_dim = config_data['out_dim']
    u = u.reshape(batchsize, nx, out_dim)

    uxx = (u[:, :-2, 0] - 2 * u[:, 1:-1, 0] + u[:, 2:, 0]) / dx ** 2
    mxx = (u[:, :-2, 1] - 2 * u[:, 1:-1, 1] + u[:, 2:, 1]) / dx ** 2

    Du1 = mxx + a[:, 1:-1, 1]
    Du2 = uxx - u[:, 1:-1, 1] / (E * a[:, 1:-1, 0])

    w, m = u[:, :, 0], u[:, :, 1]
    w0 = torch.repeat_interleave(w[:, 0], nx, dim=0).reshape((batchsize, nx))
    wL = torch.repeat_interleave(w[:, -1], nx, dim=0).reshape((batchsize, nx))

    m0 = torch.repeat_interleave(m[:, 0], nx, dim=0).reshape((batchsize, nx))
    mL = torch.repeat_interleave(m[:, -1], nx, dim=0).reshape((batchsize, nx))

    dw0 = (-1.5 * w[:, 0] + 2 * w[:, 1] - 0.5 * w[:, 2]) / dx
    dwdx0 = torch.repeat_interleave(dw0, nx, dim=0).reshape((batchsize, nx))

    dwL = (0.5 * w[:, -3] - 2 * w[:, -2] + 1.5 * w[:, -1]) / dx
    dwdxL = torch.repeat_interleave(dwL, nx, dim=0).reshape((batchsize, nx))

    dmL = (0.5 * m[:, -3] - 2 * m[:, -2] + 1.5 * m[:, -1]) / dx
    dmdxL = torch.repeat_interleave(dmL, nx, dim=0).reshape((batchsize, nx))

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

    return Du1, Du2, boundary_l, boundary_r

@pino_loss_reduced_order2_1d
def FDM_ReducedOrder2_Euler_Bernoulli_Beam_BSF(config_data, a, u, bc):
    batchsize = u.size(0)
    nx = u.size(1)
    dx = 1 / (nx - 1)
    E = config_data['E']
    BC = config_data['BC']
    out_dim = config_data['out_dim']
    u = u.reshape(batchsize, nx, out_dim)

    mxx = (u[:, :-2, 1] - 2 * u[:, 1:-1, 1] + u[:, 2:, 1]) / dx ** 2
    uxx = (u[:, :-2, 0] - 2 * u[:, 1:-1, 0] + u[:, 2:, 0]) / dx ** 2
    I = a[:, 1:-1, 0]
    q = a[:, 1:-1, 1]
    m = u[:, 1:-1, 1]

    G1, G2, d2G1dx2, d2G2dx2, boundary_l, boundary_r \
        = boundary_function(u[:, :, 0], u[:, :, 1], bc, nx, dx, BC)

    Du1 = mxx - d2G2dx2[:, 1:-1] + q
    # Du2 = E * I * (uxx - d2G1dx2[:, 1:-1]) + m - G2[:, 1:-1]
    Du2 = (uxx - d2G1dx2[:, 1:-1]) - (m - G2[:, 1:-1]) / (E * I)

    return Du1, Du2, boundary_l, boundary_r