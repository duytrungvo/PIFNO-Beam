import scipy.io
import numpy as np

try:
    from pyDOE import lhs
    # Only needed for PINN's dataset
except ImportError:
    lhs = None

import torch

class MatReader(object):
    def __init__(self, file_path, to_torch=True, to_cuda=False, to_float=True):
        super(MatReader, self).__init__()

        self.to_torch = to_torch
        self.to_cuda = to_cuda
        self.to_float = to_float

        self.file_path = file_path

        self.data = None
        self.old_mat = None
        self._load_file()

    def _load_file(self):
        self.data = scipy.io.loadmat(self.file_path)
        self.old_mat = True

    def load_file(self, file_path):
        self.file_path = file_path
        self._load_file()

    def read_field(self, field):
        x = self.data[field]

        if not self.old_mat:
            x = x[()]
            x = np.transpose(x, axes=range(len(x.shape) - 1, -1, -1))

        if self.to_float:
            x = x.astype(np.float32)

        if self.to_torch:
            x = torch.from_numpy(x)

            if self.to_cuda:
                x = x.cuda()

        return x

    def set_cuda(self, to_cuda):
        self.to_cuda = to_cuda

    def set_torch(self, to_torch):
        self.to_torch = to_torch

    def set_float(self, to_float):
        self.to_float = to_float

class Loader_1D(object):
    def __init__(self, datapath, nx=2**10+1, sub=8, in_dim=1, out_dim=1):
        dataloader = MatReader(datapath)
        self.sub = sub
        self.s = int(np.ceil(nx / sub))
        self.in_dim = in_dim
        self.out_dim = out_dim
        if len((dataloader.read_field('input')).size()) == 2:
            self.x_data = dataloader.read_field('input')[:, ::sub].unsqueeze(2)
        else:
            self.x_data = dataloader.read_field('input')[:, ::sub, :in_dim-1]
        if len((dataloader.read_field('output')).size()) == 2:
            self.y_data = dataloader.read_field('output')[:, ::sub].unsqueeze(2)
        else:
            self.y_data = dataloader.read_field('output')[:, ::sub, :out_dim]
        self.gridx = dataloader.read_field('x')[:, ::sub]

    def make_loader(self, n_sample, batch_size, start=0, train=True):
        xs = self.x_data[start:start + n_sample]
        ys = self.y_data[start:start + n_sample]

        # xs = xs.reshape(n_sample, self.s)
        # xs = torch.stack([xs, self.gridx.repeat([n_sample, 1])], dim=2)
        xs = torch.cat((xs, self.gridx.repeat([n_sample, 1]).unsqueeze(2)), 2)
        dataset = torch.utils.data.TensorDataset(xs, ys)
        if train:
            loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)
        else:
            loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
        return loader