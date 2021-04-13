import torch
import pandas as pd
import numpy as np
import glob
import os
# import random
import torch.nn.functional as F

# import matplotlib.pyplot as plt
# %matplotlib inline
import json

from torch import nn
from scipy import ndimage
from torchvision import transforms

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
dataset_dir = 'yeast_processed_log_norm_area_50'
model_dir = 'unet2D_aug_norm'
normalize = True


## TODO: refactor into a transforms file
## Define Transforms
if normalize:
    with open(os.path.join('../data', dataset_dir, 'stds.json'), 'r') as f:
        stds = json.load(f)
        
    normalize_transform = transforms.Normalize(
                                        stds[0],
                                        stds[1])
else:
    normalize_transform=None

class RandomTranspose(object):
    """Randomly transpose sample
    """

    def __call__(self, sample):
        transpose = random.random() < 0.5
        if transpose:
            # sample is CxHxW
            return torch.transpose(sample, -1, -2)
        else:
            return sample


transform_dict = {'train': transforms.Compose([
                transforms.Resize((512,512), interpolation=0),
                RandomTranspose(),
                transforms.RandomHorizontalFlip(),
                transforms.RandomVerticalFlip(),
            ]),
             'val': transforms.Compose([
                transforms.Resize((512,512), interpolation=0),
            ])}

## Define model
# TODO: Refactor into own model file
# Model
# Inspired by https://github.com/MrGiovanni/ModelsGenesis
class ContBatchNorm2d(nn.modules.batchnorm._BatchNorm):
    def _check_input_dim(self, input):

        if input.dim() != 4:
            raise ValueError('expected 4D input (got {}D input)'.format(input.dim()))
        #super(ContBatchNorm3d, self)._check_input_dim(input)

    def forward(self, input):
        self._check_input_dim(input)
        return F.batch_norm(
            input, self.running_mean, self.running_var, self.weight, self.bias,
            True, self.momentum, self.eps)


class LUConv(nn.Module):
    def __init__(self, in_chan, out_chan, act):
        super(LUConv, self).__init__()
        self.conv1 = nn.Conv2d(in_chan, out_chan, kernel_size=3, padding=1)
        self.bn1 = ContBatchNorm2d(out_chan)

        if act == 'relu':
            self.activation = nn.ReLU(out_chan)
        elif act == 'prelu':
            self.activation = nn.PReLU(out_chan)
        elif act == 'elu':
            self.activation = nn.ELU(inplace=True)
        else:
            raise

    def forward(self, x):
        out = self.activation(self.bn1(self.conv1(x)))
        return out


def _make_nConv(in_channel, depth, act, double_chnnel=False):
    if double_chnnel:
        layer1 = LUConv(in_channel, 32 * (2 ** (depth+1)),act)
        layer2 = LUConv(32 * (2 ** (depth+1)), 32 * (2 ** (depth+1)),act)
    else:
        layer1 = LUConv(in_channel, 32*(2**depth),act)
        layer2 = LUConv(32*(2**depth), 32*(2**depth)*2,act)

    return nn.Sequential(layer1,layer2)

class DownTransition(nn.Module):
    def __init__(self, in_channel,depth, act):
        super(DownTransition, self).__init__()
        self.ops = _make_nConv(in_channel, depth,act)
        self.maxpool = nn.MaxPool2d(2)
        self.current_depth = depth

    def forward(self, x):
        if self.current_depth == 3:
            out = self.ops(x)
            out_before_pool = out
        else:
            out_before_pool = self.ops(x)
            out = self.maxpool(out_before_pool)
        return out, out_before_pool

class UpTransition(nn.Module):
    def __init__(self, inChans, outChans, depth,act):
        super(UpTransition, self).__init__()
        self.depth = depth
        self.up_conv = nn.ConvTranspose2d(inChans, outChans, kernel_size=2, stride=2)
        self.ops = _make_nConv(inChans+ outChans//2,depth, act, double_chnnel=True)

    def forward(self, x, skip_x):
        out_up_conv = self.up_conv(x)
        concat = torch.cat((out_up_conv,skip_x),1)
        out = self.ops(concat)
        return out


class OutputTransition(nn.Module):
    def __init__(self, inChans, n_labels):

        super(OutputTransition, self).__init__()
        self.final_conv = nn.Conv2d(inChans, n_labels, kernel_size=1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        out = self.sigmoid(self.final_conv(x))
        return out

class UNet2D(nn.Module):
    # the number of convolutions in each layer corresponds
    # to what is in the actual prototxt, not the intent
    def __init__(self, n_class=1, n_channels=1, act='relu'):
        super(UNet2D, self).__init__()

        self.down_tr64 = DownTransition(n_channels,0,act)
        self.down_tr128 = DownTransition(64,1,act)
        self.down_tr256 = DownTransition(128,2,act)
        self.down_tr512 = DownTransition(256,3,act)

        self.up_tr256 = UpTransition(512, 512,2,act)
        self.up_tr128 = UpTransition(256,256, 1,act)
        self.up_tr64 = UpTransition(128,128,0,act)
        self.out_tr = OutputTransition(64, n_class)

    def forward(self, x):
        self.out64, self.skip_out64 = self.down_tr64(x)
        self.out128,self.skip_out128 = self.down_tr128(self.out64)
        self.out256,self.skip_out256 = self.down_tr256(self.out128)
        self.out512,self.skip_out512 = self.down_tr512(self.out256)

        self.out_up_256 = self.up_tr256(self.out512,self.skip_out256)
        self.out_up_128 = self.up_tr128(self.out_up_256, self.skip_out128)
        self.out_up_64 = self.up_tr64(self.out_up_128, self.skip_out64)
        self.out = self.out_tr(self.out_up_64)

        return self.out



print(f"Loading model at ../models/{model_dir}/{dataset_dir}/best_model.pth")
checkpoint = torch.load(f'../models/{model_dir}/{dataset_dir}/best_model.pth',
                        map_location=torch.device('cpu'))

model = checkpoint['model']
model.cuda()
model.eval()
torch.set_grad_enabled(False)

print("model loaded!")


def mask_to_box_many(a, max_label=0):
    labels, n = ndimage.measurements.label(a, np.ones((3, 3)))
              
    if max_label < n:
                    
        # find largest contiguous areas
        # print(n)
        uniques, counts = np.unique(labels, return_counts=True)
        # print(max_label)
        # print(counts)
        ind = np.argpartition(counts, -(max_label+1))[-(max_label+1):] # find indices of top max_label counts; +1 because 0 will be one of them
        # print(ind)
        
        top_uniques = uniques[ind]
        top_counts = counts[ind]
        
        # remove 0
        if 0 in top_uniques:
            rm_i = np.where(top_uniques == 0)
        else:
            # we have 1 too many, remove smallest value
            rm_i = np.argmin(top_counts)
            
        top_counts = np.delete(top_counts, rm_i)
        top_uniques = np.delete(top_uniques, rm_i)
        
#         print(top_uniques)
#         print(top_counts)
#         print(type(top_uniques))
#         print(np.vectorize(lambda label: label in top_uniques)(labels))
        
        # zero out anything that's not in top max_label uniques:
        labels = np.where(np.vectorize(lambda label: label in top_uniques)(labels), labels, 0)
        
#         print(np.unique(labels))
    
        
        
    
    objs = ndimage.find_objects(labels)

    # filter out nones
    return list(filter(None, objs)) 

def mask_to_box(a, max_label=1):
    objs = ndimage.find_objects(a, max_label=max_label)

    # Get the height and width
    return objs

###############################################
## Convert to bounding box
class DatasetBB(torch.utils.data.Dataset):
    
    def __init__(self, folder_path, transform = None, normalize_transform = None, interp_size = 512):
        super(DatasetBB, self).__init__()
        if transform:
            self.transform = transform
        else:
            self.transform = transforms.compose([
                transforms.Resize((interp_size, interp_size), interpolation=0)
            ])
        if normalize_transform:
            self.normalize_transform = normalize_transform
        else:
            self.normalize_transform = None
        
        self.interp_size = interp_size
        self.img_files = glob.glob(os.path.join(folder_path,'features','*.npy'))
        self.mask_files = []
        self.pair_names = []
        for img_path in self.img_files:
            self.mask_files.append(os.path.join(folder_path,'masks',os.path.basename(img_path)) )
            self.pair_names.append(os.path.splitext(os.path.basename(img_path))[0])

    def __getitem__(self, index):
        img_path = self.img_files[index]
        mask_path = self.mask_files[index]

        data = np.load(img_path)
        label = np.load(mask_path)
        
        data = torch.from_numpy(data).float()
        label = torch.from_numpy(label).float()
        
        samp = torch.cat((data, label))
        
        samp = self.transform(samp)
        data = samp[:-1]
        label = samp[-1].unsqueeze(0)
        
        if self.normalize_transform:
            data = self.normalize_transform(data)
            

#         data = F.interpolate(data.unsqueeze(0), (self.interp_size, self.interp_size)).squeeze(0)
#         label = F.interpolate(label.unsqueeze(0), (self.interp_size, self.interp_size)).squeeze(0)
        return data, label, self.pair_names[index]

    def __len__(self):
        return len(self.img_files)

print(f'loading data at ../data/{dataset_dir}/test')
ds_test = DatasetBB(os.path.join('../data', dataset_dir, 'test'), transform=transform_dict['val'], normalize_transform=normalize_transform)
# ds_val = DatasetBB(os.path.join('data', dataset_dir, 'val'), transform=transform_dict['val'], normalize_transform=normalize_transform)

test_dl = torch.utils.data.DataLoader(ds_test, batch_size=1, num_workers=8, shuffle=False)
# val_dl = torch.utils.data.DataLoader(ds_val, batch_size=1, num_workers=8, shuffle=False)

top_n = 3


def get_predictions_to_bbs(dataloader, model, top_n):
    
    model.cuda()
    model.eval()
    torch.set_grad_enabled(False)
    
    bb_preds = []
    bb_preds_many = []
    bb_trues = []
    pairnames = []
    
    for _, (data, mask, pairname) in enumerate(dataloader):
        data, mask = data.to(device), mask.to(device)

        pred = model(data)

        pred = pred[0,0,:,:].cpu().numpy().round().astype(int)
        mask = mask[0,0,:,:].cpu().numpy().round().astype(int)
        

        bb_preds.append(mask_to_box(pred))
        bb_preds_many.append(mask_to_box_many(pred, max_label=top_n))
        bb_trues.append(mask_to_box(mask))
        pairnames.append(pairname)
        
#         plt.matshow(pred)
        
#         break
    
    return pd.DataFrame(list(zip(pairnames, bb_preds, bb_preds_many, bb_trues)), columns=['pair', 'pred', 'pred_many', 'true'])

print('Turn predictions into bounding boxes')
test_df = get_predictions_to_bbs(test_dl, model, top_n)
# val_df = get_predictions_to_bbs(val_dl, model, top_n)

# -1 to make inclusive
deslice = lambda y: [((x[0].start, x[0].stop-1), (x[1].start, x[1].stop-1)) for x in y]

test_df[['pred', 'true', 'pred_many']] = test_df[['pred', 'true', 'pred_many']].applymap(deslice)
# val_df[['pred', 'true', 'pred_many']] = val_df[['pred', 'true', 'pred_many']].applymap(deslice)

test_df['pair'] = test_df['pair'].apply(lambda x: x[0]).str.split(pat='_')

test_df[['protein_a', 'protein_b']] = pd.DataFrame(test_df['pair'].tolist(), index=test_df.index)
test_df.drop(['pair'], inplace=True, axis=1)

print('Finding dice/iou of big bounding box')
# Formula adapted from https://stackoverflow.com/questions/25349178/calculating-percentage-of-bounding-box-overlap-for-image-detector-evaluation
def get_dice_iou(bb1, bb2):
    """
    Calculate the Intersection over Union (IoU) of two bounding boxes.

    Parameters
    ----------
    bb1 : list: [[y1, y2], [x1, x2]]
        The (x1, y1) position is at the top left corner,
        the (x2, y2) position is at the bottom right corner
    bb1 : list: [[y1, y2], [x1, x2]]
        The (x1, y1) position is at the top left corner,
        the (x2, y2) position is at the bottom right corner

    Returns
    -------
    float
        in [0, 1]
    """

    assert bb1[0][0] <= bb1[0][1]
    assert bb1[1][0] <= bb1[1][1]
    assert bb2[0][0] <= bb2[0][1]
    assert bb2[1][0] <= bb2[1][1]

    # determine the coordinates of the intersection rectangle
    x_left = max(bb1[1][0], bb2[1][0])
    y_top = max(bb1[0][0], bb2[0][0])
    x_right = min(bb1[1][1], bb2[1][1])
    y_bottom = min(bb1[0][1], bb2[0][1])

    if x_right < x_left or y_bottom < y_top:
        return 0.0, 0.0

    # The intersection of two axis-aligned bounding boxes is always an
    # axis-aligned bounding box
    intersection_area = (x_right - x_left + 1) * (y_bottom - y_top + 1)

    # compute the area of both AABBs
    bb1_area = (bb1[1][1] - bb1[1][0] + 1) * (bb1[0][1] - bb1[0][0] + 1)
    bb2_area = (bb2[1][1] - bb2[1][0] + 1) * (bb2[0][1] - bb2[0][0] + 1)

    # compute the intersection over union by taking the intersection
    # area and dividing it by the sum of prediction + ground-truth
    # areas - the interesection area
    iou = intersection_area / float(bb1_area + bb2_area - intersection_area)
    
    # dice
    dice =  2 * intersection_area / float(bb1_area + bb2_area)
    assert iou >= 0.0
    assert iou <= 1.0
    return iou, dice

test_df['iou_dice'] = test_df.apply(lambda x: get_dice_iou(x.pred[0], x.true[0]), axis=1)
test_df[['iou', 'dice']] = pd.DataFrame(test_df['iou_dice'].tolist(), index=test_df.index)
test_df.drop(['iou_dice'], inplace=True, axis=1)

# val_df['iou_dice'] = val_df.apply(lambda x: get_dice_iou(x.pred[0], x.true[0]), axis=1)
# val_df[['iou', 'dice']] = pd.DataFrame(val_df['iou_dice'].tolist(), index=val_df.index)
# val_df.drop(['iou_dice'], inplace=True, axis=1)

## Distance Measure

# The distance measure is as described in PIPE-Sites: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-225

print('finding dm of big bounding box')
def dm(pred, true, len_a, len_b):
    
    assert pred[0][0] <= pred[0][1]
    assert pred[1][0] <= pred[1][1]
    assert true[0][0] <= true[0][1]
    assert true[1][0] <= true[1][1]
    
    delta_a = max(true[0][0]-pred[0][0], pred[0][1]-true[0][1], 0)/len_a
    delta_b = max(true[1][0]-pred[1][0], pred[1][1]-true[1][1], 0)/len_b
    
    result = ((delta_a**2 + delta_b**2)**0.5)/(2**0.5)
    
    return result

test_df['dm'] = test_df.apply(lambda x: dm(x.pred[0], x.true[0], 512, 512), axis=1)

## Best out of 3 distance measure
print('Finding best out of 3 largest contiguous regions distance measure')
def best_dm_idx(preds, true, len_a, len_b):
    dms = np.asarray([dm(pred, true, len_a, len_b) for pred in preds])
    
    if dms.min() > 0.:
        return dms.argmin()
    else:
        idxs = np.argwhere(dms > 0.)
        areas = np.asarray([(pred[0][1]-pred[0][0])*(pred[1][1]-pred[1][0]) for pred in preds])
        # no perfect dm, zero out area
        areas[idxs] = 0
        
        return areas.argmax()


test_df['best_dm_idx'] = test_df.apply(lambda x: best_dm_idx(x.pred_many, x.true[0], 512, 512), axis=1)

test_df['best_dm'] = test_df.apply(lambda x: dm(x.pred_many[x.best_dm_idx], x.true[0], 512, 512), axis=1)

test_df['best_iou_dice'] = test_df.apply(lambda x: get_dice_iou(x.pred_many[x.best_dm_idx], x.true[0]), axis=1)
test_df[['best_iou', 'best_dice']] = pd.DataFrame(test_df['best_iou_dice'].tolist(), index=test_df.index)
test_df.drop(['best_iou_dice'], inplace=True, axis=1)


print(f"Whole region bounding box scores: \n \
    iou: {test_df['iou'].mean()}, dice: {test_df['dice'].mean()}, distance metric: {test_df['dm'].mean()}")
print(f"Best of 3 contiguous regions: \n \
    iou: {test_df['best_iou'].mean()}, dice: {test_df['best_dice'].mean()}, distance metric: {test_df['best_dm'].mean()}")

