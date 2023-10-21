import torch
from torch import Tensor
from torchvision.models import resnet50

import torch.nn as nn
import torch.nn.functional as F

from typing import Any, List, Tuple

from ctran import ctranspath

class SaveFeatures:
    def __init__(self, module):
        self.hook = module.register_forward_hook(self.hook_fn)
        self.features = None
        
    def hook_fn(self, module, input, output):
        self.features = output
        
    def close(self):
        self.hook.remove()

class BasicBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride):
        super(BasicBlock, self).__init__()
        
        self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.relu = nn.ReLU(inplace=False)
        
        self.conv2 = nn.Conv2d(out_channels, out_channels, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_channels)
        
        self.downsample = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
            nn.BatchNorm2d(out_channels)
        )
        
    def forward(self, x):
        identity = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)
        
        identity = self.downsample(x)

        out += identity
        out = self.relu(out)

        return out

class AdaptationHead(nn.Module):
    def __init__(self, input_channels=3):
        super(AdaptationHead, self).__init__()
        
        self.initial = nn.Sequential(
            nn.Conv2d(input_channels, 64, kernel_size=7, stride=2, padding=3, bias=False),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=False),
            nn.MaxPool2d(kernel_size=3, stride=2, padding=1)
        )
        
        self.block1 = BasicBlock(64, 128, stride=2)
        self.block2 = BasicBlock(128, 256, stride=2)

        # A 1x1 conv to reduce channel dimensions from 256 to 3
        self.channel_reducer = nn.Conv2d(256, 3, kernel_size=1, stride=1, bias=False)
        
        self.final_pool = nn.AdaptiveAvgPool2d((224, 224))

    def forward(self, x):
        x = self.initial(x)
        x = self.block1(x)
        x = self.block2(x)

        x = self.channel_reducer(x)
        x = self.final_pool(x)

        return x

class ClassificationHead(nn.Module):
    def __init__(self, n_features, nhid, nlayers, n_output, dropout=0.5):
        super(ClassificationHead, self).__init__()

        # Define hidden layers
        self.hidden_layers = []
        for _ in range(nlayers - 2):
            self.hidden_layers.append(nn.Linear(nhid, nhid))
            self.hidden_layers.append(nn.ReLU())
            self.hidden_layers.append(nn.Dropout(dropout))

        # Define model layers
        self.layers = nn.Sequential(
            nn.Linear(n_features, nhid),
            nn.ReLU(),
            nn.Dropout(dropout),
            *self.hidden_layers,
            nn.Linear(nhid, n_output)
        )

    def forward(self, x):
        prediction = self.layers(x)
        return prediction

class ClassiModel(nn.Module):
    def __init__(self,
                pretrain_path: str = "",
                encode_features:int = 768,
                prediction_layers: int = 3,
                num_hidden_features: int = 96,
                num_classes: int = 3,
                drop_rate: float = 0.5):

        super(ClassiModel, self).__init__()
    
        self.adaptor = AdaptationHead(input_channels=3)

        ## Grad-CAM target layer
        self.adaptor_features = SaveFeatures(self.adaptor.block2)

        self.extractor = ctranspath()
        self.extractor.head = nn.Identity()
        td = torch.load(pretrain_path)
        self.extractor.load_state_dict(td['model'], strict=True)

        self.predictor = ClassificationHead(n_features=encode_features,
        nhid=num_hidden_features,
        nlayers=prediction_layers,
        n_output=num_classes,
        dropout=drop_rate)
    
    def forward(self, x: Tensor) -> Tensor:
        x = self.adaptor(x)
        features = self.extractor(x)
        out = self.predictor(features)
        return out

    def get_gradcam(self, x: Tensor, target_class: int):
        # Forward
        self.eval()
        out = self.forward(x)
        self.adaptor_features.features.register_hook(lambda grad: self.store_grad(grad.clone()))

        # Zero grads
        self.zero_grad()
        
        # Backward pass with respect to the target class
        one_hot_output = torch.FloatTensor(1, out.size()[-1]).zero_()
        one_hot_output[0][target_class] = 1
        out.backward(gradient=one_hot_output.to(x.device), retain_graph=True)
        
        # Get feature gradients and features
        grads_val = self.grads_val  # gradients
        target = self.adaptor_features.features  # features
        
        weights = torch.mean(grads_val, axis=(2, 3), keepdim=True)
        grad_cam = torch.sum((weights*target), axis=1, keepdim=True)
        grad_cam = F.relu(grad_cam)
        
        # Up-sample to input size
        grad_cam = F.interpolate(grad_cam, (x.size(2), x.size(3)), mode='bilinear', align_corners=False)
        
        return grad_cam

    def store_grad(self, grad):
        self.grads_val = grad