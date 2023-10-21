import os
import numpy as np
import pandas as pd
import math
import pickle
import random
import datetime
import torch
from torch.utils.data import DataLoader
import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler
from torchvision import transforms

from F_dataset import *
from F_utils import *
from F_model import *
from F_train import *

seed=619
all_train_images_path="/mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data"
num_classes=4
stats_file="/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/stats.pkl"
resize=1024
center_crop=1024
batch_size=32
num_workers=8
pretrain_path="./pretrainModel/ctranspath.pth"
encode_features=768
prediction_layers=3
num_hidden_features=96
drop_rate=0.5
learning_rate=0.001
weight_decay=1e-4
lrf=0.1
epochs=200
weights_dir="./weights/"

torch.manual_seed(seed)
np.random.seed(seed)
random.seed(seed)
torch.cuda.manual_seed_all(seed)

# Load data
train_images_path, train_images_label, val_images_path, val_images_label = read_split_data(all_train_images_path, num_categories=num_classes)

with open(stats_file, 'rb') as f:
    transform_status = pickle.load(f)
    
for key in transform_status.keys():
    transform_status[key]["mean"] = np.array([0.485, 0.456, 0.406])
    transform_status[key]["std"] = np.array([0.229, 0.224, 0.225])
    
trnsfrms_val = transforms.Compose(
        [
            transforms.Resize(resize), 
            transforms.CenterCrop(center_crop),
            transforms.ToTensor()
        ]
    )
    
# Create Dataloader
train_dataset = MyDataSet_Train(images_path=train_images_path, images_class=train_images_label, transform=trnsfrms_val, stats_file=transform_status)
val_dataset = MyDataSet_Train(images_path=val_images_path, images_class=val_images_label, transform=trnsfrms_val, stats_file=transform_status)

nw = min([os.cpu_count(), batch_size if batch_size > 1 else 0, num_workers])  
print('Using {} dataloader workers every process'.format(nw))
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, pin_memory=True, num_workers=nw, collate_fn=train_dataset.collate_fn)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, pin_memory=True, num_workers=nw, collate_fn=val_dataset.collate_fn)

# Initiate model
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = ClassiModel(pretrain_path=pretrain_path, encode_features=encode_features, prediction_layers=prediction_layers, num_hidden_features=num_hidden_features, num_classes=num_classes, drop_rate=drop_rate)

for name, para in model.named_parameters():
    if "extractor" in name:
        para.requires_grad_(False)

model = nn.DataParallel(model).to(device)
pg = [p for p in model.parameters() if p.requires_grad]

# Load Pretrain model
if LoadLastModel:
    model.load_state_dict(torch.load("/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/weights/2023-09-07_14-06-08/model-150.pth"))

optimizer = optim.Adam(pg, lr=learning_rate, weight_decay=weight_decay)

lrf = lrf
lf = lambda x: ((1 + math.cos(x * math.pi / epochs)) / 2) * (1 - lrf) + lrf  
scheduler = lr_scheduler.LambdaLR(optimizer, lr_lambda=lf)

## saving model weights according to timepoint
current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
weights_dir = os.path.join(weights_dir,current_time)
os.makedirs(weights_dir, exist_ok=True)  # Make sure the directory exists


val_info = {}  # Initialize an empty dictionary to store validation metrics

for epoch in range(1, epochs + 1):
    # Train
    mean_loss, train_acc = train_one_epoch(model=model, optimizer=optimizer, data_loader=train_loader, device=device, epoch=epoch)

    scheduler.step()
    print("[epoch {}] train accuracy: {}, mean loss: {}".format(epoch, round(train_acc, 3), mean_loss))


    # Validate
    if epoch % val_freq == 0:
        val_accuracy, val_pre, val_rec, val_f1 = evaluate(model=model, data_loader=val_loader, device=device)

        print("[epoch {}] val accuracy: {}, precision: {}, recall: {}, f1-score: {}".format(epoch, round(val_accuracy, 3), round(val_pre, 3), round(val_rec, 3), round(val_f1, 3)))

        torch.save(model.state_dict(), os.path.join(weights_dir, f"model-{epoch}.pth"))

        # Save the validation metrics to the dictionary
        val_info[epoch] = {
            "accuracy": val_accuracy,
            "precision": val_pre,
            "recall": val_rec,
            "f1-score": val_f1
        }

    if epoch == epochs:
        # Save the dictionary to val.pkl after the last epoch
        with open(os.path.join(weights_dir, "val.pkl"), 'wb') as f:
            pickle.dump(val_info, f)


# nohup \
# python train.py --seed 619 --all_train_images_path /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data \
# --stats_file /mnt/data/lyx/IMC/WSI_HE/HE_classification_2/stats.pkl \
# --resize 1024 --center_crop 1024 --batch_size 32 \
# --num_hidden_features 96 --num_classes 4 --epochs 200 --val_freq 20 \
# > train.log 2>&1 &

# nohup \
# python train.py --seed 619 --all_train_images_path /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data \
# --stats_file /mnt/data/lyx/IMC/WSI_HE/HE_classification_2/stats.pkl \
# --resize 1024 --center_crop 1024 --batch_size 32 \
# --num_hidden_features 96 --num_classes 4 --epochs 50 --LoadLastModel --val_freq 5 \
# > train.log 2>&1 &