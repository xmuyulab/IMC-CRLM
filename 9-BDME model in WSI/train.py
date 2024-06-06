import os
import numpy as np
import pandas as pd
import math
import pickle
import random
import datetime
import argparse
import torch
from torch.utils.data import DataLoader
import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler
from torchvision import transforms

from F_dataset import *
from F_utils import *
from F_model import *
from F_train import *

def main(args):
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)
    torch.cuda.manual_seed_all(args.seed)

    # Load data
    train_images_path, train_images_label, val_images_path, val_images_label = read_split_data(args.all_train_images_path, num_categories=args.num_classes)

    with open(args.stats_file, 'rb') as f:
        transform_status = pickle.load(f)
        
    for key in transform_status.keys():
        transform_status[key]["mean"] = np.array([0.485, 0.456, 0.406])
        transform_status[key]["std"] = np.array([0.229, 0.224, 0.225])
        
    trnsfrms_val = transforms.Compose(
            [
                transforms.Resize(args.resize), 
                transforms.CenterCrop(args.center_crop),
                transforms.ToTensor()
            ]
        )
        
    # Create Dataloader
    train_dataset = MyDataSet_Train(images_path=train_images_path, images_class=train_images_label, transform=trnsfrms_val, stats_file=transform_status)
    val_dataset = MyDataSet_Train(images_path=val_images_path, images_class=val_images_label, transform=trnsfrms_val, stats_file=transform_status)

    nw = min([os.cpu_count(), args.batch_size if args.batch_size > 1 else 0, args.num_workers])  
    print('Using {} dataloader workers every process'.format(nw))
    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, pin_memory=True, num_workers=nw, collate_fn=train_dataset.collate_fn)
    val_loader = DataLoader(val_dataset, batch_size=args.batch_size, shuffle=False, pin_memory=True, num_workers=nw, collate_fn=val_dataset.collate_fn)

    # Initiate model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = ClassiModel(pretrain_path=args.pretrain_path, encode_features=args.encode_features, prediction_layers=args.prediction_layers, num_hidden_features=args.num_hidden_features, num_classes=args.num_classes, drop_rate=args.drop_rate)

    for name, para in model.named_parameters():
        if "extractor" in name:
            para.requires_grad_(False)

    model = nn.DataParallel(model).to(device)
    pg = [p for p in model.parameters() if p.requires_grad]

    # Load Pretrain model
    if args.LoadLastModel:
        model.load_state_dict(torch.load("/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/weights/2023-09-07_14-06-08/model-150.pth"))

    optimizer = optim.Adam(pg, lr=args.learning_rate, weight_decay=args.weight_decay)

    lrf = args.lrf
    lf = lambda x: ((1 + math.cos(x * math.pi / args.epochs)) / 2) * (1 - lrf) + lrf  
    scheduler = lr_scheduler.LambdaLR(optimizer, lr_lambda=lf)

    ## saving model weights according to timepoint
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    weights_dir = os.path.join(args.weights_dir,current_time)
    os.makedirs(weights_dir, exist_ok=True)  # Make sure the directory exists


    val_info = {}  # Initialize an empty dictionary to store validation metrics
    
    for epoch in range(1, args.epochs + 1):
        # Train
        mean_loss, train_acc = train_one_epoch(model=model, optimizer=optimizer, data_loader=train_loader, device=device, epoch=epoch)

        scheduler.step()
        print("[epoch {}] train accuracy: {}, mean loss: {}".format(epoch, round(train_acc, 3), mean_loss))


        # Validate
        if epoch % args.val_freq == 0:
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

        if epoch == args.epochs:
            # Save the dictionary to val.pkl after the last epoch
            with open(os.path.join(weights_dir, "val.pkl"), 'wb') as f:
                pickle.dump(val_info, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments for training and inference.')

    parser.add_argument('--seed', type=int, default=619, help='Seed value for initializing random generators.')
    parser.add_argument('--all_train_images_path', default="/mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data/", help='Path to the directory containing all training images.')
    
    parser.add_argument('--stats_file', default="stats.pkl", help='Path to the file containing statistics for data transformation.')
    parser.add_argument('--resize', type=int, default=1024, help='Size to which the images are resized.')
    parser.add_argument('--center_crop', type=int, default=1024, help='Size of the central crop taken from the image after resizing.')

    parser.add_argument('--batch_size', type=int, default=32, help='Number of samples per gradient update.')
    parser.add_argument('--num_workers', type=int, default=8, help='Number of subprocesses to use for data loading.')
    parser.add_argument('--pretrain_path', default="./pretrainModel/ctranspath.pth", help='Path to the pre-trained model.')

    parser.add_argument('--encode_features', type=int, default=768, help='Number of features for encoding layers.')
    parser.add_argument('--prediction_layers',  type=int, default=3, help='Number of features for prediction layers.')
    parser.add_argument('--num_hidden_features', type=int, default=64, help='Number of hidden features in the model.')
    parser.add_argument('--num_classes', type=int, default=4, help='Number of classes for classification.')

    parser.add_argument('--drop_rate', type=float, default=0.5, help='Dropout rate for regularization.')
    parser.add_argument('--learning_rate', type=float, default=0.001, help='Learning rate for optimizer.')
    parser.add_argument('--weight_decay', type=float, default=1e-4, help='Weight decay (L2 penalty) for the optimizer.')
    parser.add_argument('--lrf', type=float, default=0.1, help='Learning rate factor for the scheduler.')
    parser.add_argument('--epochs', type=int, default=30, help='Number of times the model will iterate over the entire training dataset.')

    parser.add_argument('--LoadLastModel', action='store_true', help='Whether to load the trained model.')

    parser.add_argument('--val_freq', type=int, default=5, help='Frequency for performing validation.')
    parser.add_argument('--weights_dir', default="./weights/", help='Directory where model weights are saved.')

    args = parser.parse_args()
    
    main(args)

# nohup \
# python train.py \
# --seed 619 \
# --all_train_images_path /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data \
# --stats_file /mnt/data/lyx/IMC/WSI_HE/HE_classification_2/stats.pkl \
# --resize 1024 \
# --center_crop 1024 \
# --batch_size 32 \
# --num_hidden_features 96 \
# --num_classes 4 \
# --epochs 200 \
# --learning_rate 0.001 \
# --val_freq 20
# > train.log 2>&1 &

# nohup \
# python train.py --seed 619 --all_train_images_path /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data \
# --stats_file /mnt/data/lyx/IMC/WSI_HE/HE_classification_2/stats.pkl \
# --resize 1024 --center_crop 1024 --batch_size 32 \
# --num_hidden_features 96 --num_classes 4 --epochs 50 --LoadLastModel --val_freq 5 \
# > train.log 2>&1 &