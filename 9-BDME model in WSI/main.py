import os
import numpy as np
import pandas as pd
import math
import pickle
import random
from collections import Counter

import torch
from torch.utils.data import DataLoader
import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler

from torchvision import transforms

from F_dataset import *
from F_utils import *
from F_model import *
from F_train import *

seed = 619
torch.manual_seed(seed)
np.random.seed(seed)
random.seed(seed)
torch.cuda.manual_seed_all(seed)

## load data
all_train_images_path = "/mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data/"

train_images_path, train_images_label, val_images_path, val_images_label = read_split_data(all_train_images_path, num_categories=3)

## Create Dataloader
# mean = (0.485, 0.456, 0.406)
# std = (0.229, 0.224, 0.225)
# trnsfrms_val = transforms.Compose(
#         [
#             transforms.Resize(1024), ## 224
#             transforms.CenterCrop(1024), ## 224
#             transforms.ToTensor(),
#             transforms.Normalize(mean = mean, std = std)
#         ]
#     )

# train_dataset = MyDataSet_Train(images_path = train_images_path,images_class = train_images_label,transform = trnsfrms_val)
# val_dataset = MyDataSet_Train(images_path = val_images_path,images_class = val_images_label,transform = trnsfrms_val)

f = open("stats.pkl", 'rb')
transform_status = pickle.load(f)
f.close()

trnsfrms_val = transforms.Compose(
        [
            transforms.Resize(1024), ## 224
            transforms.CenterCrop(1024), ## 224
            transforms.ToTensor()
        ]
    )

train_dataset = MyDataSet_Train(images_path = train_images_path,images_class = train_images_label, transform=trnsfrms_val, stats_file = transform_status)
val_dataset = MyDataSet_Train(images_path = val_images_path,images_class = val_images_label, transform=trnsfrms_val,stats_file = transform_status)

batch_size = 32
nw = min([os.cpu_count(), batch_size if batch_size > 1 else 0, 8])  # number of workers
print('Using {} dataloader workers every process'.format(nw))
train_loader = DataLoader(train_dataset,
                                            batch_size=batch_size,
                                            shuffle=True,
                                            pin_memory=True,
                                            num_workers=nw,
                                            collate_fn=train_dataset.collate_fn)

val_loader = DataLoader(val_dataset,
                                            batch_size=batch_size,
                                            shuffle=False,
                                            pin_memory=True,
                                            num_workers=nw,
                                            collate_fn=val_dataset.collate_fn)

## Initiate model
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = ClassiModel(
    pretrain_path='/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/pretrainModel/ctranspath.pth',
    encode_features=768,
    prediction_layers=3,
    num_hidden_features=96,
    num_classes=3,
    drop_rate=0.5
)

for name, para in model.named_parameters():
            ## freeze
            if "extractor" in name:
                para.requires_grad_(False)

# model = densenet121(num_classes=3)

model = nn.DataParallel(model).to(device)

# for name, param in model.named_parameters():
#     if param.requires_grad:
#         print(name)

pg = [p for p in model.parameters() if p.requires_grad]
# optimizer = optim.SGD(pg, lr=0.001, momentum=0.9, weight_decay=1E-4, nesterov=True)
optimizer = optim.Adam(pg, lr=0.001,weight_decay=1E-4)

lrf = 0.1
epochs = 30

lf = lambda x: ((1 + math.cos(x * math.pi / epochs)) / 2) * (1 - lrf) + lrf  # cosine
scheduler = lr_scheduler.LambdaLR(optimizer, lr_lambda=lf)

for epoch in range(1,epochs+1):
    # train
    mean_loss, train_acc = train_one_epoch(
        model=model,
        optimizer=optimizer,
        data_loader=train_loader,
        device=device,
        epoch=epoch)

    scheduler.step()
    print("[epoch {}] train accuracy: {}".format(epoch, round(train_acc, 3)))

    # validate
    if epoch % 5 == 0:
        val_accuracy, val_pre, val_rec, val_f1 = evaluate(model=model,
                        data_loader=val_loader,
                        device=device)
        
        print("[epoch {}] val accuracy: {}, precision: {}, recall: {}, f1-score: {}".format(epoch, round(val_accuracy, 3), round(val_pre, 3), round(val_rec, 3), round(val_f1, 3)))

        torch.save(model.state_dict(), "./weights/model-{}.pth".format(epoch))

## Grad-CAM
# do_GradCAM = False
do_GradCAM = True

if do_GradCAM:
    model = ClassiModel(
    pretrain_path='/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/pretrainModel/ctranspath.pth',
    encode_features=768,
    prediction_layers=3,
    num_hidden_features=96,
    num_classes=3,
    drop_rate=0.5
    )

    model = nn.DataParallel(model).to(device)
    model.load_state_dict(torch.load("/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/weights/model-30.pth"))

    saveBasePath = "/mnt/data/lyx/IMC/WSI_HE/WSI_Data/GradCAM"
    if not os.path.exists(saveBasePath):
        os.mkdir(saveBasePath)

    ori_img_path = "/mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data"
    # for dirtype in os.listdir(ori_img_path):
    for dirtype in ["val"]:
        catetypes = os.listdir(os.path.join(ori_img_path,dirtype))
        for catetype in catetypes:
            imgs = os.listdir(os.path.join(ori_img_path,dirtype,catetype))

            for img in imgs:
                if not os.path.exists(os.path.join(saveBasePath,dirtype,catetype,img)):
                    os.makedirs(os.path.join(saveBasePath,dirtype,catetype,img))

                for target_class in catetypes:
                    visualize_gradcam(
                            img_path=os.path.join(ori_img_path,dirtype,catetype,img) ,
                            transform=trnsfrms_val, mean_ = transform_status[catetype]["mean"], std_ = transform_status[catetype]["std"],
                            model=model,
                            target_class=int(target_class),device=device,
                            save_name=os.path.join(saveBasePath,dirtype,catetype,img,"pred_"+target_class+".png") 
                            )
    


## Infer on WSI
WSIPath = "/mnt/data/lyx/IMC/WSI_HE/WSI_Trainging/WSI_crop"
PIDs = os.listdir(WSIPath)

pred = {}

for PID in PIDs:
    infer_images_dir = os.path.join(WSIPath,PID)
    infer_images_path = [os.path.join(infer_images_dir,x) for x in os.listdir(infer_images_dir)]
    mean_, std_ = compute_mean_std(infer_images_dir)

    trnsfrms_val = transforms.Compose(
            [
                transforms.Resize(1024), ## 224
                transforms.CenterCrop(1024), ## 224
                transforms.ToTensor(),
                transforms.Normalize(mean = mean_, std = std_)
            ]
        )
    
    infer_dataset = MyDataSet_Infer(images_path = infer_images_path,transform = trnsfrms_val)
    infer_dataloader=DataLoader(infer_dataset,
                                batch_size=96,
                                shuffle=False,
                                pin_memory=True,
                                num_workers=nw)

    predTensor = infer(model = model, data_loader=infer_dataloader,device=device)
    # predTensor.shape
    pred_patch = torch.argmax(predTensor, dim=1)
    pred_patch_counts = Counter(pred_patch.tolist())
    predTensormean = np.array(torch.mean(predTensor, dim=0))

    pred[PID] = {
        "pred_patch": pred_patch_counts,
        "predTensormean": predTensormean
    }

with open('predictions.pkl', 'wb') as f:
    pickle.dump(pred, f)

if True:
    # Load data from pickle file
    with open("predictions.pkl", 'rb') as f:
        pred = pickle.load(f)

    # Extract data from the pred dictionary and flatten it
    data = []
    for PID, values in pred.items():
        row = {"PID": PID}

        # For 'pred_patch'
        for key, count in values['pred_patch'].items():
            row[f"pred_patch_{key}"] = count

        # For 'predTensormean'
        for idx, val in enumerate(values['predTensormean']):
            row[f"predTensormean_{idx}"] = val

        data.append(row)

    # Convert to DataFrame
    df = pd.DataFrame(data).fillna(0)

    # Calculating the ratio for 'pred_patch'
    df["num_patch"] = np.sum(df.filter(like="pred_patch_"), axis=1)
    df["ratio_1_to_all"] = df["pred_patch_1"] / df["num_patch"]
    df["ratio_0_to_all"] = df["pred_patch_0"] / df["num_patch"]

    df["num_peritumor_patch"] = np.sum(df[["pred_patch_1", "pred_patch_0"]], axis=1)
    df["ratio_1_to_peritumor"] = df["pred_patch_1"] / df["num_peritumor_patch"]
    df["ratio_0_to_peritumor"] = df["pred_patch_0"] / df["num_peritumor_patch"]

    df["ratio_1_to_0"] = df["pred_patch_1"] / df["pred_patch_0"]

    # Save the dataframe
    df.to_csv("test.csv")

    # Sorting the DataFrame based on the ratio
    df_sorted = df.sort_values(by="ratio_0_to_all", ascending=False)

    print(df_sorted)

    df_sorted.loc[np.where(df_sorted["ratio_1_to_0"]>1)].shape


for catetype in catetypes:
    visualize_gradcam(
        img_path= "/mnt/data/lyx/IMC/WSI_HE/WSI_Trainging/WSI_crop/B2/15_13.tif" ,
        transform=trnsfrms_val,model=model,
        target_class=int(catetype),device=device,
        save_name=os.path.join("./"+catetype+".png")
        )

# Test the function
test = False
if test:
    folder1 = "/mnt/data/lyx/IMC/WSI_HE/WSI_Data/Training_TIFF_norm/1"
    folder2 = "/mnt/data/lyx/IMC/WSI_HE/WSI_Data/Training_TIFF_norm/3"

    mean1, std1 = compute_mean_std(folder1)
    mean2, std2 = compute_mean_std(folder2)

    trnsfrms_1 = transforms.Compose(
            [
                transforms.Resize(1024), ## 224
                transforms.CenterCrop(1024), ## 224
                transforms.ToTensor(),
                transforms.Normalize(mean = mean1, std = std1)
            ]
        )

    trnsfrms_2 = transforms.Compose(
            [
                transforms.Resize(1024), ## 224
                transforms.CenterCrop(1024), ## 224
                transforms.ToTensor(),
                transforms.Normalize(mean = mean2, std = std2)
            ]
        )

    trnsfrms_x = transforms.Compose(
            [
                transforms.Resize(1024), ## 224
                transforms.CenterCrop(1024), ## 224
                transforms.ToTensor(),
                transforms.Normalize(mean = (0.485, 0.456, 0.406), std = (0.229, 0.224, 0.225))
            ]
        )

    # visualize_batch_effect(folder1, folder2, transform = trnsfrms_val,method='PCA')  # You can change to method='t-SNE' for t-SNE visualization
    ### deal with batch effect
    visualize_batch_effect(folder1, folder2, transform1 = trnsfrms_1, transform2 = trnsfrms_2,method='PCA',sample_fraction = 0.25, save_path="./re-train_1 train_3.png")  # You can change to method='t-SNE' for t-SNE visualization
    visualize_batch_effect(folder1, folder2, transform1 = trnsfrms_x, transform2 = trnsfrms_x,method='PCA',sample_fraction = 0.25, save_path="./train_1 B2.png")  # You can change to method='t-SNE' for t-SNE visualization
