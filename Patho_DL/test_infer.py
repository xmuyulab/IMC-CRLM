import os
import numpy as np
import pandas as pd
import random
import pickle
import argparse
from collections import Counter
import torch
from torch.utils.data import DataLoader
from torchvision import transforms
from F_dataset import *
from F_utils import *
from F_model import *
from F_train import *

do_GradCAM_on = "test"
saveBasePath="/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/test_GradCAM"
ori_img_path="/mnt/data/lyx/IMC/WSI_HE/WSI_Testing/WSI_crop_norm"
seed = 619
WSIPath = '/mnt/data/lyx/IMC/WSI_HE/WSI_Trainging/WSI_crop'
stats_file="/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/stats.pkl"
resize = 1024
center_crop = 1024
num_hidden_features = 96
num_classes = 4
heatmap_save_path = '/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/heatmap'
pred_file = '/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/predictions.pkl'
weights_dir = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/weights/2023-09-11_08-10-49"
drop_rate = 0.5
batch_size = 32
num_workers = 8
pretrain_path = "./pretrainModel/ctranspath.pth"
encode_features = 768
prediction_layers = 3

epochs = 200

def GradCam(args):
    # Grad-CAM
    device = 'cuda' if torch.cuda.is_available() else "cpu"

    model = ClassiModel(pretrain_path=pretrain_path, encode_features=encode_features, prediction_layers=prediction_layers, num_hidden_features=num_hidden_features, num_classes=num_classes, drop_rate=drop_rate)
    model = nn.DataParallel(model).to(device)
    model.load_state_dict(torch.load(os.path.join(weights_dir, "model-{}.pth".format(epochs))))

    catetypes = [str(x) for x in range(num_classes)]

    trnsfrms_val = transforms.Compose([
        transforms.Resize(resize),
        transforms.CenterCrop(center_crop),
        transforms.ToTensor()
    ])

    with open(stats_file, 'rb') as f:
        transform_status = pickle.load(f)
        
    for key in transform_status.keys():
        transform_status[key]["mean"] = np.array([0.485, 0.456, 0.406])
        transform_status[key]["std"] = np.array([0.229, 0.224, 0.225])

    if not os.path.exists(saveBasePath):
        os.mkdir(saveBasePath)

    if do_GradCAM_on == "val":
        for dirtype in ["val"]:

            for catetype in catetypes:
                imgs = os.listdir(os.path.join(ori_img_path, dirtype, catetype))

                for img in imgs:
                    if not os.path.exists(os.path.join(saveBasePath,dirtype,catetype,img)):
                        os.makedirs(os.path.join(saveBasePath,dirtype,catetype,img))

                    for target_class in catetypes:
                        visualize_gradcam(
                                img_path=os.path.join(ori_img_path,dirtype,catetype,img) ,
                                transform=trnsfrms_val, mean_ = transform_status[catetype]["mean"], std_ = transform_status[catetype]["std"],
                                model=model,
                                target_class=int(target_class),device=device,
                                save_name=os.path.join(saveBasePath,dirtype,catetype,img,"pred_"+target_class+".png"))

    if do_GradCAM_on == "test":
        for PID in os.listdir(ori_img_path):
            imgs = os.listdir(os.path.join(ori_img_path, PID))

            for img in imgs:
                if not os.path.exists(os.path.join(saveBasePath,PID,img)):
                    os.makedirs(os.path.join(saveBasePath,PID,img))

                for target_class in catetypes:
                    visualize_gradcam(
                            img_path=os.path.join(ori_img_path,PID,img) ,
                            transform=trnsfrms_val, mean_ = np.array([0.485, 0.456, 0.406]), std_ = np.array([0.229, 0.224, 0.225]),
                            model=model,
                            target_class=int(target_class),device=device,
                            save_name=os.path.join(saveBasePath,PID,img,"pred_"+target_class+".png"))

def main(args):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.cuda.manual_seed_all(seed)

    if not os.path.exists(WSIPath):
        raise ValueError(f"The specified WSIPath {WSIPath} does not exist.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    nw = min([os.cpu_count(), batch_size if batch_size > 1 else 0, num_workers])

    # Load the model
    model = ClassiModel(
        pretrain_path=pretrain_path,
        encode_features=encode_features,
        prediction_layers=prediction_layers,
        num_hidden_features=num_hidden_features,
        num_classes=num_classes,
        drop_rate=drop_rate)
    model = nn.DataParallel(model).to(device)
    model.load_state_dict(torch.load(os.path.join(weights_dir, f"model-{epochs}.pth")))

    PIDs = os.listdir(WSIPath)
    pred = {}

    PID = "B2"

    for PID in PIDs:
        infer_images_dir = os.path.join(WSIPath, PID)
        infer_image_files = os.listdir(infer_images_dir)
        infer_images_path = [os.path.join(infer_images_dir, x) for x in infer_image_files]

        # mean_, std_ = compute_mean_std(infer_images_dir)
        mean_ = np.array([0.485, 0.456, 0.406])
        std_ = np.array([0.229, 0.224, 0.225])
        
        trnsfrms_val = transforms.Compose([
            transforms.Resize(resize),
            transforms.CenterCrop(center_crop),
            transforms.ToTensor(),
            transforms.Normalize(mean=mean_, std=std_)
        ])

        infer_dataset = MyDataSet_Infer(images_path=infer_images_path, transform=trnsfrms_val)
        infer_dataloader = DataLoader(infer_dataset, batch_size=batch_size, shuffle=False, pin_memory=True, num_workers=nw)

        predTensor = infer(model=model, data_loader=infer_dataloader, device=device)

        coords = [tuple(map(int, fname.split('.')[0].split('_'))) for fname in infer_image_files]
        data = []
        for (row, col), patch_pred_probs in zip(coords, predTensor):
            patch_pred_probs = patch_pred_probs.tolist()
            data.append([row, col] + patch_pred_probs)

        categories = [str(x) for x in range(num_classes)]
        columns_ = ["Col","Row"]
        columns_.extend([f"Class_{i}_Prob" for i in categories])

        df = pd.DataFrame(data, columns=columns_)

        save_path = os.path.join(heatmap_save_path, PID)
        os.makedirs(save_path, exist_ok=True)
        plot_class_heatmap(df=df, categories=categories, PID=PID, save_path=save_path)

        pred[PID] = {"pred_patches": df}

    with open(pred_file, 'wb') as f:
        pickle.dump(pred, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments for training and inference.')

    parser.add_argument('--seed', type=int, default=619, help='Seed value for initializing random generators.')
    
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
    parser.add_argument('--epochs', type=int, default=30, help='Number of times the model will iterate over the entire training dataset.')

    parser.add_argument('--weights_dir', default="./weights/", help='Directory where model weights are saved.')

    parser.add_argument('--WSIPath', default="/mnt/data/lyx/IMC/WSI_HE/WSI_Trainging/WSI_crop/", help='Path to Whole Slide Images (WSI) for inference.')
    parser.add_argument('--heatmap_save_path', default="./pred_heatmap", help='Path to save predictions heatmap after inference.')
    parser.add_argument('--pred_file', default="predictions.pkl", help='File to save predictions after inference.')

    args = parser.parse_args()
    
    main(args)

Test = False
if Test:
    pred_file = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/predictions.pkl"
    with open(pred_file, 'rb') as f:
        pred = pickle.load(f)
        
    all_classes = []
    PIDs = [x for x in pred.keys()]
    for PID in PIDs:
        df_ = pred[PID]["pred_patches"]
        # Extract the relevant columns (containing probabilities)
        prob_columns = [col for col in df_.columns if 'Class' in col]

        # Find the class with the max probability for each row
        df_['Max_Class'] = df_[prob_columns].idxmax(axis=1)

        # Find the max probability value for each row
        df_['Max_Prob'] = df_[prob_columns].max(axis=1)
        
        pred[PID]["pred_classes"] = dict(Counter(df_['Max_Class']))
        
        # Class1_Class0_ = pred[PID]["pred_classes"]["Class_1_Prob"] / pred[PID]["pred_classes"]["Class_0_Prob"]
        # pred[PID]["pred_classes"]["Class1_Class0"] = Class1_Class0_
        # print(f"The {PID} class1:class0 = {Class1_Class0_}")
        
        all_classes.extend([(PID, key, value) for key, value in pred[PID]["pred_classes"].items() if "Class1_Class0" not in key])

    # Convert the list into a dataframe
    all_classes_df = pd.DataFrame(all_classes, columns=['PID', 'Class', 'Count'])
    all_classes_df = all_classes_df.pivot(index='PID', columns='Class', values='Count')
        
    # with open(pred_file, 'wb') as f:
    #     pickle.dump(pred, f)
        
        
        
