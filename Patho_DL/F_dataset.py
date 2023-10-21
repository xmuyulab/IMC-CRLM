import os

from PIL import Image
import torch
from torch.utils.data import Dataset
from torchvision import transforms

def read_split_data(root: str, num_categories: int):
    train_images_path = []
    train_images_label = []
    val_images_path = []
    val_images_label = []

    # Iterate over the categories
    for cat in range(num_categories):
        train_path = os.path.join(root, f'train/{cat}/')
        val_path = os.path.join(root, f'val/{cat}/')

        # Train data
        for i in os.listdir(train_path):
            train_images_path.append(os.path.join(train_path, i))
            train_images_label.append(cat)

        # Validation data
        for i in os.listdir(val_path):
            val_images_path.append(os.path.join(val_path, i))
            val_images_label.append(cat)

    print("{} images for training.".format(len(train_images_path)))
    print("{} images for validation.".format(len(val_images_path)))

    return train_images_path, train_images_label, val_images_path, val_images_label


# class MyDataSet_Train(Dataset):

#     def __init__(self, images_path: list, images_class: list, transform=None):
#         self.images_path = images_path
#         self.images_class = images_class
#         self.transform = transform

#     def __len__(self):
#         return len(self.images_path)

#     def __getitem__(self, item):
#         img = Image.open(self.images_path[item])
#         if img.mode != 'RGB':
#             raise ValueError("image: {} isn't RGB mode.".format(self.images_path[item]))
#         label = self.images_class[item]

#         if self.transform is not None:
#             img = self.transform(img)

#         return img, label

#     @staticmethod
#     def collate_fn(batch):
#         # 官方实现的default_collate可以参考
#         # https://github.com/pytorch/pytorch/blob/67b7e751e6b5931a9f45274653f4f653a4e6cdf6/torch/utils/data/_utils/collate.py
#         images, labels = tuple(zip(*batch))

#         images = torch.stack(images, dim=0)
#         labels = torch.as_tensor(labels)
#         return images, labels

class MyDataSet_Train(Dataset):

    def __init__(self, images_path: list, images_class: list, transform=None, stats_file=None):
        self.images_path = images_path
        self.images_class = images_class
        self.transform = transform
        self.stats = stats_file

    def __len__(self):
        return len(self.images_path)

    def __getitem__(self, item):
        img = Image.open(self.images_path[item])
        if img.mode != 'RGB':
            raise ValueError("image: {} isn't RGB mode.".format(self.images_path[item]))
        label = self.images_class[item]
        
        # Based on the class, get the specific normalization statistics
        mean = self.stats[str(label)]['mean']
        std = self.stats[str(label)]['std']

        normalize = transforms.Normalize(mean=mean, std=std)
        
        # If the transform is not None, apply it
        if self.transform:
            img = self.transform(img)
        
        # Apply the class-specific normalization
        img = normalize(img)

        return img, label

    @staticmethod
    def collate_fn(batch):
        # 官方实现的default_collate可以参考
        # https://github.com/pytorch/pytorch/blob/67b7e751e6b5931a9f45274653f4f653a4e6cdf6/torch/utils/data/_utils/collate.py
        images, labels = tuple(zip(*batch))

        images = torch.stack(images, dim=0)
        labels = torch.as_tensor(labels)
        return images, labels

class MyDataSet_Infer(Dataset):

    def __init__(self, images_path: list, transform=None):
        self.images_path = images_path
        self.transform = transform

    def __len__(self):
        return len(self.images_path)

    def __getitem__(self, item):
        img = Image.open(self.images_path[item])
        if img.mode != 'RGB':
            raise ValueError("image: {} isn't RGB mode.".format(self.images_path[item]))

        if self.transform is not None:
            img = self.transform(img)

        return img

    @staticmethod
    def collate_fn(batch):
        # 官方实现的default_collate可以参考
        # https://github.com/pytorch/pytorch/blob/67b7e751e6b5931a9f45274653f4f653a4e6cdf6/torch/utils/data/_utils/collate.py
        images, labels = tuple(zip(*batch))

        images = torch.stack(images, dim=0)
        return images