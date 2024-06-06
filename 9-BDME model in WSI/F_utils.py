import os
import json
import pickle
import random
import sys
import shutil
import cv2
import torch
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from PIL import Image
from tqdm import tqdm
from torchvision import transforms
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

## .tiff to .jpeg
def convert_tif_to_jpeg_in_subdir(input_dir, output_dir, to_type = "PNG"):
    """
    Convert .tif images in subdirectories of a given directory to .jpeg format,
    preserving the subdirectory structure.

    Parameters:
    - input_dir (str): Path to the parent directory containing subdirectories with .tif images.
    - output_dir (str): Base path where the converted .jpeg images should be saved, maintaining the subdirectory structure.

    Returns:
    None
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Iterate over all subdirectories in the given directory
    for subdir in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, subdir)
        
        # Check if it's indeed a directory
        if os.path.isdir(subdir_path):
            # Create corresponding output subdirectory
            output_subdir = os.path.join(output_dir, subdir)
            os.makedirs(output_subdir, exist_ok=True)
            
            images = os.listdir(subdir_path)
            images = tqdm(images, file=sys.stdout)

            # Convert each .tif file in this subdirectory
            for file in images:
                if file.lower().endswith('.tif'):
                    input_file_path = os.path.join(subdir_path, file)
                    output_file_path = os.path.join(output_subdir, os.path.splitext(file)[0] + '.' + to_type.lower())
                    with Image.open(input_file_path) as img:
                        # print(img.shape)
                        img.convert('RGB').save(output_file_path, to_type)


# Function to retrieve category data
def get_category_data(path, categories):
    data = {}
    for category in categories:
        # Assuming the data in each category directory is a list of files
        category_path = os.path.join(path, category)
        data[category] = os.listdir(category_path)
    return data

# Split and save function
class AddGaussianNoise(object):
    def __init__(self, mean=0., std=0.01):
        self.std = std
        self.mean = mean
        
    def __call__(self, tensor):
        return tensor + torch.randn(tensor.size()) * self.std + self.mean
    
    def __repr__(self):
        return self.__class__.__name__ + '(mean={0}, std={1})'.format(self.mean, self.std)

def split_and_save(data, origin_path, save_path, fraction=0.9):
    train_data = {}
    val_data = {}

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    for category, items in data.items():
        # Shuffle the data
        random.shuffle(items)
        
        # Split the data into fraction% train and (1-fraction)% validation
        n_train = int(fraction * len(items))
        
        train_data[category] = items[:n_train]
        val_data[category] = items[n_train:]
        
        # Save the training and validation data to the specified paths
        train_path = os.path.join(save_path, 'train', category)
        val_path = os.path.join(save_path, 'val', category)
        
        os.makedirs(train_path, exist_ok=True)
        os.makedirs(val_path, exist_ok=True)
        
        # Copy items to train and val folders instead of moving
        for item in train_data[category]:
            shutil.copy2(os.path.join(origin_path, category, item), os.path.join(train_path, item))
        
        for item in val_data[category]:
            shutil.copy2(os.path.join(origin_path, category, item), os.path.join(val_path, item))
                
    return None


# Function to generate data
def generate_training_data(images_path, categories = ['0', '1']):
    for cat in categories:
        img_list = os.listdir(os.path.join(images_path, cat))
        outfile = os.path.join(images_path, cat)

        if not os.path.exists(outfile):
            os.mkdir(outfile)

        img_list = tqdm(img_list, file=sys.stdout)

        for img_name in img_list:
            im = Image.open(os.path.join(images_path, cat, img_name))

            # Original
            new_name = '0_' + img_name
            im.save(os.path.join(outfile, new_name))

            # Horizontal Flip
            new_name = '1_' + img_name
            new_im = transforms.RandomHorizontalFlip(p=1)(im)
            new_im.save(os.path.join(outfile, new_name))

            # Vertical Flip
            new_name = '2_' + img_name
            new_im = transforms.RandomVerticalFlip(p=1)(im)
            new_im.save(os.path.join(outfile, new_name))

            # Rotation
            new_name = '3_' + img_name
            new_im = transforms.RandomVerticalFlip(p=1)(im)
            new_im = transforms.RandomHorizontalFlip(p=1)(new_im)
            new_im.save(os.path.join(outfile, new_name))
            
            # Gaussian Noise
            new_name = '4_' + img_name
            noisy_transform = transforms.Compose([
            transforms.ToTensor(),
            AddGaussianNoise(0., 0.001),  # you can adjust the mean and std deviation
            transforms.ToPILImage(),
            ])
            new_im = noisy_transform(im)
            new_im.save(os.path.join(outfile, new_name))
            
            os.remove(os.path.join(images_path, cat, img_name))

    return None

def plot_data_loader_image(data_loader):
    batch_size = data_loader.batch_size
    plot_num = min(batch_size, 4)

    json_path = './class_indices.json'
    assert os.path.exists(json_path), json_path + " does not exist."
    json_file = open(json_path, 'r')
    class_indices = json.load(json_file)

    for data in data_loader:
        images, labels = data
        for i in range(plot_num):
            # [C, H, W] -> [H, W, C]
            img = images[i].numpy().transpose(1, 2, 0)
            # 反Normalize操作
            img = (img * [0.229, 0.224, 0.225] + [0.485, 0.456, 0.406]) * 255
            label = labels[i].item()
            plt.subplot(1, plot_num, i+1)
            plt.xlabel(class_indices[str(label)])
            plt.xticks([])  # 去掉x轴的刻度
            plt.yticks([])  # 去掉y轴的刻度
            plt.imshow(img.astype('uint8'))
        plt.show()


def write_pickle(list_info: list, file_name: str):
    with open(file_name, 'wb') as f:
        pickle.dump(list_info, f)


def read_pickle(file_name: str) -> list:
    with open(file_name, 'rb') as f:
        info_list = pickle.load(f)
        return info_list


def getModelSize(model):
    param_size = 0
    param_sum = 0
    for param in model.parameters():
        param_size += param.nelement() * param.element_size()
        param_sum += param.nelement()
    buffer_size = 0
    buffer_sum = 0
    for buffer in model.buffers():
        buffer_size += buffer.nelement() * buffer.element_size()
        buffer_sum += buffer.nelement()
    all_size = (param_size + buffer_size) / 1024 / 1024
    print('Model Size: {:.3f}MB'.format(all_size))
    return None

def unnormalize(tensor, mean, std):
    """Revert normalization on a tensor."""
    for t, m, s in zip(tensor, mean, std):
        t.mul_(s).add_(m)
    return tensor


def visualize_gradcam(img_path, transform, mean_, std_, model, target_class, device, save_name: str = None):
    ## Load image
    assert os.path.exists(img_path), "file: '{}' dose not exist.".format(img_path)
    img = Image.open(img_path)
    img_tensor = transform(img)

    # Based on the class, get the specific normalization statistics
    normalize = transforms.Normalize(mean=mean_, std=std_)
    img_tensor = normalize(img_tensor)

    # expand batch dimension
    input_tensor = torch.unsqueeze(img_tensor, dim=0).to(device)

    gradcam = model.module.get_gradcam(input_tensor, target_class=target_class)  # Assuming you want the heatmap for class 0

    # Normalize the Grad-CAM between 0 and 1
    gradcam = (gradcam - gradcam.min()) / (gradcam.max() - gradcam.min())

    # Resize Grad-CAM to image size
    gradcam = cv2.resize(gradcam[0, 0].cpu().detach().numpy(), (img_tensor.shape[1], img_tensor.shape[2]))

    # Convert image tensor to numpy and bring it to range [0, 255]
    img_tensor = unnormalize(img_tensor,mean = mean_, std = std_)
    img_np = img_tensor.squeeze().permute(1, 2, 0).cpu().detach().numpy() * 255.0
    img_np = img_np.astype(np.uint8)

    # Convert Grad-CAM into colormap
    heatmap = cv2.applyColorMap(np.uint8(255 * gradcam), cv2.COLORMAP_JET)

    # Superimpose the Grad-CAM on the image
    superimposed_img = heatmap * 0.4 + img_np  # Adjust the factor (0.4 here) to control the intensity of the heatmap
    superimposed_img = superimposed_img.astype(np.uint8)
    image = Image.fromarray(superimposed_img)
    image.save(save_name)


def center_crop_img(img: np.ndarray, size: int):
    h, w, c = img.shape

    if w == h == size:
        return img

    if w < h:
        ratio = size / w
        new_w = size
        new_h = int(h * ratio)
    else:
        ratio = size / h
        new_h = size
        new_w = int(w * ratio)

    img = cv2.resize(img, dsize=(new_w, new_h))

    if new_w == size:
        h = (new_h - size) // 2
        img = img[h: h+size]
    else:
        w = (new_w - size) // 2
        img = img[:, w: w+size]

    return img

def load_images_from_folder(folder,transform):
    """Load images from a folder into a list."""
    images = []
    for filename in os.listdir(folder):
        if filename.endswith(".jpg") or filename.endswith(".png") or filename.endswith(".jpeg") or filename.endswith(".tif"):
            img = Image.open(os.path.join(folder, filename))
            img = transform(img)
            img_array = np.array(img).flatten()  # Convert image to flattened array
            images.append(img_array)
    return images

def visualize_batch_effect(folder1, folder2, transform1, transform2, method='PCA',sample_fraction=1,save_path="./test.png"):
    """Visualize batch effect using PCA or t-SNE."""
    
    # Load images
    images1 = load_images_from_folder(folder1,transform1)
    images2 = load_images_from_folder(folder2,transform2)
    
    # If sample_fraction < 1.0, sample the points randomly
    if 0 < sample_fraction < 1.0:
        num_samples1 = int(sample_fraction * len(images1))
        num_samples2 = int(sample_fraction * len(images2))

        indices1 = np.random.choice(len(images1), num_samples1, replace=False)
        indices2 = np.random.choice(len(images2), num_samples2, replace=False)

        images1 = [images1[i] for i in indices1]
        images2 = [images2[i] for i in indices2]

    # Concatenate and label the datasets
    all_images = np.array(images1 + images2)
    labels = ['Folder1']*len(images1) + ['Folder2']*len(images2)
    
    # Apply dimension reduction
    if method == 'PCA':
        pca = PCA(n_components=2)
        points = pca.fit_transform(all_images)
    elif method == 't-SNE':
        tsne = TSNE(n_components=2, random_state=0)
        points = tsne.fit_transform(all_images)
    else:
        raise ValueError("Method should be either 'PCA' or 't-SNE'")
    
    # Plot the 2D points
    plt.figure(figsize=(10, 7))
    # Plot the points from folder1
    for i in range(len(images1)):
        plt.scatter(points[i, 0], points[i, 1], color='r', label='Folder1' if i == 0 else "")
        
    # Plot the points from folder2
    for i in range(len(images1), len(images1) + len(images2)):
        plt.scatter(points[i, 0], points[i, 1], color='b', label='Folder2' if i == len(images1) else "")
    
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.title(f'Batch Effect Visualization using {method}')
    plt.legend()
    plt.grid(True)
    plt.savefig(save_path)
    plt.show()

    return None

def compute_mean_std(folder_path):
    """
    Compute mean and std for all images in the given folder.
    
    Parameters:
    - folder_path: path to the folder containing images.

    Returns:
    - mean, std: mean and standard deviation of images in the folder.
    """

    # Define transformation to convert images to tensors
    to_tensor = transforms.ToTensor()

    # List to store all image tensors
    tensors = []

    # Loop over all files in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        
        # Open image and convert to tensor
        with Image.open(file_path) as img:
            tensor = to_tensor(img)

        # If there are already tensors processed, check if the current tensor's shape matches the first tensor's shape.
        # If it doesn't match, print the problematic file and continue to the next iteration without appending the tensor.
        if tensors and tensor.shape != tensors[0].shape:
            print(f"Skipped file {file_path} due to mismatched shape {tensor.shape}")
            continue

        # Append the tensor to the list
        tensors.append(tensor)

    # Stack tensors along a new dimension
    stacked_tensors = np.stack(tensors, axis=0)
    
    # Calculate mean and std along the (height, width, channel) dimensions
    mean = np.mean(stacked_tensors, axis=(0, 2, 3))
    std = np.std(stacked_tensors, axis=(0, 2, 3))

    return mean, std

def visualize_stats(stats_file_path, save_path = "./status.png"):
    # Load statistics from pickle file
    with open(stats_file_path, 'rb') as f:
        stats = pickle.load(f)

    # Extracting data for plotting
    categories = list(stats.keys())
    means = np.array([stats[category]['mean'] for category in categories])
    stds = np.array([stats[category]['std'] for category in categories])

    # Set up the figure and axes
    fig, ax = plt.subplots(2, 1, figsize=(10, 8))

    # Bar width and positions
    bar_width = 0.2
    positions = np.arange(3)  # 3 channels

    # Plot mean for each category for each channel
    for i, category in enumerate(categories):
        ax[0].bar(positions + i * bar_width, means[i], bar_width, label=category)

    # Plot std for each category for each channel
    for i, category in enumerate(categories):
        ax[1].bar(positions + i * bar_width, stds[i], bar_width, label=category)

    # Set the properties for the plots
    channel_names = ['Channel 1', 'Channel 2', 'Channel 3']

    ax[0].set_title('Mean values for each channel')
    ax[0].set_xticks(positions + bar_width * (len(categories) - 1) / 2)  # Center the xticks
    ax[0].set_xticklabels(channel_names)
    ax[0].set_ylabel('Mean Value')
    ax[0].set_xlabel('Channel')
    ax[0].legend(title='Category')

    ax[1].set_title('Standard Deviation values for each channel')
    ax[1].set_xticks(positions + bar_width * (len(categories) - 1) / 2)  # Center the xticks
    ax[1].set_xticklabels(channel_names)
    ax[1].set_ylabel('Standard Deviation Value')
    ax[1].set_xlabel('Channel')
    ax[1].legend(title='Category')

    # Show the plot
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()

    return None

def plot_heatmap_one_class(data, title, save_path):
    """
    Plot and save a heatmap based on the given data.
    
    Args:
    - data (np.array): 2D array containing heatmap values.
    - title (str): Title for the heatmap plot.
    - save_path (str): File path to save the heatmap plot as a PDF.
    """
    plt.figure(figsize=(10, 10))  # Adjust the size as per your requirement
    sns.heatmap(data, cmap='viridis', square=True)
    plt.title(title)
    plt.savefig(save_path)
    plt.show()

def plot_class_heatmap_one_class(df, categories, PID, save_path):
    """
    Plot and save heatmaps for each category based on DataFrame values.
    
    Args:
    - df (pd.DataFrame): DataFrame containing patch coordinates and probabilities.
    - categories (list): List of categories (e.g., ['0', '1', '2']).
    - PID (str): WSI ID.
    - save_path (str): Directory path to save the heatmap plots.
    """
    max_row = df['Row'].max() + 1
    max_col = df['Col'].max() + 1
    
    df['Row'] = df['Row'].astype(int)
    df['Col'] = df['Col'].astype(int)

    for category in categories:
        # Create an empty heatmap array
        heatmap_ = -np.ones((max_row, max_col))
        
        # Populate the heatmap using DataFrame values
        for index, row in df.iterrows():
            heatmap_[int(row['Row']), int(row['Col'])] = row[f'Class_{category}_Prob']

        # Generate the heatmap plot and save it
        plot_heatmap_one_class(
            data=heatmap_, 
            title=f'Class {category} Probability Heatmap for WSI: {PID}', 
            save_path=os.path.join(save_path, f"{PID} heatmap of category {category}.png")
        )

def plot_heatmap_multiclass(data, categories, title, save_path):
    """
    Plot and save a heatmap based on the given data.
    
    Args:
    - data (np.array): 2D array containing heatmap values.
    - categories (list): List of categories to define the color scale.
    - title (str): Title for the heatmap plot.
    - save_path (str): File path to save the heatmap plot as a PDF.
    """
    plt.figure(figsize=(10, 10))
    
    # Set -1 class to gray, and use the tab10 colormap for the other classes
    cmap = ['white'] + sns.color_palette("Paired", len(categories))
    
    sns.heatmap(data, cmap=cmap, square=True, vmin=-1, vmax=len(categories)-1, cbar_kws=dict(ticks=np.arange(-1, len(categories))), annot=False)
    
    colorbar = plt.gca().collections[0].colorbar
    colorbar.set_ticks(np.arange(-1, len(categories)))
    colorbar.set_ticklabels(['-1'] + categories)
    
    plt.title(title)
    plt.savefig(save_path)
    plt.show()


def plot_class_heatmap_multiclass(df, categories, PID, save_path):
    """
    Plot and save heatmap with patches colored by their highest probability class.
    
    Args:
    - df (pd.DataFrame): DataFrame containing patch coordinates and probabilities.
    - categories (list): List of categories (e.g., ['0', '1', '2']).
    - PID (str): WSI ID.
    - save_path (str): Directory path to save the heatmap plots.
    """
    max_row = df['Row'].max() + 1
    max_col = df['Col'].max() + 1
    df['Row'] = df['Row'].astype(int)
    df['Col'] = df['Col'].astype(int)

    # Create an empty heatmap array with default value -1
    heatmap_ = -np.ones((max_row, max_col))
    
    # Assign each patch to a class based on highest probability
    for index, row in df.iterrows():
        class_probs = [row[f'Class_{category}_Prob'] for category in categories]
        heatmap_[int(row['Row']), int(row['Col'])] = np.argmax(class_probs)

    # Generate the heatmap plot and save it
    plot_heatmap_multiclass(
        data=heatmap_, 
        categories=categories,
        title=f'Class Distribution Heatmap for WSI: {PID}', 
        save_path=os.path.join(save_path, f"{PID} combined heatmap.pdf")
    )