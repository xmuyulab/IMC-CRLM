import os
import pickle
import random
import argparse

from F_utils import *

def main(args):
    random.seed(args.seed)
    
    if args.convert:
        convert_tif_to_jpeg_in_subdir(input_dir=args.input_directory, output_dir=args.ConverSavePath, to_type="JPEG")

    ## split and generate data set
    # Get category list
    categories = os.listdir(args.ConverSavePath)
    data = get_category_data(args.ConverSavePath, categories)

    # Compute mean and std for each category and save in a pickle file
    if args.RGB_normalize:
        stats = {}
        for category in categories:
            folder_path = os.path.join(args.ConverSavePath, category)
            mean, std = compute_mean_std(folder_path)
            stats[category] = {'mean': mean, 'std': std}

        with open('stats.pkl', 'wb') as f:
            pickle.dump(stats, f)

        visualize_stats('stats.pkl')

    split_and_save(data=data, origin_path=args.ConverSavePath, save_path=args.SplitSavePath, fraction=args.split_fraction)
    
    generate_training_data(images_path=args.GeneratePath, categories=args.lager_categories)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process images and generate training data")
    
    # Add arguments
    parser.add_argument('--seed', type=int, default=619, help='Random seed')
    parser.add_argument('--convert', action='store_true', help='Flag to convert images from .tiff to .jpeg')
    parser.add_argument('--input_directory', type=str, default='/mnt/data/lyx/IMC/WSI_HE/WSI_Data/Training_TIFF', help='Input directory for .tiff images')
    parser.add_argument('--RGB_normalize', action='store_true', help='Flag to perform normalization on different dictionary.')
    parser.add_argument('--ConverSavePath', type=str, default='/mnt/data/lyx/IMC/WSI_HE/WSI_Data/ori_data_to_train', help='Path to save converted images')
    parser.add_argument('--SplitSavePath', type=str, default='/mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data', help='Path to save split datasets')
    parser.add_argument('--GeneratePath', type=str, default='/mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data/train', help='Path to generate training data')
    parser.add_argument('--lager_categories', nargs='+', default=['0','1'], help='List of categories to generate training data')
    parser.add_argument('--split_fraction', type=float, default=0.8, help='Fraction of data to use for training')

    # Parse arguments and call the main function
    args = parser.parse_args()
    main(args)

# nohup \
# python initial.py --seed 619 --RGB_normalize \
# --ConverSavePath /mnt/data/lyx/IMC/WSI_HE/WSI_Data/ori_data_to_train \
# --SplitSavePath /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data \
# --GeneratePath /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data/train \
# --lager_categories "0" "1" "2" \
# --split_fraction 0.8 \
# &

# nohup \
# python initial.py --seed 619 \
# --ConverSavePath /mnt/data/lyx/IMC/WSI_HE/WSI_Data/ori_data_to_train \
# --SplitSavePath /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data \
# --GeneratePath /mnt/data/lyx/IMC/WSI_HE/WSI_Data/split_train_data/train \
# --lager_categories "0" "1" \
# --split_fraction 0.8 \
# &