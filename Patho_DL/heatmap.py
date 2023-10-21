import os
import pickle
from F_dataset import *
from F_utils import *
from F_model import *
from F_train import *

predPath = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/val_predictions.pkl"
with open(predPath, 'rb') as f:
    pred = pickle.load(f)

heatmap_save_path = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/val_heatmap"

categories = [str(x) for x in range(4)]

for PID in pred.keys():
    save_path = os.path.join(heatmap_save_path, PID)
    os.makedirs(save_path, exist_ok=True)
    df  = pred[PID]["pred_patches"]
    plot_class_heatmap_multiclass(df=df, categories=categories, PID=PID, save_path=save_path)


## Predict with differential grade
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

def plot_pred_with_clinical_group(pred_filepath, clinical_df, img_savepath):
    pred_class = pd.read_csv(pred_filepath,index_col=0)
    pred_class.index = pred_class.index.str.replace('-', '_')

    idx = pred_class.index.intersection(clinical_df.index)
    pred_class = pred_class.loc[idx,]
    clinical_df = clinical_df.loc[idx,]

    pred_class["Class_1_0_Prob"] = pred_class["Class_1_Prob"] / pred_class["Class_0_Prob"]
    pred_class["Differential_grade"] = clinical_df.loc[:,"Differential_grade"]
    pred_class["Lymph_positive"] = clinical_df.loc[:,"Lymph_positive"]
    pred_class["VI"] = clinical_df.loc[:,"VI"]

    # Sample dataframe
    df = pd.DataFrame(pred_class)
    df = df.dropna(subset=['Class_1_0_Prob'])

    # Box Plot with statistical annotation
    fig, axs = plt.subplots(2, 3, figsize=(20, 15))  # Create a 3x1 grid

    # Differential_grade
    sns.boxplot(x=df['Differential_grade'], y=df['Class_0_Prob'], palette="deep", ax=axs[0,0])
    group0 = df[df['Differential_grade'] == 0]['Class_0_Prob']
    group1 = df[df['Differential_grade'] == 1]['Class_0_Prob']
    stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    axs[0,0].set_title('Box plot of Class_0_Prob by Differential_grade')
    axs[0,0].set_ylabel('Class_0_Prob value')
    axs[0,0].set_xlabel('Differential_grade')
    axs[0,0].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    sns.boxplot(x=df['Differential_grade'], y=df['Class_1_Prob'], palette="deep", ax=axs[0,1])
    group0 = df[df['Differential_grade'] == 0]['Class_1_Prob']
    group1 = df[df['Differential_grade'] == 1]['Class_1_Prob']
    stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    axs[0,1].set_title('Box plot of Class_1_Prob by Differential_grade')
    axs[0,1].set_ylabel('Class_1_Prob value')
    axs[0,1].set_xlabel('Differential_grade')
    axs[0,1].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    sns.boxplot(x=df['Differential_grade'], y=df['Class_1_0_Prob'], palette="deep", ax=axs[0,2])
    group0 = df[df['Differential_grade'] == 0]['Class_1_0_Prob']
    group1 = df[df['Differential_grade'] == 1]['Class_1_0_Prob']
    stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    axs[0,2].set_title('Box plot of Class_1_0_Prob by Differential_grade')
    axs[0,2].set_ylabel('Class_1_0_Prob value')
    axs[0,2].set_xlabel('Differential_grade')
    axs[0,2].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    # Lymph_positive
    sns.boxplot(x=df['Lymph_positive'], y=df['Class_0_Prob'], palette="deep", ax=axs[1,0])
    group0 = df[df['Lymph_positive'] == 0]['Class_0_Prob']
    group1 = df[df['Lymph_positive'] == 1]['Class_0_Prob']
    stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    axs[1,0].set_title('Box plot of Class_0_Prob by Lymph_positive')
    axs[1,0].set_ylabel('Class_0_Prob value')
    axs[1,0].set_xlabel('Lymph_positive')
    axs[1,0].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    sns.boxplot(x=df['Lymph_positive'], y=df['Class_1_Prob'], palette="deep", ax=axs[1,1])
    group0 = df[df['Lymph_positive'] == 0]['Class_1_Prob']
    group1 = df[df['Lymph_positive'] == 1]['Class_1_Prob']
    stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    axs[1,1].set_title('Box plot of Class_1_Prob by Lymph_positive')
    axs[1,1].set_ylabel('Class_1_Prob value')
    axs[1,1].set_xlabel('Lymph_positive')
    axs[1,1].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    sns.boxplot(x=df['Lymph_positive'], y=df['Class_1_0_Prob'], palette="deep", ax=axs[1,2])
    group0 = df[df['Lymph_positive'] == 0]['Class_1_0_Prob']
    group1 = df[df['Lymph_positive'] == 1]['Class_1_0_Prob']
    stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    axs[1,2].set_title('Box plot of Class_1_0_Prob by Lymph_positive')
    axs[1,2].set_ylabel('Class_1_0_Prob value')
    axs[1,2].set_xlabel('Lymph_positive')
    axs[1,2].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    # sns.boxplot(x=df['VI'], y=df['Class_0_Prob'], palette="deep", ax=axs[1,0])
    # group0 = df[df['VI'] == 0]['Class_0_Prob']
    # group1 = df[df['VI'] == 1]['Class_0_Prob']
    # stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    # axs[1,0].set_title('Box plot of Class_0_Prob by VI')
    # axs[1,0].set_ylabel('Class_0_Prob value')
    # axs[1,0].set_xlabel('VI')
    # axs[1,0].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    # sns.boxplot(x=df['VI'], y=df['Class_1_Prob'], palette="deep", ax=axs[1,1])
    # group0 = df[df['VI'] == 0]['Class_1_Prob']
    # group1 = df[df['VI'] == 1]['Class_1_Prob']
    # stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    # axs[1,1].set_title('Box plot of Class_1_Prob by VI')
    # axs[1,1].set_ylabel('Class_1_Prob value')
    # axs[1,1].set_xlabel('VI')
    # axs[1,1].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    # sns.boxplot(x=df['VI'], y=df['Class_1_0_Prob'], palette="deep", ax=axs[1,2])
    # group0 = df[df['VI'] == 0]['Class_1_0_Prob']
    # group1 = df[df['VI'] == 1]['Class_1_0_Prob']
    # stat, p_value = mannwhitneyu(group0, group1, alternative='two-sided')
    # axs[1,2].set_title('Box plot of Class_1_0_Prob by VI')
    # axs[1,2].set_ylabel('Class_1_0_Prob value')
    # axs[1,2].set_xlabel('VI')
    # axs[1,2].annotate(f'p-value = {p_value:.3e}', xy=(0.3, 0.95), xycoords='axes fraction')

    # Adjust the layout
    plt.tight_layout()

    # Save the figure
    plt.savefig(img_savepath)

    # Display the merged figure
    plt.show()

    return None


val_clinical = pd.read_csv("/mnt/data/lyx/IMC/clinical.csv",index_col=0)
val_clinical.head

plot_pred_with_clinical_group(pred_filepath = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/val_all_classes_df.csv", clinical_df  = val_clinical, img_savepath = "./val_pred_with_diff_lympho_boxplot.png")

test_clinical = pd.read_csv("/mnt/data/lyx/IMC/WSI_HE/WSI_Testing/WSI_Test_clinical.csv",index_col=0)
test_clinical.head

plot_pred_with_clinical_group(pred_filepath = "/mnt/data/lyx/IMC/WSI_HE/HE_classification_2/test_all_classes_df.csv", clinical_df  = test_clinical, img_savepath = "./test_pred_with_diff_lympho_boxplot.png")
