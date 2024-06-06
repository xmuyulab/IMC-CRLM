import os
import pickle
import matplotlib.pyplot as plt
import argparse

def main(args):
    # Base path
    base_path = args.base_path

    # Fetch all directories under base path
    all_dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d))]

    results = {
        'accuracy': [],
        'precision': [],
        'recall': [],
        'f1-score': []
    }

    # Extracting data from each directory
    for d in all_dirs:
        val_file_path = os.path.join(base_path, d, "val.pkl")
        
        if os.path.exists(val_file_path):
            with open(val_file_path, 'rb') as f:
                val_data = pickle.load(f)
                
                metrics = val_data[args.last_epoch]
                for metric, value in metrics.items():
                    results[metric].append(value)

    # Splitting results for plotting
    relabels = ["Iter_"+str(i) for i in range(len(all_dirs))]

    # Plotting
    fig, ax1 = plt.subplots(figsize=(15, 6))

    # Bar plots for each metric
    width = 0.1  # width of a bar
    ind = range(len(all_dirs))  # x locations for the groups

    for idx, (metric, values) in enumerate(results.items()):
        ax1.bar([i + idx * width for i in ind], values, width, alpha=0.6, label=metric)

    ax1.set_xlabel('Directory')
    ax1.set_ylabel('Metric Value')
    ax1.set_xticks([i + 1.5 * width for i in ind])
    ax1.set_xticklabels(relabels, rotation=45)
    ax1.legend()

    # Display the plot
    plt.title('Metrics from val.pkl across different runs')
    plt.tight_layout()    
    plt.savefig(args.output)
    plt.show()

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Visualize metrics from different runs.")
    parser.add_argument("--base_path", default="./weights/", type=str,
                        help="Base path to the directories containing 'val.pkl'.")
    parser.add_argument("--last_epoch", required=True, type=int,
                        help="The epoch number to extract from val.pkl.")
    parser.add_argument("--output", default="test.png", type=str,
                        help="Output path for the saved plot.")
    args = parser.parse_args()

    main(args)

# python robustness_eval.py --last_epoch 5 --base_path "./weights/" --output "./barplot of re-run 20 times.png"
