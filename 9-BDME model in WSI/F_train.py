import sys
import matplotlib.pyplot as plt
from tqdm import tqdm

import torch
import torch.nn.functional as F

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

def train_one_epoch(model, optimizer, data_loader, device, epoch):
    model.train()
    loss_function = torch.nn.CrossEntropyLoss()
    mean_loss = torch.zeros(1).to(device)
    optimizer.zero_grad()

    # 验证样本总个数
    total_num = len(data_loader.dataset)
    # 用于存储预测正确的样本个数
    sum_num = torch.zeros(1).to(device)

    data_loader = tqdm(data_loader, file=sys.stdout)

    for step, data in enumerate(data_loader):
        # step, data = next(enumerate(data_loader))
        images, labels = data
        images, labels = images.to(device), labels.to(device)

        pred = model(images)

        sum_num += torch.eq(torch.max(pred, dim=1)[1], labels).sum()

        loss = loss_function(pred, labels)
        loss.backward()
        mean_loss = (mean_loss * step + loss.detach()) / (step + 1)  # update mean losses

        data_loader.desc = "[epoch {}] mean loss {}".format(epoch, round(mean_loss.item(), 3))

        if not torch.isfinite(loss):
            print('WARNING: non-finite loss, ending training ', loss)
            sys.exit(1)

        optimizer.step()
        optimizer.zero_grad()

    return mean_loss.item(), sum_num.item() / total_num


@torch.no_grad()
def evaluate(model, data_loader, device):
    model.eval()

    # data_loader = tqdm(data_loader, file=sys.stdout)

    y_pred = []
    y_ture = []

    for _, data in enumerate(data_loader):        

        images, labels = data
        pred = model(images.to(device))
        pred = torch.max(pred, dim=1)[1]

        y_ture.extend(labels.cpu().detach().numpy().tolist())
        y_pred.extend(pred.cpu().detach().numpy().tolist())

    accuracy = accuracy_score(y_ture, y_pred)
    precision = precision_score(y_ture, y_pred,average='macro')
    recall = recall_score(y_ture, y_pred,average='macro')
    f1score = f1_score(y_ture,y_pred,average='macro')

    return accuracy, precision, recall, f1score

@torch.no_grad()
def infer(model, data_loader, device):
    model.eval()

    # List to store predictions
    all_preds = []

    data_loader = tqdm(data_loader, file=sys.stdout)

    for step, data in enumerate(data_loader):
        # step, data = next(enumerate(data_loader))

        images = data
        pred = F.softmax(model(images.to(device)), dim=1)  # make sure you specify the dimension in softmax
        all_preds.append(pred.cpu())  # append the batch predictions to the list

    # Concatenate all predictions into a single tensor
    all_preds_tensor = torch.cat(all_preds, dim=0)
    
    return all_preds_tensor

def plot_losses(train_losses, val_losses, save_path):
    """Plot training and validation losses."""
    plt.figure(figsize=(10, 6))
    plt.plot(train_losses, label='Training Loss', color='blue')
    plt.plot(val_losses, label='Validation Loss', color='red')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Training vs Validation Loss')
    plt.legend()
    plt.grid(True)
    plt.savefig(save_path)
    plt.show()

    return None