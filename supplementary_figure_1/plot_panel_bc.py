import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

base = "PATH/scripts/06_revision/sandra_traelseq/src/models"
data = {
	i: pd.read_csv(f"{base}/{i}", sep="\t")
	for i in os.listdir(base) if i.endswith(".log")
	}

fig, axes = plt.subplots(2, 3, figsize=(10, 6))

for log in [key for key in data.keys() if "schick" not in key]:#data.keys():
	sns.lineplot(data=data[log], x='epoch', y='loss', ax=axes[0, 0], alpha=0.1, color="black")
	sns.lineplot(data=data[log], x='epoch', y='r2_score', ax=axes[0, 1], alpha=0.1, color="black")
	sns.lineplot(data=data[log], x='epoch', y='root_mean_squared_error', ax=axes[0, 2], alpha=0.1, color="black")
	sns.lineplot(data=data[log], x='epoch', y='val_loss', ax=axes[1, 0], alpha=0.1, color="black")
	sns.lineplot(data=data[log], x='epoch', y='val_r2_score', ax=axes[1, 1], alpha=0.1, color="black")
	sns.lineplot(data=data[log], x='epoch', y='val_root_mean_squared_error', ax=axes[1, 2], alpha=0.1, color="black")

plt.suptitle(f"Training History of Human TrAEL-Seq models (HCT116)\nn={len([key for key in data.keys() if "schick" not in key])}", fontsize=16)

axes[0, 0].set_xlabel("")
axes[0, 1].set_xlabel("")
axes[0, 2].set_xlabel("")
axes[1, 0].set_xlabel("")
axes[1, 1].set_xlabel("Epoch")
axes[1, 2].set_xlabel("")

axes[0, 0].set_ylabel("Training")
axes[0, 1].set_ylabel("")
axes[0, 2].set_ylabel("")
axes[1, 0].set_ylabel("Validation")
axes[1, 1].set_ylabel("")
axes[1, 2].set_ylabel("")

axes[0, 0].set_title("Loss")
axes[0, 1].set_title("R2 Score")
axes[0, 2].set_title("RMSE")
axes[1, 0].set_title("")
axes[1, 1].set_title("")
axes[1, 2].set_title("")

plt.tight_layout()
plt.savefig("./training_history_hct116.pdf", dpi=300)
plt.show()

epochs = []
for log in [key for key in data.keys() if "schick" not in key]:
	epochs.append(data[log]['epoch'].max())

sns.histplot(epochs, binwidth=1, color="black")
plt.axvline(x=300, color='red', linestyle='dashed', linewidth=1)
plt.axvline(x=0, color='darkblue', linestyle='dotted', linewidth=1)
plt.title(f"Distribution of Training Epochs for HCT116 models\nn={len([key for key in data.keys() if "schick" not in key])}", fontsize=16)
plt.xlabel("Number of Training Epochs")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig("./training_epochs_hct116.pdf", dpi=300)
plt.show()
