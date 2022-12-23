import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from matplotlib import pyplot as plt
from umap import UMAP
import os

my_file = open("sample_names", "r")
content = my_file.read().replace('\n', '')
SAMPLE = content.split(",")

df = pd.read_csv("output/sorc/clusters.tsv", sep="\t")
df["sum_prob"] = df.iloc[:, 3:].sum(axis=1)
df.iloc[:, 3:-1] = df.iloc[:, 3:-1].div(df.sum_prob, axis=0) * 4
x = df.iloc[:, 3:-1]
x = StandardScaler().fit_transform(x)
df['new_labels'] = df["assignment"]
df.loc[(df["status"] == 'doublet'), 'new_labels'] = 'doublet'
df.loc[(df["status"] == 'unassigned'), 'new_labels'] = 'unassigned'

umap = UMAP(n_components=2, init='random', random_state=0)
proj = umap.fit_transform(x)
df_umap = pd.DataFrame(data = proj, columns = ["UMAP1", "UMAP2"])
df_final_umap = pd.concat([df_umap, df[["new_labels"]]], axis = 1)
palette_colors = sns.color_palette('tab10')
palette_dict = {'3': palette_colors[0], 'unassigned': palette_colors[1], 'doublet': palette_colors[2], '1': palette_colors[3], '2': palette_colors[4], '4': palette_colors[7], '5': palette_colors[6], '0': palette_colors[5]}
plt.figure(figsize=(9,8))
plot = sns.scatterplot(data=df_final_umap, x="UMAP1", y="UMAP2", hue="new_labels", palette=palette_dict)
fig = plot.get_figure()
fig.savefig("../pictures/" + SAMPLE[0] + "_umap_after_souporcell.png", dpi = 300) 

