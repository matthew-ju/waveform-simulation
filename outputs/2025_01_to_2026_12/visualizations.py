import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


df = pd.read_csv("summary.csv")

plt.figure(figsize=(7, 5))
counts = df['Within5?'].value_counts()
counts = counts.reindex(['YES', 'NO'], fill_value=0)
total = counts.sum()
bars = plt.bar(counts.index, counts.values, color=['green', 'red'])
plt.title('Proportion Within 5 Degrees', fontsize=10)
plt.xlabel('Within 5 Degrees?', fontsize=8)
plt.ylabel('Number of Stations', fontsize=8)

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 5, f'{yval}\n({yval/total*100:.2f}%)', ha='center', va='bottom', fontsize=15)

plt.ylim(0, total * 1.1)
plt.savefig('pass_fail_bar.pdf')
plt.close()

# ==============================================================================

plt.subplot(1, 2, 1)
plt.hist(df['AbsOffset'], bins=20, color='0.85', edgecolor='black')
plt.xlabel('Absolute Offset in Degrees', fontsize=8)
plt.ylabel('Frequency (Log)', fontsize=8)
plt.yscale('log')
plt.grid(axis='y', linestyle='--')

plt.subplot(1, 2, 2)
plt.boxplot(df['AbsOffset'], vert=True, patch_artist=True, boxprops=dict(facecolor='0.85'), medianprops=dict(color='black'))
plt.xlabel('All Stations', fontsize=8)
plt.ylabel('Absolute Offset (Degrees)', fontsize=8)
plt.xticks([1], ['All Stations'])
plt.suptitle('Distribution of Errors', fontsize=10)
plt.savefig('dist_err_histo_box.pdf')
plt.close()

# ==============================================================================
failing_stations = df[df['Within5?'] == 'NO'].sort_values(by='AbsOffset', ascending=False)
n_failing = len(failing_stations)

if n_failing > 0:
    plt.figure(figsize=(10, max(5, n_failing * 0.4)))
    labels = failing_stations['Network'] + '.' + failing_stations['Station']
    plt.barh(labels, failing_stations['AbsOffset'], color='0.85')
    plt.title(f'Absolute Offset for the {n_failing} Failing Stations (Offset > 5 Degrees)', fontsize=10)
    plt.xlabel('Absolute Offset in Degrees', fontsize=8)
    plt.ylabel('Station ID', fontsize=8)
    plt.gca().invert_yaxis()
    plt.grid(axis='x', linestyle='--', alpha=0.35)

    for i, row in failing_stations.head(min(n_failing, n_failing)).iterrows():
        plt.text(row['AbsOffset'], labels.loc[i], f" {row['AbsOffset']:.1f}°", va='center', ha='left', fontsize=9, color='black')
    plt.savefig('fail_sta_bar.pdf')
    plt.close()


