import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Set up paths relative to the project root
file_path = os.path.realpath(__file__)
src_directory = os.path.dirname(file_path)
project_root = os.path.dirname(src_directory)
data_file = os.path.join(project_root, "data", "summary.csv")
output_dir = os.path.join(project_root, "assets", "images")

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

if not os.path.exists(data_file):
    print(f"Data file not found: {data_file}")
    exit(1)

df = pd.read_csv(data_file)

# --- Proportion Pass/Fail ---
plt.figure(figsize=(7, 5))
counts = df['Within5?'].value_counts()
counts = counts.reindex(['YES', 'NO'], fill_value=0)
total = counts.sum()
bars = plt.bar(counts.index, counts.values, color=['#2ecc71', '#e74c3c']) # More professional colors
plt.title('Proportion of Stations Within 5 Degrees', fontsize=12)
plt.xlabel('Within 5 Degrees?', fontsize=10)
plt.ylabel('Number of Stations', fontsize=10)

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + (total * 0.02), f'{yval}\n({yval/total*100:.1f}%)', ha='center', va='bottom', fontsize=12)

plt.ylim(0, total * 1.15)
plt.savefig(os.path.join(output_dir, 'pass_fail_summary.png'), dpi=300, bbox_inches='tight')
plt.close()

# --- Distribution of Errors ---
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.hist(df['AbsOffset'], bins=20, color='#3498db', edgecolor='black', alpha=0.7)
plt.xlabel('Absolute Offset (Degrees)', fontsize=10)
plt.ylabel('Frequency (Log Scale)', fontsize=10)
plt.yscale('log')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.title('Histogram of Orientation Offsets')

plt.subplot(1, 2, 2)
plt.boxplot(df['AbsOffset'], vert=True, patch_artist=True, 
            boxprops=dict(facecolor='#3498db', alpha=0.7), 
            medianprops=dict(color='black', linewidth=2))
plt.ylabel('Absolute Offset (Degrees)', fontsize=10)
plt.xticks([1], ['All Stations'])
plt.title('Distribution of Offsets')

plt.suptitle('Distribution of Orientation Errors across Network', fontsize=14)
plt.savefig(os.path.join(output_dir, 'error_distribution.png'), dpi=300, bbox_inches='tight')
plt.close()

# --- Failing Stations Detail ---
failing_stations = df[df['Within5?'] == 'NO'].sort_values(by='AbsOffset', ascending=False)
n_failing = len(failing_stations)

if n_failing > 0:
    plt.figure(figsize=(10, max(6, n_failing * 0.4)))
    labels = failing_stations['Network'] + '.' + failing_stations['Station']
    plt.barh(labels, failing_stations['AbsOffset'], color='#e67e22', alpha=0.8)
    plt.title(f'Orientation Offsets for Failing Stations (>5° Error)', fontsize=12)
    plt.xlabel('Absolute Offset (Degrees)', fontsize=10)
    plt.ylabel('Station ID', fontsize=10)
    plt.gca().invert_yaxis()
    plt.grid(axis='x', linestyle='--', alpha=0.35)

    for i, row in failing_stations.iterrows():
        plt.text(row['AbsOffset'] + 0.5, list(labels).index(f"{row['Network']}.{row['Station']}"), 
                 f"{row['AbsOffset']:.1f}°", va='center', ha='left', fontsize=10)
    
    plt.savefig(os.path.join(output_dir, 'failing_stations_detail.png'), dpi=300, bbox_inches='tight')
    plt.close()
