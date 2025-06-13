import pandas as pd
import matplotlib.pyplot as plt
# Load the two datasets
apo_df = pd.read_csv("apo.csv")
holo_df = pd.read_csv("holo.csv")

apo_df.head(), holo_df.head()


# Extract residue names and RMSF values
residues = apo_df.columns[1:]  # Exclude the 'label' column
apo_rmsf = apo_df.iloc[0, 1:].values.astype(float)*10 # convert to angstroms
holo_rmsf = holo_df.iloc[0, 1:].values.astype(float)*10 # convert to angstroms

# Define bar positions
x = range(len(residues))

# Define residue groups
glp1_residues = {"PHE12", "VAL16", "LEU20"}
glp1r_residues = {"LEU142", "TYR145", "LYS202"}

# Create the bar plot
fig, ax = plt.subplots(figsize=(12, 6))
bar_width = 0.35

bars1 = ax.bar([i - bar_width / 2 for i in x], apo_rmsf, width=bar_width, label='$\mathbf{\it{apo}}$: GLP-1-R + GLP-1', color='orange')
bars2 = ax.bar([i + bar_width / 2 for i in x], holo_rmsf, width=bar_width, label='$\mathbf{\it{holo}}$: GLP-1-R + GLP-1 + PAM', color='blue')

#set limits
ax.set_ylim(0,2)

# Axis labeling and formatting
ax.set_xlabel('Residue',labelpad=20,fontweight="bold",fontsize=18)
ax.set_ylabel('RMSF ($\AA$)',labelpad=5,fontweight="bold",fontsize=18)
#ax.set_title('RMSF Comparison Between $\mathbf{\it{apo}}$ and $\mathbf{\it{holo}}$ simulations',pad=20,fontweight="bold")
ax.set_xticks(list(x))
ax.set_xticklabels(residues,fontsize=14)
ax.legend(loc="upper left",prop={"style":"italic",'size':14})



plt.show()

# Save plots at A0 and A1 resolutions
#a0_dpi = 300  # High resolution for A0
#a1_dpi = 200  # Medium resolution for A1

a0_path = "rmsf_comparison.png"
fig.savefig(a0_path, dpi=300, bbox_inches='tight')

plt.close(fig)
