import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("contacts_all.csv")

# Convert units
df['Time_ns'] = df['Time_ps'] / 1000
df['Total_Surface_Area_A2'] = df['Total_Surface_Area_nm2'] * 100  # 1 nm^2 = 100 Å^2

# Filter for first 500 ns
df_filtered = df[df['Time_ns'] <= 500]

#Get 10 highest surface contact area
df_sorted = df_filtered.sort_values(by='Total_Surface_Area_A2',ascending=False)
df_top10 = df_sorted.head(10)
print(df_top10)

# Get unique simulation IDs
sim_ids = sorted(df_filtered['simulation_id'].unique())

# Determine common axis limits
xlim = (0, 500)
ylim=(400,1000)
'''
ylim = (
    df_filtered['Total_Surface_Area_A2'].min()*0.9,
    df_filtered['Total_Surface_Area_A2'].max()*1.1
)
'''

# Create 2x2 subplots with no spacing
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0, hspace=0)
# Plot data
for idx, sim_id in enumerate(sim_ids):
    ax = axes[idx // 2, idx % 2]
    sim_data = df_filtered[df_filtered['simulation_id'] == sim_id]
    ax.plot(sim_data['Time_ns'], sim_data['Total_Surface_Area_A2'],label=f"Replicate {sim_id}")
    df_i=df_top10[df_top10['simulation_id']==sim_id]
    if len(df_i)>0:
        ax.scatter(df_i['Time_ns'], df_i['Total_Surface_Area_A2'],color='red',s=10,marker='X',label="Selected for docking")
    ax.legend(prop={"style":"italic",'size':14},loc="lower left")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    #ax.set_title(f"Simulation ID: {sim_id}")
    
    # Global x and y labels for the entire figure
fig.text(0.5, 0.04, 'Time (ns)', ha='center', va='center', fontsize=18,fontweight="bold")
fig.text(0.04, 0.5, 'Contact Surface Area ($\AA^2$)', ha='center', va='center', rotation='vertical', fontsize=18,fontweight="bold")


plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.show()

a0_path = "contact_surface.png"
fig.savefig(a0_path, dpi=300, bbox_inches='tight')

plt.close(fig)
