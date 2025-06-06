import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

nj=20
#fout="ukbb-s0-3methods.pdf"
fout="test-pm-95.pdf"
dfs=[]
#for k in ["bvsr","glmnet","combined"]:
#    for j in ["modelx","koscreen","reko-npc100"]:
#        fin = f"power.{j}.{k}.s0"
#        df1 = pd.read_csv(fin, delim_whitespace=True, header=None)
#        dfs.append(df1)
for j in ["pm90","pm95","pm99"]:
        fin = f"{j}/power2.modelx.txt"
        df1 = pd.read_csv(fin, delim_whitespace=True, header=None)
        dfs.append(df1)
for j in ["susie","koscreen","reko"]:
        fin = f"pm90/power2.{j}.txt"
        df1 = pd.read_csv(fin, delim_whitespace=True, header=None)
        dfs.append(df1)
#dfs.append(df1)
df=pd.concat(dfs, ignore_index=True)
print(df.shape)
df.columns=['V1','V2','V3']


# Number of data blocks (assuming 4 rows per block)
blocks = df.shape[0] // 20

# Create figure with two subplots (Power and Type I Error) with specific height ratio
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 12), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

# Set color palette for the blocks
# First two blocks will have distinct colors
#base_colors = sns.color_palette("Set2", n_colors=2)  # First two blocks
base_colors=["orange"]*3+["magenta"]*3+["pink"]*3+["green","blue","red"]
#base_colors=["orange","blue","red"]*3
#    remaining_colors = sns.color_palette("cubehelix", n_colors=5)  # Remaining 10 blocks (same type, different hues)

my_linestyle=['-','--',':']*3 +['-']*3
#my_linestyle=[':']*3+['--']*3+['-']*3
#mylabels = ["ModelX:glmnet", "KoScreen:glmnet", "ReKo:glmnet", "ModelX:fastBVSR", "KoScreen:fastBVSR", "ReKo:fastBVSR", "ModelX:combined", "Koscreen:combined", "ReKo:combined"]
#mylabels = ["ModelX"]*3+["KoScreen"]*3+["ReKo"] * 3
mylabels = ["glment", "fastBVSR", "combined"] *3 + ["SuSiE","KS", "ReKo"]  
#mylabels = ["0.90","0.95","0.99"]*3
#mylabels = ["0.90"]*3+["0.95"]*3+["0.99"]*3
print(mylabels)
#mylist = [0,3,6,1,4,7,2,5,8]
mylist = range(0,6)
print(mylist)

for i in [0,1,2,3,4,5,6,7,8,9,10,11]:
    # Get the current block (4 rows)
    block_data = df.iloc[nj*i : nj*(i+1)]
    
    # Calculate Power and Type I Error
    alpha = block_data['V1'] 
    power = block_data['V2'] 
    type_i_error = block_data['V3'] 
    
    color = base_colors[i]
    
    # Plot Power on the first subplot (ax1)
    ax1.plot(alpha, power, label=mylabels[i], color=color, linewidth=2,marker=None, linestyle=my_linestyle[i], markersize=6)
    
    # Plot Type I Error on the second subplot (ax2)
    ax2.plot(alpha, type_i_error, label=mylabels[i], color=color, linewidth=2,marker=None, linestyle=my_linestyle[i], markersize=6)


# Customize Power plot (ax1)
#    ax1.set_title('Power and realized FDR', fontsize=18)
ax1.set_xlabel('Nominal FDR', fontsize=18)
ax1.set_ylabel('Power', fontsize=18)
ax1.set_ylim(0.0, 0.6)
ax1.set_xlim(-0.01, 0.22)
ax1.grid(True)
ax1.tick_params(axis='both', which='major', labelsize=16); 
#ax1.legend(title=None, loc='lower right', frameon=False,fontsize=16,ncol=2)


# Customize Type I Error plot (ax2)
#ax2.set_title('Type I Error vs. FDR Threshold', fontsize=14)
#    ax2.xticks(fontsize=16); 
#    ax2.yticks(fontsize=16); 
ax2.set_xlabel('Nominal FDR', fontsize=18)
ax2.set_ylabel('Realized FDR', fontsize=18)
ax2.set_ylim(-0.0, 0.3)  # Limit the Type I Error plot range to 0 to 0.3
ax1.set_xlim(-0.01, 0.22)
ax2.grid(True)
ax2.axline((0.01,0.01),(0.2,0.2),color='gray',linewidth=4)
#ax2.legend(title=None, loc='upper left', fontsize=14,ncol=3)
ax2.tick_params(axis='both', which='major', labelsize=16); 


handles, labels = ax2.get_legend_handles_labels()
leg=ax2.legend(handles,labels, ncol=4,
#                loc='lower right', 
                loc='upper left', 
#                bbox_to_anchor=(0.5, -0.15),
#                columnspacing=1.0,
#                handlelength=3,
                fontsize=14,
                frameon=False)

fig.canvas.draw()

legend_width = leg.get_window_extent(fig.canvas.get_renderer()).width
column_width = legend_width / 4
print(column_width)

# Calculate the left position of the legend in figure coordinates
leg_pos = leg.get_window_extent(fig.canvas.get_renderer())
leg_pos = leg_pos.transformed(fig.transFigure.inverted())
left_x = leg_pos.x0 + 0.0
print(leg_pos.x0,leg_pos.x1) 
print(leg_pos.y0,leg_pos.y1) 
column_width=(leg_pos.x1-leg_pos.x0)/4
plt.figtext(left_x + column_width*0.5, leg_pos.y1 - 0.08,
        'ModelX:0.90', ha='center', fontweight='normal',fontsize=14,color="orange")
plt.figtext(left_x + column_width*1.5, leg_pos.y1 - 0.08,
        'ModelX:0.95', ha='center', fontweight='normal',fontsize=14,color="magenta")
plt.figtext(left_x + column_width*2.5, leg_pos.y1 - 0.08,
        'ModelX:0.99', ha='center', fontweight='normal',fontsize=16,color="pink")

# Adjust layout to avoid overlap
plt.tight_layout()

# Save the plot as a PDF file
plt.savefig(fout, format="pdf")

# Close the plot to free up memory
plt.close()

