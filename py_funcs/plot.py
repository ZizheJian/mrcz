from matplotlib import pyplot as plt

def plot_curves(index:int,curves:list):
    plt.figure(figsize=(10,6))
    for label,br_list,psnr_list in curves:
        plt.plot(br_list,psnr_list,marker='o',label=label)
    plt.xlabel("Bit Rate")
    plt.ylabel("PSNR")
    plt.legend()
    plt.xlim(left=0)
    plt.tight_layout()
    plt.savefig(f"history/fig_{index:04d}.png")
    plt.close()
