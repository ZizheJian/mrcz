from matplotlib import pyplot as plt

def plot_curves(index:int,curves:list):
    plt.figure(figsize=(6,6))
    for label,br_dict,psnr_dict in curves:
        plt.plot(list(br_dict.values()),list(psnr_dict.values()),marker='o',label=label)
    plt.xlabel("Bit Rate")
    plt.ylabel("PSNR")
    plt.legend()
    plt.xlim(left=0)
    plt.tight_layout()
    plt.savefig(f"history/fig_{index:04d}.png")
    plt.close()
