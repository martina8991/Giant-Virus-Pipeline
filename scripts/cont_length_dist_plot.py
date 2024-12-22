import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np


"""This script takes a text file with contig_id and the respective contig length as input and outputs a contig length distribution plot. To prevent outliers from skewering the distribution plot, values with +/- 2 std are removed."""

contig_len_file = sys.argv[1]
output_dist_plot = sys.argv[2]


def get_contig_lengths(contig_len_file):
    """convert contig lengths from input file into a list"""
    contig_lengths = []
    with open(contig_len_file, "r") as file:
        for line in file:
            contig_id, length = line.strip().split()
            contig_lengths.append(int(length))
    return contig_lengths


def remove_outliers(contig_length_list):
    """remove contigs lengths with std +/-2"""
    contig_lengths_without_outlier = []
    for x in contig_length_list:
        z_score = (x - mean) / std
        if abs(z_score) < 2:
            contig_lengths_without_outlier.append(x)
    return contig_lengths_without_outlier


def comma_formatter(x, pos):
    """Add commas to the distribution plot."""
    return f"{int(x):,}"

contig_lengths = get_contig_lengths(contig_len_file)
mean = np.mean(np.array(contig_lengths))
std = np.std(np.array(contig_lengths))
contig_lengths_without_outlier = remove_outliers(contig_lengths)

plt.hist(contig_lengths_without_outlier, bins=100, color="skyblue", edgecolor="black")
plt.subplots_adjust(left=0.15)
plt.gca().xaxis.set_major_formatter(FuncFormatter(comma_formatter))
plt.gca().yaxis.set_major_formatter(FuncFormatter(comma_formatter))
plt.xlabel("Contig Length")
plt.ylabel("Frequency")
plt.title("Contig Length Distribution Plot")
plt.savefig(output_dist_plot)
