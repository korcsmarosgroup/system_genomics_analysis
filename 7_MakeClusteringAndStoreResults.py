import sys
import argparse
import pandas as pd


from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial import distance
from scipy.cluster.hierarchy import set_link_color_palette


def parse_args(argv):
    """ Command line interface for the the module """
    parser = argparse.ArgumentParser()

    # Input file path to learn from
    parser.add_argument("-i", "--input",
                        help="<path to the input file name> [mandatory]",
                        dest="input_file",
                        action="store",
                        required=True)
    parser.add_argument("-o", "--output",
                        help="<path to the output file> [mandatory]",
                        dest="output_file",
                        action="store",
                        required=True)

    parser.add_argument("-mt", "--metric",
                        help="Distance metric for clustering. It can be jaccard, hamming or euclidean. Default: hamming",
                        dest="metric",
                        action="store",
                        required=False)
    parser.add_argument("-me", "--method",
                        help="Method of calculating clustering. Default: average.",
                        dest="method",
                        action="store",
                        required=False)

    parser.add_argument("-th", "--treshold",
                        help="Choose the clustering threshold. It does not chooses automatically. Default = 0.4 ",
                        dest="treshold",
                        action="store",
                        required=False)

    parser.add_argument("-k", "--clusters",
                        help="Number of clusters. Default: 4",
                        dest="clusters",
                        action="store",
                        required=False)

    parser.add_argument("-t", "--title",
                        help="Title of the figure. Default: 'Hierarchical Clustering of Ulcerative Colitis Patients'",
                        dest="figtitle",
                        action="store",
                        required=False)

    parser.add_argument("-xl", "--xlabel",
                        help="X label, DEfault: Patient label",
                        dest="xlabel",
                        action="store",
                        required=False)
    parser.add_argument("-p", "--palette",
                        help = "Palette for thg dendogram. By default the matplotlib colours."
                               "You can use 'tol or 'wong' for colour blind friendly colours or a list of hex values.",
                        dest = "palette",
                        action= "store",
                        required = False)
    results = parser.parse_args(argv)

    return results


def make_dendogram(panda_df, file_name, metric="hamming", method="average", tereshold=0.4, k=4, figtitle ="", xlabel="",
                   palette = "default"):
    """
    It makes a dendogram using Hamming distance from the previously setted data using avarage method.
    :param xlabel: X label of the figure
    :param figtitle: title of the clustering figure
    :param k: number of clusters what we want to keep. It can change.
    :param tereshold: The threshold for the figure to start the clustering . It is not necessarily the same as the
    we see on the clustering outcome. You need to be careful with that.
    :param method: The clustering method default is average. For details see:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    :param metric: The clustering metric. The default is hamming, which is good for bitwise clusters.
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
    :param panda_df: the source data frame
    :param file_name: the output file name
    """

    data = panda_df.values
    data = data.T
    Z = linkage(data, method=method, metric=metric)
    Y = distance.pdist(data, metric=metric)
    W = distance.pdist(data.T, metric=metric)
    plt.figure(figsize=(25, 10))
    plt.title(figtitle, fontsize=20, fontname="Arial")
    plt.xlabel(xlabel, fontsize=20, fontname="Arial")
    plt.ylabel(metric.capitalize() +' Distance', fontsize=20, fontname="Arial")

    label_colors = {}
    labels = panda_df.columns.tolist()
    for label in labels:
        print(label)
        label_colors[label] = (0, 0, 0)
    tol_palette = ["#4477AA", "#66CCEE","#228833","#CCBB44","#EE6677","AA3377"]
    wong_palette = ["#E69D00", "#56B3E9",  "#0072B2", "#D55E00", "#CC79A7", "#009E73", "#F0E442"]
    if palette == "tol":
        set_link_color_palette(tol_palette)
    elif palette == "wong":
        set_link_color_palette(wong_palette)
    elif type(palette) == list:
        try:
            set_link_color_palette(palette)
        except:
            print("The given colours are not recognisable. Using default.")
            set_link_color_palette(None)
    else:
        set_link_color_palette(None)

    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=4.0,  # font size for the x axis labels
        color_threshold=tereshold,  # 6, 0.10, #4.5, #150, #0.6,
        labels=panda_df.columns.tolist())
    print(Z[0])
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(label_colors[lbl.get_text()])
    plt.savefig(file_name, dpi=600, format="png", transparent=True)

    clusters = fcluster(Z, k, criterion='maxclust')
    out = open(file_name.replace(".png", "_cluster.txt"), "w")
    out2 = open(file_name.replace(".png", "_distance_matrix.txt"), "w")
    Y = distance.squareform(Y)
    for k in range(len(clusters)):
        out.write(str(clusters[k]) + "\t" + str(panda_df.columns.tolist()[k]) + "\n")
        for i in range(k):
            out2.write(str(panda_df.columns.tolist()[k]) + " " + str(panda_df.columns.tolist()[i]) + " ")
            out2.write(str(1- Y[i, k]))
            out2.write("\n")

    out3 = open(file_name.replace(".png", "protein_distance_matrix.txt"), "w")
    W = distance.squareform(W)
    for k in range(len(panda_df.index.tolist())):
        for i in range(k):
            out3.write(str(panda_df.index.tolist()[k]) + " " + str(panda_df.index.tolist()[i]) + " ")
            out3.write(str(1.0 - (W[i, k])))
            out3.write("\n")
    out.close()
    out2.close()
    out3.close()


def main(argv):
    args = parse_args(argv)
    input_df = pd.read_csv(args.input_file, index_col=0, sep="\t")
    input_df_T = input_df

    if args.metric:
        metric = args.metric
    else:
        metric = "hamming"

    if args.treshold:
        treshold = float(args.treshold)
    else:
        treshold = 0.4

    if args.clusters:
        k = args.clusters
    else:
        k = 4

    if args.figtitle:
        figtitle = args.figtitle
    else:
        figtitle = "Hierarchical Clustering of Ulcerative Colitis Patients"

    if args.xlabel:
        xlabel = args.xlabel
    else:
        xlabel = "Patient Index"

    if args.method:
        method = args.method
    else:
        method = "average"

    if args.palette:
        palette = args.palette
    else:
        palette = None

    make_dendogram(input_df_T, args.output_file, metric=metric, method= method,
                   tereshold=treshold, k=k, figtitle=figtitle, xlabel=xlabel,
                   palette=palette)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))