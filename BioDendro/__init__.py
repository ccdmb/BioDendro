"""
"""

__name__ = "BioDendro"
__version__ = "0.0.1"


import argparse

import plotly

from BioDendro.preprocess import MGF
from BioDendro.preprocess import SampleRecord
from BioDendro.preprocess import split_msms_title
from BioDendro.preprocess import remove_redundancy
from BioDendro.cluster import Tree


def pipeline(
        mgf_path,
        components_path,
        neutral=False,
        cutoff=0.6,
        bin_threshold=8e-4,
        clustering_method="jaccard",
        processed="processed.xlsx",
        results_dir="results",
        out_html="simple_dendrogram.html",
        width=900,
        height=400,
        quiet=False,
        ):
    """ Runs the default BioDendro pipeline. """

    if quiet:
        printer = lambda *x: None
    else:
        printer = lambda *s: print(*s)

    printer((
        "Running {} v{}\n\n".format(__name__, __version__)
        "- input mgf file = {}".format(mgf_path)
        "- input components file = {}".format(components_path)
        "- neutral = {}\n".format(neutral)
        "- cutoff = {}\n".format(cutoff)
        "- bin_threshold = {}\n".format(bin_threshold)
        "- clustering_method = {}\n".format(clustering_method)
        "- output processed file = {}\n".format(processed)
        "- output results directory = {}\n".format(results_dir)
        "- output html dendrogram = {}\n".format(out_html)
        "- dendrogram figure width = {}\n".format(width)
        "- dendrogram figure height = {}\n".format(height)
        "\n"
    ))

    #Open the trigger data <file>.msg
    with open(mgf, 'r') as handle:
        mgf = MGF.parse(handle)

    # Customised MGF title handler.
    # TODO: This title filter will fail for some mgf title fields.
    for rec in mgf.records:
        rec.title = split_msms_title(rec.title)

    #Open the sample list <file>.csv
    with open(components_path, 'r') as handle:
        components = SampleRecord.parse(handle)

    #Now remove redundancy and print best trigger ion list
    printer("Processing inputs")
    table = remove_redundancy(components, mgf, neutral)

    # Write out an excel file too
    table.to_excel(processed, index=False)

    printer("Binning and clustering\nThis may take some time...")
    tree = Tree(bin_threshold, clustering_method, cutoff)
    tree.fit(table)

    printer("Writing per-cluster summaries")
    tree.write_summaries(path=results_dir)

    printer("Writing output html dendrogram")
    _ = tree.iplot(filename=out_html, x=width, y=height)

    printer("Finished")
    return

def main():
    parser = argparse.ArgumentParser(
        description="Run the BioDendro pipeline."
        )

    parser.add_argument(
        "mgf",
        dest="mgf_path",
        help="MGF input file.",
        type=str
        )
    parser.add_argument(
        "components",
        dest="components_path",
        help="Listed components file.",
        type=str
        )

    parser.add_argument(
        "-n", "--neutral",
        help="Apply neutral loss.",
        action="store_true",
        default=False
        )

    parser.add_argument(
        "-c", "--cutoff",
        help="Distance threshold for selecting clusters from tree.",
        type=float,
        default=0.6
        )

    parser.add_argument(
        "-b", "--bin-threshold",
        help="Threshold for binning m/z values prior to clustering.",
        dest="bin_threshold",
        type=float,
        default=8e-4
        )

    parser.add_argument(
        "-d", "--cluster-method",
        help="The distance metric used during tree construction.",
        dest="clustering_method",
        default="jaccard",
        choices=["jaccard", "braycurtis"],
        )

    parser.add_argument(
        "-p", "--processed",
        default="processed.xlsx",
        help="Path to write preprocessed output to."
        )

    parser.add_argument(
        "-o", "--output",
        dest="out_html",
        default="simple_dendrogram.html",
        help="Path to write interactive html plot to."
        )

    parser.add_argument(
        "-r", "--results-dir",
        dest="results_dir",
        default=None,
        help=("Directory to write per-cluster plots and table to." 
              "Default is to use 'results_20180606112200' where the number is"
              "the current datetime.")
        )

    parser.add_argument(
        "-x", "--width",
        dest="width_px",
        help="The width of the dendrogram plot in pixels.",
        type=int,
        default=900
        )

    parser.add_argument(
        "-y", "--height",
        dest="height_px",
        help="The height of the dendrogram plot in pixels.",
        type=int,
        default=400
        )

    parser.add_argument(
        "-q", "--quiet",
        help="Suppress status notifications written to stdout.",
        action="store_true",
        type=bool,
        default=False
        )


    args = parser.parse_args()
    pipeline(**args.__dict__)
    return
