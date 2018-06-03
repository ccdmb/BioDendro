"""
"""

import argparse

import plotly

from BioDendro.preprocess import MGFRecord
from BioDendro.preprocess import get_csv_record
from BioDendro.preprocess import remove_redundancy
#from BioDendro.plot import Dendrogram


def pipeline(
        mgf,
        components,
        neutral=False,
        cutoff=0.6,
        bin_threshold=8e-4,
        clustering_method="jaccard",
        matrix="",
        processed="",
        out_html="",
        results_dir="results",
        width=900,
        height=400,
        **kwargs
        ):

    #Open the trigger data <file>.msg
    with open(mgf, 'r') as handle:
        rec = get_record(handle)


    #Open the sample list <file>.csv
    with open(components, 'r') as handle:
        ndic = get_csv_record(handle, rec)

    #Now remove redundancy and print best trigger ion list
    table = remove_redundancy(ndic, neutral)

    # Write out an excel file too
    table.to_excel(processed, index=False)

    """
    data = Dendrogram(
        table,
        bin_threshold=bin_threshold,
        cutoff=cutoff,
        clustering_method=clustering_method
        )"""

    data = Dendrogram.from_xlsx(processed)
    data.clusterize()
    data.generate_linkage()
    k = data.visualize(filename=out_html, x=width, y=height)
    return

def main():
    parser = argparse.ArgumentParser(
        description="Process MGF file into appropriate input for biodendro."
        )

    parser.add_argument(
        "mgf",
        help="MGF input file (file1.mgf)",
        type=str
        )
    parser.add_argument(
        "components",
        help="Listed components file (file2.txt)",
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
        type=float,
        default=0.6
        )

    parser.add_argument(
        "-b", "--bin-threshold",
        dest="bin_threshold",
        type=float,
        default=8e-4
        )

    parser.add_argument(
        "-d", "--cluster-method",
        dest="clustering_method",
        default="jaccard",
        choices=["jaccard", "braycurtis"],
        )

    parser.add_argument(
        "-m", "--matrix",
        default="full_matrix.xlsx",
        )

    parser.add_argument(
        "-p", "--processed",
        default="out.xlsx",
        help=("Path to write output to.")
        )

    parser.add_argument(
        "-o", "--output",
        dest="out_html",
        default="simple_dendrogram.html",
        help=("Path to write output to.")
        )

    parser.add_argument(
        "-r", "--results-dir",
        dest="results_dir",
        help=""
        )

    parser.add_argument(
        "-x", "--width",
        dest="width_px",
        help="",
        type=int,
        default=900
        )
    parser.add_argument(
        "-y", "--height",
        dest="height_px",
        help="",
        type=int,
        default=400
        )


    args = parser.parse_args()
    pipeline(**args.__dict__)

    return
