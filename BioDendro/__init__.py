"""
"""

__name__ = "BioDendro"
__version__ = "0.0.1"

import os
from os.path import join as pjoin
import argparse
from datetime import datetime

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
    results_dir=None,
    out_html="simple_dendrogram.html",
    width=900,
    height=800,
    quiet=False,
    scaling=False,
    filtering=False,
    eps=0.6,
    mz_tol=0.002,
    retention_tol=5,
    **kwargs
):
    """ Runs the default BioDendro pipeline.

    Minimal example:
    >>> tree = BioDendro.pipeline("MSMS.mgf", "component_list.txt")

    Required parameters:
    mgf_path           Name of .mgf file 
    components_path    Name of .txt file

    Optional parameters:
    parameter          description
                       default value
                       user options

    neutral            convert MSMS spectra to neutral loss spectra
                       False
                       True or False

    cutoff             the y-axis value to cut the dendogram and form clusters
                       0.6
                       From 0 to 1

    bin_threshold      maximum difference between two consecutive masses before a new mass bin is created
                       0.0008
                       Typical range 0.00001-1 Dalton

    clustering method  Distance metric used for clustering
                       "jaccard"
                       "jaccard" or "braycurtis"

    width              width of dendogram output in pixels
                       900
                       Recommended maximum 1200

    height             height of dendogram output in pixels. Branch labels are included in this dimension
                       800
                       Recommended maximum 1200

    quiet              suppress pipeline messages
                       False
                       True or False

    scaling            Highest m/z within an MSMS spectra is normalised to 1 and all other masses are scaled to it
                       False
                       True or False

    filtering          Ions within an MSMS spectra are removed based on scaled intensity
                       False
                       True or False

    eps                Scaled value for which ions below are removed
                       0.6
                       0.0 - 1.0

    mz_tol             m/z tolerance window to align MSMS spectra to component list
                       0.002
                       unlimited

    retention_tol      retention time tolerance window to align MSMS spectra to component list
                       5
                       unlimited
 
    results_dir        Directory to write per-cluster plots and table in results_time stamp format
                       None (will write to current directory)
                       can be user defined
 
    out_html           Directory to write interactive dendogram plot too
                       None (will write to results folder)
                       can be user defined
    """

    if quiet:
        printer = lambda *x: None
    else:
        printer = lambda *s: print(*s)

    if results_dir is None:
        dt = datetime.now()
        results_dir = "results_{}".format(dt.strftime("%Y%m%d%H%M%S"))

    printer((
        "Running {name} v{version}\n\n"
        "- input mgf file = {mgf}\n"
        "- input components file = {comp}\n"
        "- neutral = {neu}\n"
        "- cutoff = {cut}\n"
        "- bin_threshold = {bin}\n"
        "- clustering_method = {clu}\n"
        "- output processed file = {pro}\n"
        "- output results directory = {res}\n"
        "- output html dendrogram = {html}\n"
        "- dendrogram figure width = {x}\n"
        "- dendrogram figure height = {y}\n"
        "- scaling = {scal}\n"
        "- filtering = {fil}\n"
        "- eps = {eps}\n"
        "- mz_tolerance = {mz_tol}\n"
        "- retention_tolerance = {retention_tol}\n"
        "\n"
    ).format(
        name=__name__,
        version=__version__,
        mgf=mgf_path,
        comp=components_path,
        neu=neutral,
        cut=cutoff,
        bin=bin_threshold,
        clu=clustering_method,
        res=results_dir,
        pro=pjoin(results_dir, processed),
        html=pjoin(results_dir, out_html),
        x=width,
        y=height,
        scal=scaling,
        fil=filtering,
        eps=eps,
        mz_tol=mz_tol,
        retention_tol=retention_tol,
    ))

    params = [
        ("version", __version__),
        ("input mgf file", mgf_path),
        ("input components file", components_path),
        ("neutral", neutral),
        ("cutoff", cutoff),
        ("bin_threshold", bin_threshold),
        ("clustering_method", clustering_method),
        ("output results directory", results_dir),
        ("output processed file", pjoin(results_dir, processed)),
        ("output html dendrogram", pjoin(results_dir, out_html)),
        ("dendrogram figure width", width),
        ("dendrogram figure height", height),
        ("scaling", scaling),
        ("filtering", filtering),
        ("eps", eps),
    ]

    # Open the trigger data <file>.msg
    with open(mgf_path, 'r') as handle:
        mgf = MGF.parse(handle, scaling=scaling, filtering=filtering, eps=eps)

    # Customised MGF title handler.
    # TODO: This title filter will fail for some mgf title fields.
    for rec in mgf.records:
        rec.title = split_msms_title(rec.title)

    # Open the sample list <file>.csv
    with open(components_path, 'r') as handle:
        components = SampleRecord.parse(handle)

    # Now remove redundancy and print best trigger ion list
    printer("Processing inputs")
    table = remove_redundancy(
        components,
        mgf,
        neutral=neutral,
        mz_tol=mz_tol,
        retention_tol=retention_tol
    )

    printer("Binning and clustering\nThis may take some time...")
    tree = Tree(bin_threshold, clustering_method, cutoff)
    tree.fit(table)

    printer("Writing per-cluster summaries")
    os.makedirs(results_dir, exist_ok=True)
    tree.write_summaries(path=results_dir)

    # Write out an excel file too
    table.drop(columns="mz").drop_duplicates().to_excel(
        pjoin(results_dir, processed), index=False)

    with open(pjoin(results_dir, "params.txt"), "w") as handle:
        params_file = "\n".join(["{}\t{}".format(k, v) for k, v in params])
        handle.write(params_file)

    printer("Writing output html dendrogram")

    _ = tree.plot(
        filename=pjoin(results_dir, out_html),
        width=width,
        height=height
    )

    printer("Finished")
    return tree


def main():
    parser = argparse.ArgumentParser(
        description="Run the BioDendro pipeline."
    )

    parser.add_argument(
        "mgf",
        help="MGF input file.",
    )

    parser.add_argument(
        "components",
        help="Listed components file.",
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
              "the current datetime.\n"
              "WARNING: will overwrite contents of existing directories.")
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
        default=800
    )

    parser.add_argument(
        "-q", "--quiet",
        help="Suppress status notifications written to stdout.",
        action="store_true",
        default=False
    )

    parser.add_argument(
        "-s", "--scaling",
        help="Flag to scale the m/z intensities values",
        action="store_true",
        default=False
    )

    parser.add_argument(
        "-f", "--filtering",
        help="Flag to filter the m/z intensities values",
        action="store_true",
        default=False
    )

    parser.add_argument(
        "-e", "--eps",
        help="Intensity threshold for filtering, set value between 0.0-1.0",
        type=float,
        default=0.6
    )

    parser.add_argument(
        "-mz_tol", "--mz_tol",
        help="Mass tolerance to remove redundancy, default is 0.002",
        type=float,
        default=0.002
    )

    parser.add_argument(
        "-retention_tol", "--retention_tol",
        help=("Retention time tolerance to remove redundancy, "
              "default is 5 secs"),
        type=float,
        default=5
    )

    args = parser.parse_args()

    pipeline(
        mgf_path=args.mgf,
        components_path=args.components,
        **args.__dict__
    )
    return
