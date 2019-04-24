"""
Module cluster contains a single class dedicated to hierarchical clustering
of mass spec data.
"""

from os.path import join as pjoin

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster

import matplotlib
matplotlib.use("AGG")

from matplotlib import pyplot as plt  # noqa
import plotly  # noqa

from BioDendro.plot import dendrogram  # noqa


class Tree(object):

    """
    TODO: This should contain example of how to use this.
    """

    def __init__(
        self,
        threshold=8e-4,
        clustering_method="jaccard",
        cutoff=0.6,
        sample_col="component",
        mz_col="mz",
    ):
        """ Constructs a tree object to bin and cluster mass spectra.

        Keyword arguments:
        df -- A pandas dataframe containing samples and mz values.
        threshold -- The threshold to use when dynamically binning mz's.
        clustering_method -- The distance metric used when hierarchically
        clustering samples based on mz bin presence absence.
        cutoff -- .
        sample_col -- The column name in the df to use as the samples.
        mz_col -- The column name in the df to use as the mz values.
        """

        self.threshold = threshold
        self.clustering_method = clustering_method
        self.cutoff = cutoff
        self.sample_col = sample_col
        self.mz_col = mz_col

        return

    def fit(self, df):
        """ Bins the data and generates the tree.

        Keyword arguments:
            df -- A pandas dataframe containing samples and mz values.

        Modifies:
            self.df -- Creates a copy of the input data.
            Other elements modified indirectly.
            """

        self.df = df.copy()
        threshold = self.threshold
        clustering_method = self.clustering_method
        cutoff = self.cutoff

        self._bin(threshold)
        self._hclust(clustering_method)
        self.cut_tree(cutoff)
        return

    @staticmethod
    def _bin_name(arr):
        return "{:.4f}_{:.4f}_{:.4f}".format(
            np.around(np.mean(arr), 4),
            np.around(np.min(arr), 4),
            np.around(np.max(arr), 4)
        )

    @staticmethod
    def _bin_starts(arr, threshold):
        """ Find starts of bins in a sorted array.
        Not intended for public use.

        Keyword arguments:
            arr -- A pandas Series object, in sorted order.
            threshold -- See __init__.

        returns:
        np.array
        """

        diffs = arr.diff()
        diffs[0] = threshold + 1
        bin_starts = np.where(diffs >= threshold)[0]
        return bin_starts

    @classmethod
    def _bin_names(cls, column, starts):
        """ Get names of the bins and assign to mz rows.

        Names are derived from bin members, corresponding to the
        mean, min and max, separated by underscores.

        Keyword arguments:
        column -- A list/array/series of mz values
        starts -- A list/array/series of bin start indices.
            e.g. from _bin_starts.

        Returns:
        np.array of strings, corresponding to names of bins.
        Elements in the mz array with the same bin name, belong to the same
        bin.
        """

        # Create an empty array to store python strings
        # These will be the bin names
        bins = np.zeros(len(column), dtype=object)

        # We'll keep track of the name lengths so we can convert to
        # C-strings later
        # longest_name_len = 0

        i = 0
        # Starting at 1 because we need i - 1. Which will always be 0.
        for i in range(1, len(starts)):
            # Get the slice as an object to avoid rewriting it.
            bin_range = slice(starts[i - 1], starts[i])

            # Get the mz values within this bin.
            arr = column[bin_range]

            # Compute a bin name based on its members' mz values.
            bin_name = cls._bin_name(arr)

            # Update bin name lengths
            # if len(bin_name) > longest_name_len:
            #     longest_bin_len = len(bin_name)

            # assign name to elements in this empty bin.
            bins[bin_range] = bin_name

        # Repeat procedure from above for last cluster.
        # Remembering that the starts mark, the start of the cluster, so we
        # need to go all the way to column end.
        bin_range = slice(starts[i], None)
        arr = column[bin_range]

        bin_name = cls._bin_name(arr)
        # if len(bin_name) > longest_name_len:
        #    longest_bin_len = len(bin_name)
        bins[bin_range] = bin_name

        # Return the bin names. And convert the array type to a C-string with
        # the maximum number of characters required.
        # return bins.astype("|S{}".format(longest_bin_len))

        # Uncomment above line if C string performance is necessary.
        # But you would need to decode the strings from binary to get them to
        # print normally.
        return bins

    @staticmethod
    def _pivot(df, bins, index_col):
        """ Construct a one-hot encoded dataframe of samples vs bins.

        Given any dataframe with samples, construct a wide form dataframe
        of boolean values denoting presence absense of ions in a the bin.

        Keyword arguments:
        df -- A pandas dataframe
        bins -- A list/array of ion bins to cluster the rows by. These will
        form the column names. Must be the same length as df, but may
        contain duplicates.
        index_col -- The column to use from df.
        """

        df["present"] = True
        df["bins"] = bins
        return df.pivot_table(index=index_col, columns="bins",
                              values="present", fill_value=False)

    def _bin(self, threshold=None):
        """ Get names of the bins and assign to mz rows.

        Keyword arguments:
        threshold -- See __init__. If none, inherits threshold from object.

        Uses:
        self.df
        self.mz_col
        self.sample_col

        Modifies:
        self.onehot_df -- A one-hot encoded dataframe of samples vs bins.
        """

        if threshold is None:
            threshold = self.threshold

        df = self.df.copy()
        colname = self.mz_col

        bin_starts = self._bin_starts(df[colname], threshold)
        bins = self._bin_names(df[colname], bin_starts)
        self.onehot_df = self._pivot(df, bins, self.sample_col)
        return

    def _hclust(self, clustering_method=None):
        """ Hierarchically cluster the one hot encoded dataframe.

        Keyword arguments:
            clustering_method -- The distance metric used to construct linkage
            with. May be either "jaccard" or "braycurtis". If none inherits
            from object.

        Uses:
            self.onehot_df

        Modifies:
            self.tree -- A scipy linkage array.
            """

        if clustering_method is None:
            clustering_method = self.clustering_method

        self.tree = linkage(self.onehot_df, method="complete",
                            metric=clustering_method)
        return

    def cut_tree(self, cutoff=None):
        """ Selects clusters from the clustered tree based on distance.

        Keyword arguments:
        cutoff -- A float, see __init__ for details.
        If None inherits from object.

        Uses:
        self.cutoff
        self.tree

        Modified:
        self.clusters -- An array corresponding to the clusters.
        Will be the same length as the number of unique samples.
        """

        if cutoff is None:
            cutoff = self.cutoff
        else:
            self.cutoff = cutoff

        self.clusters = fcluster(self.tree, cutoff, criterion='distance')
        self.cluster_map = dict(zip(list(self.onehot_df.index), self.clusters))
        return

    @staticmethod
    def _plot_bin_freqs(df, height=4.5, width_base=1, width_multiplier=0.2):
        """ Plots barchart frequencies of mz bins in a cluster.

        Keyword arguments:
        df -- A pandas dataframe with columns corresponding to mz bins, and
        rows corresponding to samples. Values are boolean or ints [0, 1].
        height -- The plot height in inches.
        width_base -- The basic plot width in inches.
        width_multiplier -- For each bin in the category, add this width in
        inches. Ensures x-axis labels are legible.

        Returns:
        fig -- A matplotlib figure object.
        ax -- A matplotlib axis object, containing the barchart.
        """

        frequencies = df.apply(lambda x: np.mean(x), axis=0)

        width = width_base + width_multiplier * df.shape[1]
        fig, ax = plt.subplots(figsize=(width, height))

        xticks = np.arange(frequencies.shape[0])
        xticklabels = df.columns.values

        ax.bar(xticks, frequencies)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, rotation=90)
        ax.set_ylabel("Frequency")
        ax.set_xlabel("m/z bin")

        # Helps keep the long xticklabels in the output figure.
        fig.tight_layout()
        return fig, ax

    @staticmethod
    def _exclude_false_columns(table):
        """ Filter out columns that are all false. """

        # any(axis=0) at least on sample has True value for each column.
        return table.loc[:, table.any(axis=0)]

    def write_summaries(self, path="results"):
        """ Write summary tables and plots to a directory.

        Keyword arguments:
        path -- The directory to write the output to. This directory must
        exist.
        """

        df = self.onehot_df.copy()
        clusters = self.clusters

        for cluster, subtab in df.groupby(clusters):
            nmembers = subtab.shape[0]

            # Filter out columns that are all false for ease of visualisation.
            subtab = self._exclude_false_columns(subtab)

            csv_filename = pjoin(
                path,
                "cluster_{}_{}.xlsx".format(cluster, nmembers)
            )
            subtab.to_excel(csv_filename)

            fig, ax = self._plot_bin_freqs(subtab)
            fig.suptitle("Cluster {} with {} members".format(
                cluster,
                subtab.shape[0])
            )
            plt_filename = pjoin(
                path,
                "cluster_{}_{}.png".format(cluster, nmembers)
            )
            fig.savefig(plt_filename)

            # Prevents plotting these plots in interactive mode.
            plt.close()

        filename = pjoin(path, "clusters.xlsx")
        df["cluster"] = clusters
        df = df[["cluster"] + [c for c in df.columns if c != "cluster"]]
        df.to_excel(filename)
        return

    def cluster_table(self, cluster=None, sample=None):
        """ Return a table of presence-absence metabolites for a given cluster.
        """

        if cluster is None and sample is None:
            raise ValueError("Either cluster or sample must be set.")
        elif cluster is not None and sample is not None:
            raise ValueError("Please provide a sample or cluster, not both.")
        elif sample is not None:
            try:
                cluster = self.cluster_map[sample]
            except KeyError:
                raise KeyError("The sample you provided isn't in the dataset.")

        subtab = self.onehot_df.loc[self.clusters == cluster]
        return self._exclude_false_columns(subtab)

    def cluster_hist(
        self,
        cluster=None,
        sample=None,
        height=5,
        width_base=1,
        width_multiplier=0.15
    ):
        """ Plot a histogram of metabolite frequencies in the data. """

        subtab = self.cluster_table(cluster, sample)
        return self._plot_bin_freqs(subtab, height,
                                    width_base, width_multiplier)

    def plot(
        self,
        filename=None,
        width=900,
        height=800,
        fontsize=12,
        auto_open=True,
    ):
        """ Plots an interactive tree from these data using plotly.

        Keyword arguments:
        filename -- The path to save the plotly html file. If None, don't write
        the plot.
        width -- Width of the plot in pixels.
        height -- Height of the plot in pixels.
        fontsize -- Fontsize used for xlabels. Note this doesn't currently
            change the actual fontsize, it is used to scale the bottom margin
            so that the labels don't get cut off.

        Uses:
        self.onehot_df
        self.clusters
        self.tree

        Returns:
        A dictionary suitable to be give to plotly.
        """

        title = (
            "Component clusters. method = {}, cutoff = {}, threshold = {}"
        ).format(
            self.clustering_method,
            self.cutoff,
            self.threshold,
        )

        dendro = dendrogram(
            self,
            width=width,
            height=height,
            title=title,
            xlabel="Components",
            ylabel="Distance",
            margin_scalar=fontsize,
        )

        if filename is not None:
            plotly.offline.plot(
                dendro,
                filename=filename,
            )

        return dendro
