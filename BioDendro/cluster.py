"""
Module cluster contains a single class dedicated to hierarchical clustering
of mass spec data.
"""

from os.path import join as pjoin

import pandas as pd
import numpy as np
import scipy as sp
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import pdist

import matplotlib
matplotlib.use("AGG")
from matplotlib import pyplot as plt


from BioDendro.plot import create_dendro


class Tree(object):

    """
    TODO: This should contain example of how to use this.
    """

    def __init__(
            self,
            df,
            bin_threshold=8e-4,
            clustering_method="jaccard",
            sample_col="sample",
            mz_col="mz",
            ):
        """ Constructs a tree object to bin and cluster mass spectra.

        Keyword arguments:
        df -- A pandas dataframe containing samples and mz values.
        bin_threshold -- The threshold to use when dynamically binning mz's.
        clustering_method -- The distance metric used when hierarchically
            clustering samples based on mz bin presence absence.
        sample_col -- The column name in the df to use as the samples.
        mz_col -- The column name in the df to use as the mz values.
        """
        self.df = df.copy()
        self.bin_threshold = bin_threshold
        self.clustering_method = clustering_method
        self.sample_col = sample_col
        self.mz_col = mz_col

        self.bin()
        self.hclust()
        self.select_clusters()
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
        longest_name_len = 0

        # Starting at 1 because we need i - 1. Which will always be 0.
        for i in range(1, len(starts)):
            # Get the slice as an object to avoid rewriting it.
            bin_range = slice(starts[i - 1], starts[i])

            # Get the mz values within this bin.
            arr = column[bin_range]

            # Compute a bin name based on its members' mz values.
            bin_name = cls._bin_name(arr)

            # Update bin name lengths
            if len(bin_name) > longest_name_len:
                longest_bin_len = len(bin_name)

            # assign name to elements in this empty bin.
            bins[bin_range] = bin_name

        # Repeat procedure from above for last cluster.
        # Remembering that the starts mark, the start of the cluster, so we
        # need to go all the way to column end.
        bin_range = slice(starts[i], None)
        arr = column[bin_range]

        bin_name = cls._bin_name(arr)
        if len(bin_name) > longest_name_len:
            longest_bin_len = len(bin_name)
        bins[bin_range] = bin_name

        # Return the bin names. And convert the array type to a C-string with
        # the maximum number of characters required.
        return bins.astype("|S{}".format(longest_bin_len))


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
        present = np.ones(len(df), dtype=np.bool)
        df["present"] = True
        df["bins"] = bins
        return df.pivot_table(index=index_col, columns="bins",
                              values="present", fill_value=False)


    def bin(self, bin_threshold=None):
        """ Get names of the bins and assign to mz rows.

        Keyword arguments:
        bin_threshold -- See __init__. If none, inherits threshold from object.

        Uses:
        self.df
        self.mz_col

        Modifies:
        self.onehot_df -- A one-hot encoded dataframe of samples vs bins.
        """

        if bin_threshold is None:
            bin_threshold = self.bin_threshold

        df = self.df.copy()
        colname = self.mz_col

        bin_starts = self._bin_diffs(df[colname], bin_threshold)
        bins = self._bin_names(df[colname], bin_starts)
        self.onehot_df = self._pivot(df, bins, "sample")
        return


    def hclust(self, clustering_method=None):
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


        df = self.df.copy()
        self.tree = linkage(self.onehot_df, method="complete",
                            metric=clustering_method)
        return


    def select_clusters(self, cutoff=None):
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

        self.clusters = fcluster(self.tree, cutoff, criterion='distance')
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


    def write_summaries(self, path="results"):
        """ Write summary tables and plots to a directory.

        Keyword arguments:
        path -- The directory to write the output to. This directory must
            exist.
        """
        df = self.df
        clusters = self.clusters

        for cluster, subtab in df.groupby(clusters):
            nmembers = subtab.shape[0]

            subtab = subtab.loc[:, subtab.any(axis=0)]
            csv_filename = pjoin(path,
                    "cluster_{}_{}.csv".format(cluster, nmembers))
            subtab.to_csv(csv_filename, sep="\t")

            fig, ax = self.plot_bin_freqs(subtab)
            fig.suptitle("Cluster {} with {} members".format(cluster,
                                                             subtab.shape[0]))
            plt_filename = pjoin(path,
                                 "cluster_{}_{}.png".format(cluster, nmembers))
            fig.savefig(plt_filename)

            # Prevents plotting these plots in interactive mode.
            plt.close()

        filename = pjoin(path, "clusters.csv")
        df["cluster"] = clusters
        df = df[["cluster"] + [c for c in df.columns if c != "cluster"]]
        df.to_csv(filename, sep="\t")
        return


    def iplot(
            self,
            filename=None,
            cutoff=None,
            clustering_method=None,
            width=900,
            height=400
            ):
        """ Plots an interactive tree from these data using plotly.

        Keyword arguments:
        filename -- The path to save the plotly html file.
        cutoff -- See __init__. If None will inherit value from object.
        clustering_method -- As for cutoff.
        width -- Width of the plot in pixels.
        height -- Height of the plot in pixels.

        Uses:
        self.onehot_df
        """

        df = self.onehot_df

        if clustering_method is None:
            clustering_method = self.clustering_method

        if cutoff is None:
            cutoff = self.cutoff


        dendro = create_dendro(
            df,
            labels=df.index,
            distfun=lambda x: pdist(x, metric=clustering_method),
            color_threshold=cutoff
            )

        dendro['layout'].update({
            'width': width,
            'height': height,
            'title': 'BioDendro',
            'xaxis': {'title': 'sample'},
            'yaxis': {'title': 'distance'},
            'hovermode': 'closest'
            })

        dendro['data'].update({'hoverinfo': 'all'})
        plotly.offline.plot(dendro, filename=filename)
        return dendro
