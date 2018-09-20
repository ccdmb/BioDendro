"""
Preprocess contains methods for parsing and manipulating mass spec files.
"""

import os
import re
from bisect import bisect_left
from collections import namedtuple

import pandas as pd


# Named tuple to represent ions and pepmass in MGF
Ion = namedtuple("Ion", ["mz", "intensity"])


class MGF(object):
    """ """

    def __init__(self, records):
        self.records = records
        self.mzs = [r.pepmass.mz for r in records]
        return


    @classmethod
    def parse(cls, handle):
        records = MGFRecord.parse(handle)
        records.sort(key=lambda x: x.pepmass.mz)
        return cls(records)


    def closest(self, mz, retention, mz_tol=0.002, retention_tol=5):
        """ Find the closest trigger match to a mz and retention value.

        Keyword arguments:
        mz --
        retention --
        mz_tol --
        retention_tol --

        Uses:
        self.records
        self.mzs

        NB self.records must be sorted by mzs, and mzs must correspond to the
        records.
        """

        closest = None

        lower_mz = mz - mz_tol
        upper_mz = mz + mz_tol

        lower_retention = retention - retention_tol
        upper_retention = retention + retention_tol

        # Initialise the minimum mass distance
        # Unrealistic number to guarantee match
        min_dist_mz = float("inf")
        # Initialise minimum retention distance
        min_dist_retention = float("inf")

        # Use a binary search to find the first passing mz value
        # Start the loop from there.
        min_bound = bisect_left(self.mzs, lower_mz)

        for record in self.records[min_bound:]:
            if record.pepmass.mz >= upper_mz:
                break
            elif (record.retention <= lower_retention
                    or record.retention >= upper_retention):
                continue

            dist_retention = abs(retention - record.retention)

            # find the closest trigger to the sample
            if dist_retention < min_dist_retention:
                min_dist_retention = dist_retention
                # set best trigger match
                closest = record

        return closest


class MGFRecord(object):
    """ Represents single MGF records intended to be used in a list. """

    def __init__(
            self,
            title,
            retention,
            pepmass,
            charge=None,
            ions=[],
            ):
        self.title = title
        self.retention = retention
        self.pepmass = pepmass
        self.charge = charge
        self.ions = ions
        return


    def __str__(self):
        cls = self.__class__.__name__
        template = ("{}(title='{}', retention={}, pepmass={}, "
                    "charge='{}', ions={})")
        return template.format(cls, self.title, self.retention, self.pepmass,
                               self.charge, self.ions)


    def __repr__(self):
        return str(self)


    @staticmethod
    def _split_kvline(key, string):
        """ Takes a field key and string and removes the key and equals sign """
        return string[len(key) + 1:].strip()


    @classmethod
    def _get_title(cls, string, key="TITLE"):
        """ Process the title field.

        Note there is no standard for storing metadata in here, so specialised
        extraction is risky for compatibility with other peoples datasets.
        """
        return cls._split_kvline(key, string)


    @classmethod
    def _get_pepmass(cls, string, key="PEPMASS"):
        """ Strips key from mass and returns an Ion object. """
        pepmass = cls._split_kvline(key, string)
        return cls._get_ion(pepmass)


    @classmethod
    def _get_retention(cls, string, key="RTINSECONDS"):
        """ Get the retention time in seconds back. """
        rtinseconds = cls._split_kvline(key, string)
        return float(rtinseconds)

    @classmethod
    def _get_charge(cls, string, key="CHARGE"):
        """ Not fully implemented. Just basic support. """
        return cls._split_kvline(key, string)


    @classmethod
    def _get_ion(cls, string):
        """ Convert a string into an Ion named tuple. """
        string = string.split()

        mz = float(string[0].strip())

        # Intensity/peak area field is optional.
        if len(string) > 1:
            intensity = float(string[1].strip())
        else:
            intensity = None

        # NB any other values in here are ignored by specification.
        return Ion(mz, intensity)


    @classmethod
    def _read(cls, lines):
        title = None
        retention = None
        pepmass = None
        charge = None
        ions = []

        for line in lines:
            if line.startswith("TITLE"):
                title = cls._get_title(line)

            elif line.startswith("RTINSECONDS"):
                retention = cls._get_retention(line)

            elif line.startswith("PEPMASS"):
                pepmass = cls._get_pepmass(line)

            elif line.startswith("CHARGE"):
                charge = cls._get_charge(line)

            # Skip lines that we don't explicitly handle.
            elif "=" in line:
                continue

            else:
                # Eventually need to wrap this in try... except
                ion = cls._get_ion(line)
                ions.append(ion)

        return cls(title, retention, pepmass, charge, ions)


    @classmethod
    def parse(cls, handle):
        """ Parses an MGF file into a list of MGF objects.

        keyword arguments:
        handle -- a file like object or list of strings representing the mgf.
        """

        output = []

        in_block = False
        block = []
        for line in handle:
            if line.startswith("END"):
                output.append(cls._read(block))
                block = []
                in_block = False

            elif line.startswith("BEGIN"):
                assert len(block) == 0
                in_block = True

            elif in_block:
                block.append(line.strip())

        return output


def split_msms_title(line):
    """ Split a title line into a useful format. """

    regex = re.compile(r"\\|/")
    sline = line.split(" ")


    # Using a regex to spit on file paths to handle
    # different OS's
    filename = regex.split(sline[1])[-1]
    basename = os.path.splitext(filename)[0]
    return basename


class SampleRecord(object):

    def __init__(self, mz, retention):
        """ A simple class to store 'real samples'. """
        self.mz = mz
        self.retention = retention
        return


    @classmethod
    def _read(cls, line, sep="_"):
        """ Read a line and construct new object. """
        sline = line.strip().split(sep)

        # Get real sample mass
        mz = float(sline[3].lstrip("m/z"))

        # Get real sample retention time in seconds.
        retention = float(sline[4].lstrip('RT')) * 60
        return cls(mz, retention)


    @classmethod
    def parse(cls, handle):
        """ Parse lines in a file-like object and return list of objects. """

        output = []

        for line in handle:
            # Skip any lines with "Components"
            if 'Components' in line:
                continue

            output.append(cls._read(line))

        return output


def remove_redundancy(samples, mgf, mz_tol=0.002, retention_tol=5,
                      neutral=False):
    """ Selects the closest trigger mass to the real sample mass
    Prints the best trigger id and ion list
    """

    output = []

    # Looping through real samples.
    for sample in samples:
        # Find all close triggers in the MGF.
        trigger = mgf.closest(sample.mz, sample.retention, mz_tol,
                              retention_tol)

        if trigger is None:
            print("NUP", sample.__dict__)
            continue

        # Add all of the ion masses
        for ion in trigger.ions:
            if neutral:
                # get neutral loss
                ion_mz = round(ion.mz - trigger.pepmass.mz, 5)
            else:
                ion_mz = ion.mz

            record = (
                "{}_{}_{}".format(trigger.title, trigger.pepmass.mz, trigger.retention),
                ion_mz
                )
            output.append(record)

    # Return the table, sorted by mz
    table = pd.DataFrame(output, columns=['sample', 'mz'])
    table.sort_values(by='mz', inplace=True)
    table.reset_index(drop=True, inplace=True)
    return table
