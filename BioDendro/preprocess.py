
#Author: P Moolhuijzen
#Date  : 26 February 2018
#Inputs: <input1>.mgf and <input1>.csv text files
#Ouputs: <msms_nonredundant_list>.txt

import os
import re
import sys
import argparse
import logging
from collections import defaultdict
from collections import namedtuple

import xlsxwriter
import csv
import operator
import pandas as pd


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

logger.debug("Loaded module: preprocess")

# Named tuple to represent ions and pepmass in MGF
Ion = namedtuple("Ion", ["mz", "intensity"])

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

        # Possibly should assert that some fields not none.
        logger.debug("Created record: {}".format(self))
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
        logger.debug("Parsing MGF file")

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


def get_csv_record(handle, mgf, mz_tol=0.002, retention_tol=5):
    """ Collects all the csv list (real samples list) matches to the trigger
    data within a mass of 0.002 and retention time of 5 secs return a
    dictionary key real sample id and value of triggers (id and ions list)
    """

    output = {}

    for line in handle:
        # Skip any lines with "Components"
        if 'Components' in line:
            continue

        # real sample mass
        line = line.rstrip('\n')
        sline = line.split('_')
        mz = float(sline[3].lstrip('m/z'))

        # real sample retention time
        retention = float(sline[4].lstrip('RT'))
        retention = int(retention * 60)

        upper_mz = mz + 0.002
        lower_mz = mz - mz_tol

        upper_retention = retention + retention_tol
        lower_retention = retention - retention_tol

        comp = []

        # Using t_ to denote "trigger" since retention and mz are common terms.
        for trigger in mgf:
            if ((lower_mz < trigger.pepmass.mz < upper_mz) and
                    (lower_retention < trigger.retention < upper_retention)):
                comp.append(trigger)

        if len(comp) > 0:
            output[(mz, retention)] = comp

    return output


def remove_redundancy(ndic, neutral=False):
    """ Selects the closest trigger mass to the real sample mass
    Prints the best trigger id and ion list
    """

    output = []

    # Looping through real samples.
    for (mz, retention), triggers in ndic.items():
        # Initialise the minimum mass distance
        d_mz = 10 # Unrealistic number to guarantee match

        # Initialise minimum retention distance
        d_retention = 1000

        for trigger in triggers:
            this_d_mz = abs(mz - trigger.pepmass.mz)
            this_d_retention = abs(retention - trigger.retention)

            # find the closest trigger to the sample
            if this_d_mz < d_mz and this_d_retention < d_retention:
                # new distance of trigger mass to sample
                d_mz = this_d_mz
                # new distance of trigger retention time to sample
                d_retention = this_d_retention
                # set best trigger match
                best = trigger

        for ion in best.ions:
            if neutral:
                # get negative loss
                ion_mz = round(ion.mz - trigger.pepmass.mz, 5)
            else:
                ion_mz = ion.mz
            record = (
                "{}_{}_{}".format(best.title, best.pepmass.mz, best.retention),
                ion_mz
                )

            output.append(record)


    # Return the table, sorted by mz
    table = pd.DataFrame(output, columns=['sample', 'mz'])
    table.sort_values(by='mz', inplace=True)
    return table
