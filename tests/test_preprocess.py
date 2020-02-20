import pytest

from BioDendro.preprocess import split_msms_title
from BioDendro.preprocess import MGFRecord
from BioDendro.preprocess import Ion
from BioDendro.preprocess import SampleRecord


# Test MGFRecord methods
# NB _get_charge is untested since unused

def test_MGFRecord__split_kvline():
    sample = "TITLE=This; is=some text"
    expected = "This; is=some text"

    actual = MGFRecord._split_kvline("TITLE", sample)
    assert actual == expected
    return


def test_MGFRecord__get_title():
    """ Note this is intentionally the same as split kvline """
    sample = "TITLE=This; is=some text"
    expected = "This; is=some text"

    actual = MGFRecord._get_title(sample)
    assert actual == expected


@pytest.mark.parametrize("sample,expected", [
    ("PEPMASS=6.66 5.11", Ion(mz=6.66, intensity=5.11)),
    ("PEPMASS=2.22", Ion(mz=2.22, intensity=None)),
    ])
def test_MGFRecord__get_pepmass(sample, expected):
    actual = MGFRecord._get_pepmass(sample)
    assert actual == expected
    return


@pytest.mark.parametrize("sample,expected", [
    ("RTINSECONDS=6.66", 6.66),
    ("RTINSECONDS=2", 2.0),
    ])
def test_MGFRecord__get_retention(sample, expected):
    actual = MGFRecord._get_retention(sample)
    assert actual == expected
    return


@pytest.mark.parametrize("sample,expected", [
    ("6.66 5.11", Ion(mz=6.66, intensity=5.11)),
    ("2.22", Ion(mz=2.22, intensity=None)),
    ("2.22 3 40", Ion(mz=2.22, intensity=3.0)),
    ])
def test_MGFRecord__get_ion(sample, expected):
    actual = MGFRecord._get_ion(sample)
    assert actual == expected
    return


@pytest.mark.parametrize("sample,expected", [
    (
        ['TITLE=File: "\\Mac\\Home\\Desktop\\TEMP\\QE_2017_001814.raw"; SpectrumID: "2"; scans: "2"',
         'PEPMASS=102.03323 33553618.00000',
         'CHARGE=2+',
         'RTINSECONDS=0',
         'SCANS=2',
         '56.96532 908954',
         '58.06568 37828.7',
         '60.04488 38341.2'],
        MGFRecord(
            title='File: "\\Mac\\Home\\Desktop\\TEMP\\QE_2017_001814.raw"; SpectrumID: "2"; scans: "2"',
            retention=0,
            pepmass=Ion(102.03323, 33553618.0),
            charge="2+",
            ions=[
                Ion(56.96532, 908954),
                Ion(58.06568, 37828.7),
                Ion(60.04488, 38341.2),
                ]
            )
    ),
    (
        ["TITLE=scan=986 profile data",
         "RTINSECONDS=297.916",
         "PEPMASS=758.571517944336 12066.720502853394",
         "CHARGE=1+",
         "37.05708507 1.0",
         "38.08264425 1.0",
         "42.1891818 1.0"],
        MGFRecord(
            title='scan=986 profile data',
            retention=297.916,
            pepmass=Ion(758.571517944336, 12066.720502853394),
            charge="1+",
            ions=[
                Ion(37.05708507, 1.0),
                Ion(38.08264425, 1.0),
                Ion(42.1891818, 1.0),
                ]
            )
    ),
    ])
def test_MGFRecord__read(sample, expected):
    actual = MGFRecord._read(sample)

    assert actual.title == expected.title
    assert actual.retention == expected.retention
    assert actual.pepmass == expected.pepmass
    assert actual.charge == expected.charge

    for act_ion, exp_ion in zip(actual.ions, expected.ions):
        assert act_ion == exp_ion
    return


@pytest.mark.parametrize("sample,expected", [
    (
        [
            'MASS=Monoisotopic',
            'BEGIN IONS',
            'TITLE=File: "\\Mac\\Home\\Desktop\\TEMP\\QE_2017_001814.raw"; SpectrumID: "2"; scans: "2"',
            'PEPMASS=102.03323 33553618.00000',
            'CHARGE=2+',
            'RTINSECONDS=0',
            'SCANS=2',
            '56.96532 908954',
            '58.06568 37828.7',
            '60.04488 38341.2',
            'END IONS',
            'BEGIN IONS',
            "TITLE=scan=986 profile data",
            "RTINSECONDS=297.916",
            "PEPMASS=758.571517944336 12066.720502853394",
            "CHARGE=1+",
            "37.05708507 1.0",
            "38.08264425 1.0",
            "42.1891818 1.0",
            'END IONS',
        ],
        [
            MGFRecord(
                title='File: "\\Mac\\Home\\Desktop\\TEMP\\QE_2017_001814.raw"; SpectrumID: "2"; scans: "2"',
                retention=0,
                pepmass=Ion(102.03323, 33553618.0),
                charge="2+",
                ions=[
                    Ion(56.96532, 908954),
                    Ion(58.06568, 37828.7),
                    Ion(60.04488, 38341.2),
                    ]
                ),
            MGFRecord(
                title='scan=986 profile data',
                retention=297.916,
                pepmass=Ion(758.571517944336, 12066.720502853394),
                charge="1+",
                ions=[
                    Ion(37.05708507, 1.0),
                    Ion(38.08264425, 1.0),
                    Ion(42.1891818, 1.0),
                    ]
            )
        ]
    ),
    ])
def test_MGFRecord_parse(sample, expected):
    actual = MGFRecord.parse(sample)

    assert len(actual) == len(expected)

    for act_rec, exp_rec in zip(actual, expected):
        assert act_rec.title == exp_rec.title
        assert act_rec.retention == exp_rec.retention
        assert act_rec.pepmass == exp_rec.pepmass
        assert act_rec.charge == exp_rec.charge

        for act_ion, exp_ion in zip(act_rec.ions, exp_rec.ions):
            assert act_ion == exp_ion
    return


# Test MGF methods

def test_MGF_closest():
    return


# Test SampleRecord methods

@pytest.mark.parametrize("sample,expected", [
    ("Chinese Spring 1_001829_4250_m/z129.1274_RT5.4654\n", {"mz": 129.1274, "retention": 5.4654 * 60}),
    ("blk_001809_0035_m/z141.0508_RT0.5557", {"mz": 141.0508, "retention": 0.5557 * 60}),
    ("Cobra 2_001852_2090_m/z194.1174_RT4.0473", {"mz": 194.1174, "retention": 4.0473 * 60}),
    ("Espada 1_001827_1146_m/z235.1439_RT1.5430", {"mz": 235.1439, "retention": 1.5430* 60})
    ])
def test_SampleRecord__read(sample, expected):
    actual = SampleRecord._read(sample, sep="_")

    for key, exp in expected.items():
        assert getattr(actual, key) == exp
    return


@pytest.mark.parametrize("sample,expected", [
    (
        ["Chinese Spring 1_001829_4250_m/z129.1274_RT5.4654",
         "blk_001809_0035_m/z141.0508_RT0.5557",
         "Cobra 2_001852_2090_m/z194.1174_RT4.0473",
         "Espada 1_001827_1146_m/z235.1439_RT1.5430"],
        [{"mz": 129.1274, "retention": 5.4654 * 60},
         {"mz": 141.0508, "retention": 0.5557 * 60},
         {"mz": 194.1174, "retention": 4.0473 * 60},
         {"mz": 235.1439, "retention": 1.5430* 60}]
    )
    ])
def test_SampleRecord_parse(sample, expected):
    actual = SampleRecord.parse(sample)

    assert len(actual) == len(expected)
    for act, exp in zip(actual, expected):
        for key, exp_val in exp.items():
            assert getattr(act, key) == exp_val
    return


# Test standalone methods

@pytest.mark.parametrize("sample,expected", [
    ('File: "\\Mac\\Home\\Desktop\\TEMP\\QE_2017_001814.raw"; SpectrumID: "2"; scans: "2"', "QE_2017_001814"),
    (r'File: "/home/user/test_name.raw"; SpectrumID: "3"', "test_name")
    ])
def test_split_msms_title(sample, expected):
    """ Split title should return a short version of the title. """
    actual = split_msms_title(sample)

    assert actual == expected
    return


