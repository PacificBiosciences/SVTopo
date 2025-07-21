import pytest
from collections import namedtuple
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning, module="pkg_resources")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="mpl_toolkits")
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib.projections")

@pytest.fixture
def sample_event_info():
    return [
        {
            "region": {
                "start_chrom": "chr1",
                "end_chrom": "chr1",
                "start": 1000,
                "end": 2000,
                "variant_ids": []
            },
            "orientation": "+",
            "coverages": {"chr1:1000": 30, "chr1:2000": 30},
            "sample_order_index": 1
        },
        {
            "region": {
                "start_chrom": "chr1",
                "end_chrom": "chr1",
                "start": 2000,
                "end": 3000,
                "variant_ids": []
            },
            "orientation": "+",
            "coverages": {},
            "sample_order_index": 2
        },
        {
            "region": {
                "start_chrom": "chr1",
                "end_chrom": "chr1",
                "start": 3000,
                "end": 4000,
                "variant_ids": []
            },
            "orientation": "+",
            "coverages": {"chr1:3000": 30, "chr1:4000": 30},
            "sample_order_index": 3
        }
    ]

@pytest.fixture
def sample_sv_info():
    return {
        "svtopo_version": "0.2.0",
        "event_graphs": [
            [  # Simple deletion event
                {
                    "region": {"start_chrom": "chr1", "end_chrom": "chr1", "start": 1000, "end": 2000},
                    "orientation": "+",
                    "coverages": {"1000": 30, "2000": 30}
                }
            ]
        ]
    }

@pytest.fixture
def mock_region():
    Region = namedtuple('Region', ['chrom', 'start', 'end', 'size'])
    return Region(chrom="chr1", start=1000, end=2000, size=1000) 