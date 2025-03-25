import pytest
from svtopovz.svtopovz import is_filtered_out

class TestFilters:
    def test_simple_deletion(self):
        event_info = [
            {"orientation": "+", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
            {"orientation": "+", "coverages": {}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
            {"orientation": "+", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
        ]
        assert is_filtered_out(0, event_info, False) is True  # Should filter simple deletions
        assert is_filtered_out(0, event_info, True) is False  # Should keep when including simple

    def test_empty_event(self):
        assert is_filtered_out(0, [], False) is True
        assert is_filtered_out(0, [], True) is True

    @pytest.mark.parametrize("event_info,expected", [
        # Simple deletion
        ([
            {"orientation": "+", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
            {"orientation": "+", "coverages": {}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
            {"orientation": "+", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
        ], True),
        # Simple duplication
        ([
            {"orientation": "+", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
            {"orientation": "-", "coverages": {}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
            {"orientation": "+", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
        ], True),
        # Complex event
        ([
            {"orientation": "+", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr1"}},
            {"orientation": "-", "coverages": {"1": 30}, "region": {"start_chrom": "chr1", "end_chrom": "chr2"}},
        ], False),
    ])
    def test_event_types(self, event_info, expected):
        assert is_filtered_out(0, event_info, False) is expected 