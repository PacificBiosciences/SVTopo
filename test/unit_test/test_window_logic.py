import pytest
from svtopovz.svtopovz import (
    choose_window_coordinates,
    get_positions_in_window,
    get_windows_by_chromosome,
    get_windows
)

class TestWindowLogic:
    def test_get_positions_in_window(self, sample_event_info):
        positions = get_positions_in_window(sample_event_info)
        assert "chr1" in positions
        assert sorted(positions["chr1"]) == [1000, 2000, 3000, 4000]

    def test_get_positions_multi_chrom(self):
        multi_chrom_event = [
            {
                "region": {
                    "start_chrom": "chr1",
                    "end_chrom": "chr2",
                    "start": 1000,
                    "end": 2000,
                }
            }
        ]
        positions = get_positions_in_window(multi_chrom_event)
        assert "chr1" in positions and "chr2" in positions
        assert positions["chr1"] == [1000]
        assert positions["chr2"] == [2000]

    def test_get_windows_by_chromosome(self):
        positions = {"chr1": [1000, 2000, 10000]}  # Gap of 8000 between 2000 and 10000
        windows = get_windows_by_chromosome(positions, max_size=5000)
        assert len(windows["chr1"]) == 2  # Should split into two windows due to gap

    def test_choose_window_coordinates_zero_max_size(self, sample_event_info):
        windows = choose_window_coordinates(sample_event_info, 0)
        assert len(windows) > 0  # Should use BIG_NUMBER instead of 0

    @pytest.mark.parametrize("max_size,expected_windows", [
        (1000, 3),  # Small max size -> more windows
        (10000, 1), # Large max size -> fewer windows
    ])
    def test_choose_window_coordinates_different_sizes(self, sample_event_info, max_size, expected_windows):
        windows = choose_window_coordinates(sample_event_info, max_size)
        assert len(windows) == expected_windows 