#!/usr/bin/env python
"""
Comprehensive unit tests for svtopovz module
"""

import pytest
import matplotlib.pyplot as plt
from unittest.mock import Mock, patch

from svtopovz.svtopovz import (
    choose_window_coordinates,
    get_positions_in_window,
    get_windows_by_chromosome,
    get_windows,
    create_complex_sv_image,
    get_image_name,
    validate_sv_info,
    is_filtered_out,
    make_plot,
    make_plots,
    svtopovz,
    MakePlotArgs,
)
from svtopovz.utils import Region


class TestSvtopovz:
    """Comprehensive test class for svtopovz functions"""

    def setup_method(self):
        """Setup method to ensure clean matplotlib state before each test"""
        plt.close("all")  # Close all existing figures
        plt.switch_backend("Agg")  # Use non-interactive backend for testing

    def teardown_method(self):
        """Teardown method to clean up matplotlib after each test"""
        plt.close("all")  # Close all figures created during the test

    # ============================================================================
    # Window Logic Tests (from test_window_logic.py)
    # ============================================================================

    def test_get_positions_in_window(self, sample_event_info):
        """Test get_positions_in_window with sample event info"""
        positions = get_positions_in_window(sample_event_info)
        assert "chr1" in positions
        assert sorted(positions["chr1"]) == [1000, 2000, 3000, 4000]

    def test_get_positions_multi_chrom(self):
        """Test get_positions_in_window with multi-chromosome event"""
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

    def test_get_positions_in_window_empty(self):
        """Test get_positions_in_window with empty event info"""
        positions = get_positions_in_window([])
        assert positions == {}

    def test_get_windows_by_chromosome_single_position(self):
        """Test get_windows_by_chromosome with single position"""
        positions = {"chr1": [1000]}
        windows = get_windows_by_chromosome(positions, max_size=5000)
        assert len(windows["chr1"]) == 1
        window = windows["chr1"][0]
        assert window.chrom == "chr1"
        assert window.start == 1000
        assert window.end == 1000

    def test_get_windows_by_chromosome_multiple_positions_small_gap(self):
        """Test get_windows_by_chromosome with multiple positions that fit in one window"""
        positions = {"chr1": [1000, 2000, 3000]}  # Gap of 1000 between positions
        windows = get_windows_by_chromosome(positions, max_size=5000)
        assert len(windows["chr1"]) == 1  # Should fit in one window
        window = windows["chr1"][0]
        assert window.start == 1000
        assert window.end == 3000

    def test_get_windows_by_chromosome_multiple_positions_large_gap(self):
        """Test get_windows_by_chromosome with multiple positions that need multiple windows"""
        positions = {"chr1": [1000, 2000, 10000]}  # Gap of 8000 between 2000 and 10000
        windows = get_windows_by_chromosome(positions, max_size=5000)
        assert len(windows["chr1"]) == 2  # Should split into two windows due to gap

    def test_get_windows_by_chromosome_multiple_positions_two_window(self):
        """Test get_windows_by_chromosome with multiple positions that need multiple windows"""
        positions = {
            "chr9": [6713695, 6727606, 6727607, 6746138, 30688951, 30690313, 30690512]
        }
        windows = get_windows_by_chromosome(positions, max_size=500_000)
        expected_windows = {
            "chr9": [
                Region(chrom="chr9", start=6713695, end=6746138, size=32443),
                Region(chrom="chr9", start=30688951, end=30690512, size=1561),
            ]
        }
        assert windows == expected_windows

    def test_get_windows_by_chromosome_empty_positions(self):
        """Test get_windows_by_chromosome with empty positions"""
        positions = {"chr1": []}
        windows = get_windows_by_chromosome(positions, max_size=5000)
        assert len(windows) == 0

    def test_choose_window_coordinates_zero_max_size(self, sample_event_info):
        """Test choose_window_coordinates with zero max size"""
        windows = choose_window_coordinates(sample_event_info, 0)
        expected = [Region(chrom="chr1", start=0, end=5000, size=5000)]
        assert expected == windows  # Should use BIG_NUMBER instead of 0

    @pytest.mark.parametrize(
        "max_size,expected_windows",
        [
            (1100, 2),  # Small max size -> more windows
            (10000, 1),  # Large max size -> fewer windows
        ],
    )
    def test_choose_window_coordinates_different_sizes(
        self, sample_event_info, max_size, expected_windows
    ):
        """Test choose_window_coordinates with different max sizes"""
        windows = choose_window_coordinates(sample_event_info, max_size)
        assert len(windows) == expected_windows

    def test_choose_window_coordinates_empty_event(self):
        """Test choose_window_coordinates with empty event info"""
        windows = choose_window_coordinates([], 1000)
        assert windows == []

    def test_get_windows_with_padding(self):
        """Test get_windows adds padding to first and last windows"""
        windows_by_chrom = {
            "chr1": [Region("chr1", 1000, 2000, 1000), Region("chr1", 3000, 4000, 1000)]
        }
        positions = {"chr1": [1000, 2000, 3000, 4000]}
        windows = get_windows(windows_by_chrom, positions)
        assert len(windows) == 2
        # First window should have padding subtracted from start
        assert windows[0].start < 1000
        # Last window should have padding added to end
        assert windows[1].end > 4000

    def test_get_windows_empty_windows(self):
        """Test get_windows with empty windows but positions"""
        windows_by_chrom = {}
        positions = {"chr1": [1000, 2000, 3000]}
        windows = get_windows(windows_by_chrom, positions)
        assert len(windows) == 1
        # Should create a single window with padding
        assert windows[0].start < 1000
        assert windows[0].end > 3000

    # ============================================================================
    # Image Generation Tests (from test_image_generation.py)
    # ============================================================================

    def test_get_image_name(self, sample_event_info):
        """Test get_image_name with sample event info"""
        image_name = get_image_name(sample_event_info)
        assert isinstance(image_name, str)
        assert "chr1" in image_name
        assert "-" in image_name  # Should contain coordinate separators

    def test_get_image_name_empty_event(self):
        """Test get_image_name with empty event info"""
        image_name = get_image_name([])
        assert image_name == ""

    def test_get_image_name_single_block(self):
        """Test get_image_name with single block"""
        event_info = [
            {
                "coverages": {"chr1:1000": 30},
                "region": {
                    "start_chrom": "chr1",
                    "end_chrom": "chr1",
                    "start": 1000,
                    "end": 2000,
                },
            }
        ]
        image_name = get_image_name(event_info)
        assert "chr1-1000-2000" in image_name

    def test_get_image_name_multi_chrom(self):
        """Test get_image_name with multi-chromosome event"""
        event_info = [
            {
                "coverages": {"chr1:1000": 30},
                "region": {
                    "start_chrom": "chr1",
                    "end_chrom": "chr1",
                    "start": 1000,
                    "end": 2000,
                },
            },
            {
                "coverages": {"chr2:3000": 30},
                "region": {
                    "start_chrom": "chr2",
                    "end_chrom": "chr2",
                    "start": 3000,
                    "end": 4000,
                },
            },
        ]
        image_name = get_image_name(event_info)
        assert "chr1-1000-2000" in image_name
        assert "chr2-3000-4000" in image_name
        assert "__" in image_name  # Should separate chromosomes

    def test_validate_sv_info_valid(self, sample_sv_info):
        """Test validate_sv_info with valid SV info"""
        validate_sv_info(sample_sv_info)  # Should not raise exception

    def test_validate_sv_info_invalid(self):
        """Test validate_sv_info with invalid SV info"""
        with pytest.raises(AssertionError):
            validate_sv_info({"invalid": "data"})

    def test_validate_sv_info_missing_version(self):
        """Test validate_sv_info with missing version"""
        sv_info = {"event_graphs": []}
        with pytest.raises(AssertionError):
            validate_sv_info(sv_info)

    def test_validate_sv_info_missing_event_graphs(self):
        """Test validate_sv_info with missing event_graphs"""
        sv_info = {"svtopo_version": "0.3.0"}
        with pytest.raises(AssertionError):
            validate_sv_info(sv_info)

    def test_create_complex_sv_image(self, sample_event_info):
        """Test create_complex_sv_image with sample event info"""
        result = create_complex_sv_image(
            sample_event_info, gene_records={}, bed_records={}, max_window_size=1000
        )
        assert isinstance(result, list)
        assert all(isinstance(x, str) for x in result)

    def test_create_complex_sv_image_empty(self):
        """Test create_complex_sv_image with empty event info"""
        result = create_complex_sv_image(
            [], gene_records={}, bed_records={}, max_window_size=1000
        )
        assert result == []

    @pytest.mark.parametrize("window_size", [0.1, 1.0, 10.0])
    def test_different_window_sizes(self, sample_event_info, window_size):
        """Test create_complex_sv_image with different window sizes"""
        result = create_complex_sv_image(sample_event_info, {}, {}, window_size)
        assert result is not None

    def test_create_complex_sv_image_with_annotations(self, sample_event_info):
        """Test create_complex_sv_image with gene and bed annotations"""
        gene_records = {"genes.gtf": []}
        bed_records = {"annotations.bed": []}
        result = create_complex_sv_image(
            sample_event_info,
            gene_records=gene_records,
            bed_records=bed_records,
            max_window_size=1000,
        )
        assert isinstance(result, list)

    # ============================================================================
    # Filter Tests (from test_filters.py)
    # ============================================================================

    def test_simple_deletion(self):
        """Test is_filtered_out with simple deletion event"""
        event_info = [
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
            {
                "orientation": "+",
                "coverages": {},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
        ]
        assert (
            is_filtered_out(0, event_info, False) is True
        )  # Should filter simple deletions
        assert (
            is_filtered_out(0, event_info, True) is False
        )  # Should keep when including simple

    def test_empty_event(self):
        """Test is_filtered_out with empty event"""
        assert is_filtered_out(0, [], False) is True
        assert is_filtered_out(0, [], True) is True

    @pytest.mark.parametrize(
        "event_info,expected",
        [
            # Simple deletion
            (
                [
                    {
                        "orientation": "+",
                        "coverages": {"1": 30},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
                    },
                    {
                        "orientation": "+",
                        "coverages": {},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
                    },
                    {
                        "orientation": "+",
                        "coverages": {"1": 30},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
                    },
                ],
                True,
            ),
            # Simple duplication
            (
                [
                    {
                        "orientation": "+",
                        "coverages": {"1": 30},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
                    },
                    {
                        "orientation": "-",
                        "coverages": {},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
                    },
                    {
                        "orientation": "+",
                        "coverages": {"1": 30},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
                    },
                ],
                True,
            ),
            # Complex event
            (
                [
                    {
                        "orientation": "+",
                        "coverages": {"1": 30},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
                    },
                    {
                        "orientation": "-",
                        "coverages": {"1": 30},
                        "region": {"start_chrom": "chr1", "end_chrom": "chr2"},
                    },
                ],
                False,
            ),
        ],
    )
    def test_event_types(self, event_info, expected):
        """Test is_filtered_out with different event types"""
        assert is_filtered_out(0, event_info, False) is expected

    def test_simple_duplication(self):
        """Test is_filtered_out with simple duplication event"""
        event_info = [
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
            {
                "orientation": "-",
                "coverages": {},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
        ]
        assert (
            is_filtered_out(0, event_info, False) is True
        )  # Should filter simple duplications
        assert (
            is_filtered_out(0, event_info, True) is False
        )  # Should keep when including simple

    def test_nonreciprocal_translocation(self):
        """Test is_filtered_out with nonreciprocal translocation"""
        event_info = [
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
            {
                "orientation": "+",
                "coverages": {},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr2", "end_chrom": "chr2"},
            },
        ]
        assert (
            is_filtered_out(0, event_info, False) is True
        )  # Should filter nonreciprocal translocations
        assert (
            is_filtered_out(0, event_info, True) is False
        )  # Should keep when including simple

    def test_single_ended_bnd(self):
        """Test is_filtered_out with single-ended BND"""
        event_info = [
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
        ]
        assert (
            is_filtered_out(0, event_info, False) is True
        )  # Should filter single-ended BNDs
        assert (
            is_filtered_out(0, event_info, True) is False
        )  # Should keep when including simple

    def test_complex_event_not_filtered(self):
        """Test is_filtered_out with complex event that should not be filtered"""
        event_info = [
            {
                "orientation": "+",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr1"},
            },
            {
                "orientation": "-",
                "coverages": {"1": 30},
                "region": {"start_chrom": "chr1", "end_chrom": "chr2"},
            },
        ]
        assert (
            is_filtered_out(0, event_info, False) is False
        )  # Should not filter complex events

    # ============================================================================
    # Additional Tests for Untested Functions
    # ============================================================================

    @patch("svtopovz.svtopovz.create_complex_sv_image")
    @patch("matplotlib.pyplot.savefig")
    @patch("matplotlib.pyplot.close")
    def test_make_plot(
        self, mock_close, mock_savefig, mock_create_image, sample_event_info
    ):
        """Test make_plot function"""
        mock_create_image.return_value = ["chr1_1000_2000"]

        plot_args = MakePlotArgs(
            event_info=sample_event_info,
            gene_annotation_records={},
            bed_records={},
            max_gap_size_mb=1.0,
            prefix="test",
            svtopo_dir="/tmp",
            image_type="png",
        )

        with patch("os.path.join", return_value="/tmp/images/test_chr1-1000-2000.png"):
            make_plot(plot_args)

        mock_create_image.assert_called_once()
        mock_savefig.assert_called_once()
        mock_close.assert_called_once()

    @patch("svtopovz.svtopovz.make_plot")
    @patch("svtopovz.svtopovz.unpack_json")
    @patch("svtopovz.svtopovz.validate_sv_info")
    @patch("svtopovz.svtopovz.is_filtered_out")
    def test_make_plots(
        self, mock_filtered, mock_validate, mock_unpack, mock_make_plot, sample_sv_info
    ):
        """Test make_plots function"""
        mock_unpack.return_value = (sample_sv_info, "test_prefix")
        mock_filtered.return_value = False  # Don't filter out any events

        jsons = ["test.json"]
        gene_records = {}
        bed_records = {}

        make_plots(
            jsons=jsons,
            gene_annotation_records=gene_records,
            bed_records=bed_records,
            include_simple_breakpoints=True,
            max_gap_size_mb=1.0,
            svtopo_dir="/tmp",
            image_type="png",
        )

        mock_unpack.assert_called_once()
        mock_validate.assert_called_once()
        mock_make_plot.assert_called()

    @patch("svtopovz.svtopovz.make_plots")
    @patch("svtopovz.svtopovz.build_review_page")
    @patch("svtopovz.svtopovz.unpack_annotation_records")
    @patch("svtopovz.svtopovz.glob")
    def test_svtopovz(
        self, mock_glob, mock_unpack_annotations, mock_build_page, mock_make_plots
    ):
        """Test svtopovz main function"""
        mock_glob.side_effect = [["test1.json"], ["test2.json.gz"]]
        mock_unpack_annotations.return_value = {}

        # Create a mock args object
        args = Mock()
        args.svtopo_dir = "/tmp"
        args.genes = "genes.gtf"
        args.annotation_bed = ["bed1.bed", "bed2.bed"]
        args.include_simple_breakpoints = True
        args.max_gap_size_mb = 1.0
        args.image_type = "png"

        svtopovz(args)

        mock_glob.assert_called()
        mock_unpack_annotations.assert_called()
        mock_make_plots.assert_called_once()
        mock_build_page.assert_called_once()

    def test_make_plot_args_namedtuple(self):
        """Test MakePlotArgs namedtuple creation"""
        args = MakePlotArgs(
            event_info=[],
            gene_annotation_records={},
            bed_records={},
            max_gap_size_mb=1.0,
            prefix="test",
            svtopo_dir="/tmp",
            image_type="png",
        )
        assert args.event_info == []
        assert args.prefix == "test"
        assert args.image_type == "png"
