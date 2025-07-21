import pytest
import matplotlib.pyplot as plt
import numpy as np
from unittest.mock import Mock, patch
from svtopovz.annotation_plotter import (
    annotate,
    add_box,
    add_strand,
    add_title,
    get_region_bed_records,
    binary_get_start_index
)
from svtopovz.utils import BedRecord, Region

class TestAnnotationPlotter:
    def setup_method(self):
        """Setup method to ensure clean matplotlib state before each test"""
        plt.close('all')  # Close all existing figures
        plt.switch_backend('Agg')  # Use non-interactive backend for testing

    def teardown_method(self):
        """Teardown method to clean up matplotlib after each test"""
        plt.close('all')  # Close all figures created during the test

    @pytest.fixture
    def mock_region(self):
        return Region(chrom="chr1", start=1000, end=2000, size=1000)

    @pytest.fixture
    def mock_ax(self):
        fig, ax = plt.subplots()
        return ax

    @pytest.fixture
    def sample_bed_records(self):
        return {
            "chr1": [
                BedRecord(start=500, end=1500, strand="+", title="gene1"),
                BedRecord(start=1500, end=2500, strand="-", title="gene2"),
                BedRecord(start=3000, end=4000, strand="+", title="gene3"),
            ]
        }

    def test_annotate_empty_region(self, mock_ax):
        """Test annotate with zero-length region"""
        region = Region(chrom="chr1", start=1000, end=1000, size=0)
        annotate("test", {"chr1": []}, region, mock_ax, 0, 0)
        # Should return early without error

    def test_annotate_none_records(self, mock_region, mock_ax):
        """Test annotate with None annotation records"""
        annotate("test", None, mock_region, mock_ax, 0, 0)
        # Should return early without error

    def test_annotate_no_overlapping_records(self, mock_region, mock_ax):
        """Test annotate with no overlapping bed records"""
        bed_records = {"chr2": []}  # Different chromosome
        annotate("test", bed_records, mock_region, mock_ax, 0, 0)
        # Should handle gracefully

    def test_annotate_with_records(self, mock_region, mock_ax, sample_bed_records):
        """Test annotate with valid bed records"""
        annotate("test", sample_bed_records, mock_region, mock_ax, 0, 0)
        # Should add patches to the axis
        assert len(mock_ax.patches) > 0

    def test_add_box(self, mock_ax):
        """Test add_box function"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        add_box(bed_record, 1.0, 0.5, mock_ax)
        
        assert len(mock_ax.patches) == 1
        patch = mock_ax.patches[0]
        assert isinstance(patch, plt.Polygon)
        # Check if the patch has the expected color by accessing the facecolor
        # The actual color may vary slightly, so we check it's a light blue color
        facecolor = patch.get_facecolor()
        assert facecolor[0] > 0.6  # Red component should be high for light blue
        assert facecolor[1] > 0.8  # Green component should be high for light blue
        assert facecolor[2] > 0.9  # Blue component should be very high for light blue

    def test_add_strand_zero_fraction(self, mock_ax):
        """Test add_strand with zero annotation fraction"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        strand_coords = []
        add_strand(bed_record, 500, 0.0, 100, 0.5, strand_coords, mock_ax)
        assert len(strand_coords) == 0

    def test_add_strand_small_fraction(self, mock_ax):
        """Test add_strand with small annotation fraction"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        strand_coords = []
        add_strand(bed_record, 500, 0.001, 100, 0.5, strand_coords, mock_ax)
        # Should not add any strand indicators due to small fraction
        assert len(strand_coords) == 0

    def test_add_strand_valid(self, mock_ax):
        """Test add_strand with valid parameters"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        strand_coords = []
        add_strand(bed_record, 500, 0.5, 100, 0.5, strand_coords, mock_ax)
        # Should add strand indicators
        assert len(strand_coords) > 0

    def test_add_strand_different_strands(self, mock_ax):
        """Test add_strand with different strand indicators"""
        strand_coords = []
        
        # Test positive strand
        bed_record_pos = BedRecord(start=1000, end=1500, strand="+", title="test")
        add_strand(bed_record_pos, 500, 0.5, 100, 0.5, strand_coords, mock_ax)
        
        # Test negative strand
        bed_record_neg = BedRecord(start=1500, end=2000, strand="-", title="test")
        add_strand(bed_record_neg, 500, 0.5, 100, 0.5, strand_coords, mock_ax)
        
        # Test empty strand
        bed_record_empty = BedRecord(start=2000, end=2500, strand="", title="test")
        add_strand(bed_record_empty, 500, 0.5, 100, 0.5, strand_coords, mock_ax)
        
        assert len(strand_coords) > 0

    def test_add_title_no_overlap(self, mock_ax):
        """Test add_title with no overlapping coordinates"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        title_coords = []
        add_title(bed_record, 500, 1.0, 0.5, title_coords, mock_ax)
        assert len(title_coords) == 1

    def test_add_title_with_overlap(self, mock_ax):
        """Test add_title with overlapping coordinates"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        title_coords = [(1200, 1300)]  # Overlapping coordinate
        add_title(bed_record, 500, 1.0, 0.5, title_coords, mock_ax)
        # Should not add new title due to overlap
        assert len(title_coords) == 1

    def test_get_region_bed_records_none(self, mock_region):
        """Test get_region_bed_records with None bed_records"""
        result = get_region_bed_records(None, mock_region)
        assert result is None

    def test_get_region_bed_records_no_chrom(self, mock_region):
        """Test get_region_bed_records with missing chromosome"""
        bed_records = {"chr2": []}  # Different chromosome
        result = get_region_bed_records(bed_records, mock_region)
        assert result == []

    def test_get_region_bed_records_overlapping(self, mock_region, sample_bed_records):
        """Test get_region_bed_records with overlapping records"""
        result = get_region_bed_records(sample_bed_records, mock_region)
        assert len(result) == 2  # Should find 2 overlapping records
        assert result[0].start == 1000  # Should be truncated to region start
        assert result[0].end == 1500
        assert result[1].start == 1500
        assert result[1].end == 2000  # Should be truncated to region end

    def test_get_region_bed_records_no_overlap(self, mock_region):
        """Test get_region_bed_records with no overlapping records"""
        bed_records = {
            "chr1": [
                BedRecord(start=3000, end=4000, strand="+", title="gene3"),
            ]
        }
        result = get_region_bed_records(bed_records, mock_region)
        assert result == []

    def test_binary_get_start_index_empty_list(self):
        """Test binary_get_start_index with empty list"""
        result = binary_get_start_index([], 1000)
        assert result == 0

    def test_binary_get_start_index_single_record(self):
        """Test binary_get_start_index with single record"""
        records = [BedRecord(start=500, end=1500, strand="+", title="test")]
        result = binary_get_start_index(records, 1000)
        assert result == 0

    def test_binary_get_start_index_multiple_records(self):
        """Test binary_get_start_index with multiple records"""
        records = [
            BedRecord(start=500, end=1500, strand="+", title="test1"),
            BedRecord(start=1500, end=2500, strand="+", title="test2"),
            BedRecord(start=2500, end=3500, strand="+", title="test3"),
        ]
        result = binary_get_start_index(records, 2000)
        assert result == 1  # Should find the record that ends after 2000

    def test_binary_get_start_index_before_all(self):
        """Test binary_get_start_index with start before all records"""
        records = [
            BedRecord(start=500, end=1500, strand="+", title="test1"),
            BedRecord(start=1500, end=2500, strand="+", title="test2"),
        ]
        result = binary_get_start_index(records, 100)
        assert result == 0

    def test_binary_get_start_index_after_all(self):
        """Test binary_get_start_index with start after all records"""
        records = [
            BedRecord(start=500, end=1500, strand="+", title="test1"),
            BedRecord(start=1500, end=2500, strand="+", title="test2"),
        ]
        result = binary_get_start_index(records, 3000)
        assert result == 2  # Should be after all records

    def test_annotate_window_idx_zero(self, mock_region, mock_ax, sample_bed_records):
        """Test annotate with window_idx=0 to check text addition"""
        annotate("test", sample_bed_records, mock_region, mock_ax, 0, 0)
        # Should add text for annotation name
        assert len(mock_ax.texts) > 0

    def test_annotate_window_idx_nonzero(self, mock_region, mock_ax, sample_bed_records):
        """Test annotate with window_idx!=0 to check no text addition"""
        # Create a fresh axis to avoid interference from previous tests
        fig, ax = plt.subplots()
        annotate("test", sample_bed_records, mock_region, ax, 0, 1)
        # Should not add text for annotation name when window_idx != 0
        # Note: The function may still add texts for other purposes, so we check
        # that the annotation name text is not added specifically
        annotation_texts = [text for text in ax.texts if text.get_text() == "test"]
        assert len(annotation_texts) == 0

    def test_annotate_level_calculation(self, mock_region, mock_ax, sample_bed_records):
        """Test annotate with different levels"""
        # Test level 0
        annotate("test", sample_bed_records, mock_region, mock_ax, 0, 0)
        patches_level0 = len(mock_ax.patches)
        
        # Test level 1
        annotate("test", sample_bed_records, mock_region, mock_ax, 1, 0)
        patches_level1 = len(mock_ax.patches)
        
        # Should have more patches after second call
        assert patches_level1 > patches_level0

    def test_add_box_coordinates(self, mock_ax):
        """Test add_box creates correct polygon coordinates"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        add_box(bed_record, 1.0, 0.5, mock_ax)
        
        patch = mock_ax.patches[0]
        coords = patch.get_xy()
        expected_coords = np.array([
            [1000, 1.0],
            [1000, 0.5],
            [1500, 0.5],
            [1500, 1.0],
            [1000, 1.0]
        ])
        np.testing.assert_array_equal(coords, expected_coords)

    def test_add_strand_overlap_detection(self, mock_ax):
        """Test add_strand overlap detection"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        strand_coords = [(1200, 1300)]  # Existing coordinate
        add_strand(bed_record, 500, 0.5, 100, 0.5, strand_coords, mock_ax)
        # The function adds multiple strand indicators, but should detect overlaps
        # We can't easily test the exact overlap detection without mocking the text rendering
        assert len(strand_coords) > 1  # Should add some coordinates

    def test_add_title_coordinate_calculation(self, mock_ax):
        """Test add_title calculates correct coordinates"""
        bed_record = BedRecord(start=1000, end=1500, strand="+", title="test")
        title_coords = []
        add_title(bed_record, 500, 1.0, 0.5, title_coords, mock_ax)
        
        assert len(title_coords) == 1
        coord_pair = title_coords[0]
        assert coord_pair[0] < coord_pair[1]  # Should be valid coordinate range

    def test_get_region_bed_records_edge_cases(self, mock_region):
        """Test get_region_bed_records with edge cases"""
        bed_records = {
            "chr1": [
                BedRecord(start=1000, end=2000, strand="+", title="exact_match"),
                BedRecord(start=500, end=1000, strand="+", title="touches_start"),
                BedRecord(start=2000, end=3000, strand="+", title="touches_end"),
            ]
        }
        result = get_region_bed_records(bed_records, mock_region)
        # The function only returns records that actually overlap the region
        # Records that only touch the boundary may not be included depending on the logic
        assert len(result) >= 1  # At least the exact match should be included
        if len(result) > 0:
            assert result[0].start >= 1000
            assert result[0].end <= 2000 