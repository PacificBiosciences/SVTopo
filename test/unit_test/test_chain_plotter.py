#!/usr/bin/env python
"""
Unit tests for chain_plotter module
"""

import pytest
import matplotlib.pyplot as plt
import numpy as np
from unittest.mock import Mock, patch, MagicMock
import sys
import os

from svtopovz.chain_plotter import (
    plot_chain_representation,
    generate_reference_blocks,
    generate_sample_ordered_blocks,
    get_max_tied_block_count,
    add_chain_ends,
    flip_inverted_blocks,
    get_spanned_blocks,
    convert_to_letters,
    plot_reference_blocks_proportionally,
    plot_reference_blocks_same_size,
    plot_sample_paths,
    plot_block,
    get_polygon,
    plot_reference_block_sizes
)
from svtopovz.utils import Region, MIN_CLUSTER_PRECISION


class TestChainPlotter:
    """Test class for chain_plotter functions"""

    def setup_method(self):
        """Setup method to ensure clean matplotlib state before each test"""
        plt.close('all')  # Close all existing figures

    def teardown_method(self):
        """Teardown method to clean up matplotlib after each test"""
        plt.close('all')  # Close all figures created during the test

    @pytest.fixture
    def mock_ax(self):
        """Create a mock matplotlib axis with get_ylim returning a tuple"""
        ax = Mock()
        ax.get_ylim.return_value = (0, 10)
        ax.text = Mock()
        ax.add_patch = Mock()
        return ax

    @pytest.fixture
    def sample_windows(self):
        """Create sample windows for testing"""
        return [
            Region("chr1", 1000, 2000, 1000),
            Region("chr1", 2000, 3000, 1000)
        ]

    @pytest.fixture
    def sample_region_info(self):
        """Create sample region info for testing"""
        return {
            0: [{
                "region": {
                    "start_chrom": "chr1",
                    "end_chrom": "chr1",
                    "start": 1000,
                    "end": 2000
                },
                "sample_order_index": 0,
                "coverages": {"chr1:1000": 5, "chr1:1500": 10},
                "orientation": "+"
            }],
            1: [{
                "region": {
                    "start_chrom": "chr1",
                    "end_chrom": "chr1",
                    "start": 2000,
                    "end": 3000
                },
                "sample_order_index": 1,
                "coverages": {"chr1:2000": 3, "chr1:2500": 7},
                "orientation": "-"
            }]
        }

    @pytest.fixture
    def sample_ref_blocks_list(self):
        """Create sample reference blocks as a list for plotting functions"""
        return [
            ("chr1", 1000, "chr1", 2000),
            ("chr1", 2000, "chr1", 3000)
        ]

    @pytest.fixture
    def sample_ref_blocks_dict(self):
        """Create sample reference blocks as a dict for path-finding functions"""
        return {
            ("chr1", 1000, "chr1", 2000): 0,
            ("chr1", 2000, "chr1", 3000): 1
        }

    @pytest.fixture
    def sample_ref_blocks_by_start(self):
        """Create sample reference blocks by start"""
        return {
            ("chr1", 1000): ("chr1", 1000, "chr1", 2000),
            ("chr1", 2000): ("chr1", 2000, "chr1", 3000)
        }

    @pytest.fixture
    def sample_ref_blocks_by_end(self):
        """Create sample reference blocks by end"""
        return {
            ("chr1", 2000): ("chr1", 1000, "chr1", 2000),
            ("chr1", 3000): ("chr1", 2000, "chr1", 3000)
        }

    def test_get_max_tied_block_count_single(self, sample_region_info):
        """Test get_max_tied_block_count with single blocks"""
        result = get_max_tied_block_count(sample_region_info)
        assert result == 1

    def test_get_max_tied_block_count_multiple(self):
        """Test get_max_tied_block_count with multiple blocks"""
        region_info = {
            0: [
                {"coverages": {"chr1:1000": 5}},
                {"coverages": {"chr1:1500": 10}},
                {"coverages": {"chr1:2000": 3}}
            ],
            1: [
                {"coverages": {"chr1:2500": 7}}
            ]
        }
        result = get_max_tied_block_count(region_info)
        assert result == 3

    def test_get_max_tied_block_count_empty(self):
        """Test get_max_tied_block_count with empty region info"""
        result = get_max_tied_block_count({})
        assert result == 0

    def test_get_max_tied_block_count_no_coverages(self):
        """Test get_max_tied_block_count with no coverages"""
        region_info = {
            0: [{"coverages": {}}],
            1: [{"coverages": {}}]
        }
        result = get_max_tied_block_count(region_info)
        assert result == 0

    def test_generate_reference_blocks_basic(self, sample_region_info, sample_windows):
        """Test generate_reference_blocks with basic inputs"""
        ref_blocks, ref_blocks_by_start, ref_blocks_by_end = generate_reference_blocks(
            sample_region_info, sample_windows
        )
        
        assert isinstance(ref_blocks, dict)
        assert isinstance(ref_blocks_by_start, dict)
        assert isinstance(ref_blocks_by_end, dict)
        assert len(ref_blocks) > 0

    def test_generate_reference_blocks_empty(self, sample_windows):
        """Test generate_reference_blocks with empty region info"""
        ref_blocks, ref_blocks_by_start, ref_blocks_by_end = generate_reference_blocks(
            {}, sample_windows
        )
        
        assert isinstance(ref_blocks, dict)
        assert isinstance(ref_blocks_by_start, dict)
        assert isinstance(ref_blocks_by_end, dict)

    def test_generate_reference_blocks_multi_chrom(self, sample_windows):
        """Test generate_reference_blocks with multi-chromosome data"""
        region_info = {
            0: [{
                "region": {
                    "start_chrom": "chr1",
                    "end_chrom": "chr2",
                    "start": 1000,
                    "end": 2000
                }
            }]
        }
        ref_blocks, ref_blocks_by_start, ref_blocks_by_end = generate_reference_blocks(
            region_info, sample_windows
        )
        
        # Should not include blocks that span between chromosomes
        for block in ref_blocks:
            assert block[0] == block[2]  # Same chromosome

    def test_convert_to_letters_basic(self):
        """Test convert_to_letters with basic inputs"""
        assert convert_to_letters(0, []) == "A"
        assert convert_to_letters(1, []) == "B"
        assert convert_to_letters(25, []) == "Z"
        assert convert_to_letters(26, []) == "AA"

    def test_convert_to_letters_with_skipped(self):
        """Test convert_to_letters with skipped blocks"""
        assert convert_to_letters(1, [0]) == "A"  # B becomes A when 0 is skipped
        assert convert_to_letters(2, [0, 1]) == "A"  # C becomes A when 0,1 are skipped

    def test_convert_to_letters_negative(self):
        """Test convert_to_letters with negative input"""
        with pytest.raises(SystemExit):
            convert_to_letters(-1, [])

    def test_flip_inverted_blocks_none(self):
        """Test flip_inverted_blocks with None input"""
        result = flip_inverted_blocks(None)
        assert result is None

    def test_flip_inverted_blocks_no_inverted(self):
        """Test flip_inverted_blocks with no inverted blocks"""
        sample_paths = [
            [([0, 1], "+")],
            [([2, 3], "+")]
        ]
        result = flip_inverted_blocks(sample_paths)
        assert result == sample_paths

    def test_flip_inverted_blocks_with_inverted(self):
        """Test flip_inverted_blocks with inverted blocks"""
        sample_paths = [
            [([0, 1], "-")],
            [([2, 3], "+")]
        ]
        result = flip_inverted_blocks(sample_paths)
        
        # First block should be flipped
        assert result[0][0][0] == [1, 0]  # Reversed indices
        assert result[0][0][1] == "-"  # Orientation unchanged
        # Second block should be unchanged
        assert result[1][0][0] == [2, 3]
        assert result[1][0][1] == "+"

    def test_add_chain_ends_basic(self, sample_ref_blocks_list):
        """Test add_chain_ends with basic inputs"""
        sample_paths = [
            [([0, 1], "+")],
            [([1], "-")]
        ]
        # Create ref_blocks_by_idx mapping integer indices to tuples
        ref_blocks_by_idx = {0: sample_ref_blocks_list[0], 1: sample_ref_blocks_list[1]}
        result = add_chain_ends(sample_paths, 1, ref_blocks_by_idx)
        assert isinstance(result, list)
        assert len(result) == len(sample_paths)

    def test_add_chain_ends_empty_paths(self, sample_ref_blocks_list):
        sample_paths = []
        ref_blocks_by_idx = {0: sample_ref_blocks_list[0], 1: sample_ref_blocks_list[1]}
        result = add_chain_ends(sample_paths, 1, ref_blocks_by_idx)
        assert result == []

    def test_add_chain_ends_empty_blocks(self, sample_ref_blocks_list):
        sample_paths = [
            [([], "+")],
            [([], "-")]
        ]
        ref_blocks_by_idx = {0: sample_ref_blocks_list[0], 1: sample_ref_blocks_list[1]}
        result = add_chain_ends(sample_paths, 1, ref_blocks_by_idx)
        assert result is None

    def test_get_spanned_blocks_direct_unspanned(self):
        """Test get_spanned_blocks with direct unspanned connection"""
        block = {
            "region": {"start_chrom": "chr1"},
            "orientation": "+"
        }
        prev_info = {
            "is_direct_unspanned": True,
            "orientation": "+",
            "chrom": "chr1"
        }
        
        result = get_spanned_blocks(block, prev_info, 0, 1)
        
        assert isinstance(result, list)
        assert len(result) == 1
        assert result[0][1] == "+"  # Orientation

    def test_get_spanned_blocks_different_chrom(self):
        """Test get_spanned_blocks with different chromosome"""
        block = {
            "region": {"start_chrom": "chr2"},
            "orientation": "+"
        }
        prev_info = {
            "is_direct_unspanned": False,
            "orientation": "+",
            "chrom": "chr1"
        }
        
        result = get_spanned_blocks(block, prev_info, 0, 1)
        
        assert isinstance(result, list)
        assert len(result) == 1

    def test_get_spanned_blocks_unknown_orientation(self):
        """Test get_spanned_blocks with unknown orientation"""
        block = {
            "region": {"start_chrom": "chr1"},
            "orientation": ""
        }
        prev_info = {
            "is_direct_unspanned": True,
            "orientation": "+",
            "chrom": "chr1"
        }
        
        result = get_spanned_blocks(block, prev_info, 0, 1)
        
        assert isinstance(result, list)
        assert len(result) == 1
        assert result[0][1] == "+"  # Should inherit from prev_info

    def test_get_spanned_blocks_with_prev_info(self):
        """Test get_spanned_blocks with previous block info"""
        block = {
            "region": {"start_chrom": "chr1"},
            "orientation": "+"
        }
        prev_info = {
            "is_direct_unspanned": False,
            "orientation": "+",
            "start_index": 0,
            "end_index": 1,
            "chrom": "chr1"
        }
        
        result = get_spanned_blocks(block, prev_info, 2, 3)
        
        assert isinstance(result, list)
        assert len(result) == 1

    def test_generate_sample_ordered_blocks_basic(self, sample_ref_blocks_dict, sample_ref_blocks_by_start, sample_ref_blocks_by_end, sample_region_info):
        """Test generate_sample_ordered_blocks with basic inputs"""
        result = generate_sample_ordered_blocks(
            sample_ref_blocks_dict, sample_ref_blocks_by_start, sample_ref_blocks_by_end, sample_region_info
        )
        
        assert isinstance(result, list)
        assert len(result) > 0

    def test_generate_sample_ordered_blocks_too_many_tied(self, sample_ref_blocks_dict, sample_ref_blocks_by_start, sample_ref_blocks_by_end):
        """Test generate_sample_ordered_blocks with too many tied blocks"""
        region_info = {
            0: [
                {"coverages": {"chr1:1000": 5}},
                {"coverages": {"chr1:1500": 10}},
                {"coverages": {"chr1:2000": 3}},
                {"coverages": {"chr1:2500": 7}},
                {"coverages": {"chr1:3000": 2}}
            ]
        }
        
        result = generate_sample_ordered_blocks(
            sample_ref_blocks_dict, sample_ref_blocks_by_start, sample_ref_blocks_by_end, region_info
        )
        
        # Should return None when too many tied blocks
        assert result is None

    def test_generate_sample_ordered_blocks_no_matches(self, sample_ref_blocks_dict, sample_ref_blocks_by_start, sample_ref_blocks_by_end):
        """Test generate_sample_ordered_blocks with no matching blocks"""
        region_info = {
            0: [{
                "region": {
                    "start_chrom": "chr1",
                    "end_chrom": "chr1",
                    "start": 9999,  # Non-matching coordinates
                    "end": 9999
                },
                "coverages": {"chr1:9999": 5},
                "orientation": "+"
            }]
        }
        
        result = generate_sample_ordered_blocks(
            sample_ref_blocks_dict, sample_ref_blocks_by_start, sample_ref_blocks_by_end, region_info
        )
        
        # Should return None when no blocks match
        assert result is None

    @patch('svtopovz.chain_plotter.generate_reference_blocks')
    @patch('svtopovz.chain_plotter.generate_sample_ordered_blocks')
    @patch('svtopovz.chain_plotter.plot_reference_blocks_proportionally')
    @patch('svtopovz.chain_plotter.plot_reference_blocks_same_size')
    @patch('svtopovz.chain_plotter.plot_sample_paths')
    def test_plot_chain_representation(self, mock_plot_paths, mock_plot_same_size, mock_plot_prop, mock_gen_sample, mock_gen_ref, sample_region_info, sample_windows):
        """Test plot_chain_representation function"""
        mock_gen_ref.return_value = ({}, {}, {})
        mock_gen_sample.return_value = [[([0, 1], "+")]]
        mock_plot_paths.return_value = True
        
        result = plot_chain_representation(sample_region_info, sample_windows, {"breaks": [Mock()], "chains": Mock()})
        
        assert result is True
        mock_gen_ref.assert_called_once()
        mock_gen_sample.assert_called_once()
        mock_plot_prop.assert_called_once()
        mock_plot_same_size.assert_called_once()
        mock_plot_paths.assert_called_once()

    @patch('svtopovz.chain_plotter.generate_reference_blocks')
    @patch('svtopovz.chain_plotter.generate_sample_ordered_blocks')
    def test_plot_chain_representation_no_paths(self, mock_gen_sample, mock_gen_ref, sample_region_info, sample_windows):
        """Test plot_chain_representation when no paths are found"""
        mock_gen_ref.return_value = ({}, {}, {})
        mock_gen_sample.return_value = None
        
        result = plot_chain_representation(sample_region_info, sample_windows, {"breaks": [Mock()], "chains": Mock()})
        
        assert result is False

    def test_plot_reference_blocks_proportionally(self, sample_ref_blocks_dict, sample_windows, mock_ax):
        """Test plot_reference_blocks_proportionally function"""
        colors = ["red", "blue"]
        skipped_ref_blocks = []
        axes = [mock_ax, mock_ax]
        
        # Create ref_blocks_by_idx mapping integer indices to tuples
        ref_blocks_by_idx = {v: k for k, v in sample_ref_blocks_dict.items()}
        
        plot_reference_blocks_proportionally(
            ref_blocks_by_idx, colors, sample_windows, axes, skipped_ref_blocks
        )
        
        # Should add text and patches to axes
        for ax in axes:
            assert ax.text.called or ax.add_patch.called

    def test_plot_reference_blocks_same_size(self, sample_ref_blocks_list, mock_ax):
        """Test plot_reference_blocks_same_size function"""
        colors = ["red", "blue"]
        skipped_ref_blocks = []
        region_size = 2000
        
        plot_reference_blocks_same_size(
            sample_ref_blocks_list, colors, region_size, mock_ax, skipped_ref_blocks
        )
        
        # Should add patches to axis
        assert mock_ax.add_patch.called

    def test_plot_sample_paths(self, sample_ref_blocks_list, mock_ax):
        """Test plot_sample_paths function"""
        sample_paths = [[([0, 1], "+")], [([1], "-")]]
        colors = ["red", "blue"]
        skipped_ref_blocks = []
        region_size = 2000
        
        result = plot_sample_paths(
            sample_paths, sample_ref_blocks_list, colors, region_size, mock_ax, skipped_ref_blocks
        )
        
        # Should add patches to axis
        assert mock_ax.add_patch.called

    def test_plot_block(self, sample_ref_blocks_list, mock_ax):
        """Test plot_block function"""
        block_numbers = [0, 1]
        colors = ["red", "blue"]
        skipped_ref_blocks = []
        region_size = 2000
        
        plot_block(
            mock_ax, block_numbers, 0, sample_ref_blocks_list, skipped_ref_blocks,
            colors, 0.5, region_size, "+", 1000, 0
        )
        
        # Should add patches to axis
        assert mock_ax.add_patch.called

    def test_get_polygon_forward(self):
        """Test get_polygon with forward orientation"""
        result = get_polygon(1000, 0, 1000, 2000, 0.5, "+")
        
        assert isinstance(result, list)
        assert len(result) > 0
        # Should have arrowhead pointing right for forward orientation
        assert result[2][0] > result[1][0]  # Arrowhead point is right of base

    def test_get_polygon_reverse(self):
        """Test get_polygon with reverse orientation"""
        result = get_polygon(1000, 0, 1000, 2000, 0.5, "-")
        
        assert isinstance(result, list)
        assert len(result) > 0
        # Should have arrowhead pointing left for reverse orientation
        assert result[2][0] < result[1][0]  # Arrowhead point is left of base

    def test_get_polygon_zero_length(self):
        """Test get_polygon with zero length"""
        result = get_polygon(0, 0, 0, 2000, 0.5, "+")
        
        assert isinstance(result, list)
        assert len(result) > 0

    @patch('matplotlib.pyplot.legend')
    def test_plot_reference_block_sizes(self, mock_legend):
        """Test plot_reference_block_sizes function"""
        block_lens = [("A", 1000, "red"), ("B", 2000, "blue")]
        plot_reference_block_sizes(block_lens)
        # The function should call legend once
        assert mock_legend.called 