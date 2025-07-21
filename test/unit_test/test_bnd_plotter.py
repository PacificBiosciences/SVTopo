#!/usr/bin/env python
"""
Unit tests for bnd_plotter module
"""

import pytest
import matplotlib.pyplot as plt
import numpy as np
from unittest.mock import Mock, patch, MagicMock
import sys
import os

from svtopovz.bnd_plotter import (
    plot_spanned_block,
    plot_coverages,
    plot_rectangle,
    plot_arrowhead,
    scale_coverages,
    is_inverted,
    generate_midpoints,
    get_unspanned_connections,
    get_unspanned_connection_candidates,
    refine_unspanned_connections,
    plot_unspanned_block,
    reset_coordinates,
    plot_unspanned_blocks,
    plot_spanned_blocks,
    plot_breaks
)
from svtopovz.utils import Region, MAX_ALLOWED_COVERAGE, MIN_CLUSTER_PRECISION


class TestBndPlotter:
    """Test class for bnd_plotter functions"""

    def setup_method(self):
        """Setup method to ensure clean matplotlib state before each test"""
        plt.close('all')  # Close all existing figures

    def teardown_method(self):
        """Teardown method to clean up matplotlib after each test"""
        plt.close('all')  # Close all figures created during the test

    @pytest.fixture
    def mock_ax(self):
        """Create a mock matplotlib axis"""
        return Mock()

    @pytest.fixture
    def sample_window(self):
        """Create a sample window for testing"""
        return Region("chr1", 1000, 2000, 1000)

    @pytest.fixture
    def sample_coverages(self):
        """Create sample coverage data"""
        return {
            "chr1:1000": 5,
            "chr1:1500": 10,
            "chr1:2000": 3
        }

    @pytest.fixture
    def sample_block(self):
        """Create a sample block for testing"""
        return {
            "region": {
                "start_chrom": "chr1",
                "end_chrom": "chr1",
                "start": 1000,
                "end": 2000
            },
            "sample_order_index": 0,
            "orientation": "+",
            "coverages": {
                "chr1:1000": 5,
                "chr1:1500": 10,
                "chr1:2000": 3
            }
        }

    @pytest.fixture
    def sample_spanned_block_ends(self):
        """Create sample spanned block ends data"""
        return {
            "starts": {
                (1000, 0): {"region": {"start": 1000, "end": 1500}, "orientation": "+"},
                (2000, 1): {"region": {"start": 2000, "end": 2500}, "orientation": "-"}
            },
            "ends": {
                (1500, 0): {"region": {"start": 1000, "end": 1500}, "orientation": "+"},
                (2500, 1): {"region": {"start": 2000, "end": 2500}, "orientation": "-"}
            }
        }

    def test_scale_coverages_below_max(self, sample_coverages):
        """Test scale_coverages when coverage is below maximum"""
        result = scale_coverages(sample_coverages.copy())
        assert result == sample_coverages
        assert max(result.values()) <= MAX_ALLOWED_COVERAGE

    def test_scale_coverages_above_max(self):
        """Test scale_coverages when coverage is above maximum"""
        high_coverages = {
            "chr1:1000": 50,
            "chr1:1500": 100,
            "chr1:2000": 75
        }
        result = scale_coverages(high_coverages.copy())
        assert max(result.values()) <= MAX_ALLOWED_COVERAGE
        # Check that relative ratios are preserved
        original_ratios = [high_coverages[k] / max(high_coverages.values()) for k in sorted(high_coverages.keys())]
        result_ratios = [result[k] / max(result.values()) for k in sorted(result.keys())]
        np.testing.assert_array_almost_equal(original_ratios, result_ratios, decimal=5)

    def test_scale_coverages_empty(self):
        """Test scale_coverages with empty coverage dict"""
        with pytest.raises(ValueError):
            scale_coverages({})

    def test_is_inverted_normal(self, sample_block):
        """Test is_inverted with normal (non-inverted) block"""
        assert not is_inverted(sample_block)

    def test_is_inverted_inverted(self):
        """Test is_inverted with inverted block"""
        inverted_block = {
            "region": {
                "start": 2000,
                "end": 1000
            }
        }
        assert is_inverted(inverted_block)

    def test_generate_midpoints(self):
        """Test generate_midpoints function"""
        pointa = (0, 0)
        pointb = (10, 10)
        midpoints = generate_midpoints(pointa, pointb)
        
        assert len(midpoints) == 10
        assert midpoints[0] == pointa
        # The function generates 10 points, but the last one is at 9/10 of the way
        assert midpoints[-1] == (9.0, 9.0)
        
        # Check that midpoints are evenly distributed
        for i in range(1, 9):
            expected_x = i
            expected_y = i
            assert midpoints[i] == (expected_x, expected_y)

    def test_generate_midpoints_negative_coords(self):
        """Test generate_midpoints with negative coordinates"""
        pointa = (-10, -5)
        pointb = (10, 5)
        midpoints = generate_midpoints(pointa, pointb)
        
        assert len(midpoints) == 10
        assert midpoints[0] == pointa
        # The function generates 10 points, but the last one is at 9/10 of the way
        assert midpoints[-1] == (8.0, 4.0)

    def test_get_unspanned_connection_candidates(self, sample_block, sample_spanned_block_ends):
        """Test get_unspanned_connection_candidates"""
        prev_connections, next_connections = get_unspanned_connection_candidates(
            sample_block, sample_spanned_block_ends
        )
        
        assert isinstance(prev_connections, list)
        assert isinstance(next_connections, list)
        # Should find connections within MIN_CLUSTER_PRECISION
        for connection, sample_idx in prev_connections + next_connections:
            assert isinstance(connection, dict)
            assert isinstance(sample_idx, int)

    def test_get_unspanned_connection_candidates_no_matches(self, sample_block):
        """Test get_unspanned_connection_candidates with no matching connections"""
        empty_spanned_ends = {"starts": {}, "ends": {}}
        prev_connections, next_connections = get_unspanned_connection_candidates(
            sample_block, empty_spanned_ends
        )
        
        assert prev_connections == []
        assert next_connections == []

    def test_refine_unspanned_connections_single(self):
        """Test refine_unspanned_connections with single connection"""
        connections = [({"region": {"start": 1000}}, 0)]
        result_connection, result_idx = refine_unspanned_connections(connections, 1)
        
        assert result_connection == connections[0][0]
        assert result_idx == 0

    def test_refine_unspanned_connections_multiple(self):
        """Test refine_unspanned_connections with multiple connections"""
        connections = [
            ({"region": {"start": 1000}}, 2),
            ({"region": {"start": 1500}}, 1),
            ({"region": {"start": 2000}}, 3)
        ]
        result_connection, result_idx = refine_unspanned_connections(connections, 1)
        
        # Should select the one closest to current_sample_index (1)
        assert result_idx == 1

    def test_refine_unspanned_connections_empty(self):
        """Test refine_unspanned_connections with empty connections"""
        result_connection, result_idx = refine_unspanned_connections([], 1)
        
        assert result_connection is None
        assert result_idx is None

    def test_get_unspanned_connections(self, sample_block, sample_spanned_block_ends):
        """Test get_unspanned_connections"""
        result = get_unspanned_connections(sample_block, sample_spanned_block_ends)
        
        assert isinstance(result, dict)
        assert "prev_block" in result
        assert "next_block" in result
        assert "prev_sample_idx" in result
        assert "next_sample_idx" in result

    def test_reset_coordinates(self, sample_window, sample_block):
        """Test reset_coordinates function"""
        result = reset_coordinates(sample_window, sample_block.copy())
        
        assert result["region"]["start"] >= sample_window.start
        assert result["region"]["end"] <= sample_window.end

    def test_reset_coordinates_outside_window(self, sample_window):
        """Test reset_coordinates with block outside window"""
        outside_block = {
            "region": {
                "start_chrom": "chr1",
                "end_chrom": "chr1",
                "start": 500,
                "end": 2500
            }
        }
        result = reset_coordinates(sample_window, outside_block.copy())
        
        # The function doesn't clip coordinates, it only handles chromosome mismatches
        # For same chromosome, coordinates should remain unchanged
        assert result["region"]["start"] == 500
        assert result["region"]["end"] == 2500

    @patch('svtopovz.bnd_plotter.get_unspanned_connections')
    @patch('svtopovz.bnd_plotter.reset_coordinates')
    def test_plot_unspanned_block(self, mock_reset, mock_get_connections, mock_ax, sample_window, sample_block, sample_spanned_block_ends):
        """Test plot_unspanned_block function"""
        mock_get_connections.return_value = {
            "prev_block": None,
            "next_block": None,
            "prev_sample_idx": None,
            "next_sample_idx": None
        }
        mock_reset.return_value = sample_block
        
        # The function returns None when no connections are found
        result = plot_unspanned_block(
            sample_block, 0.5, sample_window, sample_spanned_block_ends, 0, mock_ax
        )
        
        assert result is None
        mock_get_connections.assert_called_once()
        mock_reset.assert_called_once()

    def test_plot_coverages(self, mock_ax, sample_coverages):
        """Test plot_coverages function"""
        sorted_coverages = sorted(sample_coverages.keys())
        
        plot_coverages(mock_ax, sample_coverages, sorted_coverages, 1, "chr1", 0.5)
        
        # Should call hlines for coverage plotting
        assert mock_ax.hlines.called

    def test_plot_coverages_empty(self, mock_ax):
        """Test plot_coverages with empty coverages"""
        with pytest.raises(IndexError):
            plot_coverages(mock_ax, {}, [], 1, "chr1", 0.5)

    def test_plot_rectangle(self, mock_ax):
        """Test plot_rectangle function"""
        plot_rectangle(mock_ax, 1, 1000, 2000, 0.5)
        
        # Should add a polygon patch
        assert mock_ax.add_patch.called

    def test_plot_rectangle_multiple_matches(self, mock_ax):
        """Test plot_rectangle with multiple matched order count"""
        plot_rectangle(mock_ax, 3, 1000, 2000, 0.5)
        
        # Should add a polygon patch with different color
        assert mock_ax.add_patch.called

    def test_plot_arrowhead_forward(self, mock_ax, sample_window):
        """Test plot_arrowhead with forward orientation"""
        start_coords, end_coords = plot_arrowhead(
            mock_ax, sample_window.size, 1, 1000, 2000, "+", 0, 0.5
        )
        
        assert isinstance(start_coords, list)
        assert isinstance(end_coords, list)
        assert mock_ax.add_patch.called

    def test_plot_arrowhead_reverse(self, mock_ax, sample_window):
        """Test plot_arrowhead with reverse orientation"""
        start_coords, end_coords = plot_arrowhead(
            mock_ax, sample_window.size, 1, 1000, 2000, "-", 0, 0.5
        )
        
        assert isinstance(start_coords, list)
        assert isinstance(end_coords, list)
        assert mock_ax.add_patch.called

    def test_plot_arrowhead_no_orientation(self, mock_ax, sample_window):
        """Test plot_arrowhead with no orientation"""
        start_coords, end_coords = plot_arrowhead(
            mock_ax, sample_window.size, 1, 1000, 2000, "", 0, 0.5
        )
        
        assert isinstance(start_coords, list)
        assert isinstance(end_coords, list)
        # Should not add patch when no orientation
        assert not mock_ax.add_patch.called

    @patch('svtopovz.bnd_plotter.plot_unspanned_block')
    def test_plot_unspanned_blocks(self, mock_plot_unspanned, mock_ax, sample_window):
        """Test plot_unspanned_blocks function"""
        region_info = {
            0: [{
                "region": {"start": 1000, "end": 2000},
                "sample_order_index": 0,
                "coverages": {}
            }]
        }
        spanned_block_ends = {"starts": {}, "ends": {}}
        
        plot_unspanned_blocks(
            region_info, 0, 1, spanned_block_ends, sample_window, 0, mock_ax
        )
        
        # Should call plot_unspanned_block for each unspanned block
        assert mock_plot_unspanned.called

    @patch('svtopovz.bnd_plotter.plot_spanned_block')
    def test_plot_spanned_blocks(self, mock_plot_spanned, mock_ax, sample_window):
        """Test plot_spanned_blocks function"""
        region_info = {
            0: [{
                "region": {"start": 1000, "end": 2000},
                "sample_order_index": 0,
                "coverages": {"chr1:1000": 5},
                "orientation": "+"
            }]
        }
        
        # Mock the return value to avoid unpacking error
        mock_plot_spanned.return_value = ([], [])
        
        plot_spanned_blocks(region_info, 0, 1, sample_window, 0, mock_ax)
        
        # Should call plot_spanned_block for each spanned block
        assert mock_plot_spanned.called

    @patch('svtopovz.bnd_plotter.plot_spanned_blocks')
    @patch('svtopovz.bnd_plotter.plot_unspanned_blocks')
    def test_plot_breaks(self, mock_plot_unspanned, mock_plot_spanned, mock_ax, sample_window):
        """Test plot_breaks function"""
        region_info = {
            0: [{
                "region": {"start": 1000, "end": 2000},
                "sample_order_index": 0,
                "coverages": {"chr1:1000": 5},
                "orientation": "+"
            }]
        }
        
        plot_breaks(region_info, 0, sample_window, 1, 1, 0, mock_ax)
        
        # Should call both plotting functions
        assert mock_plot_spanned.called or mock_plot_unspanned.called

    def test_plot_spanned_block_basic(self, mock_ax, sample_window, sample_coverages):
        """Test plot_spanned_block with basic inputs"""
        start_coords, end_coords = plot_spanned_block(
            sample_coverages, 0.5, sample_window, "+", 1, 0, 0, mock_ax
        )
        
        assert isinstance(start_coords, list)
        assert isinstance(end_coords, list)

    def test_plot_spanned_block_empty_coverages(self, mock_ax, sample_window):
        """Test plot_spanned_block with empty coverages"""
        with pytest.raises(ValueError):
            plot_spanned_block(
                {}, 0.5, sample_window, "+", 1, 0, 0, mock_ax
            )

    def test_plot_spanned_block_wrong_chromosome(self, mock_ax, sample_coverages):
        """Test plot_spanned_block with wrong chromosome"""
        wrong_chrom_window = Region("chr2", 1000, 2000, 1000)
        start_coords, end_coords = plot_spanned_block(
            sample_coverages, 0.5, wrong_chrom_window, "+", 1, 0, 0, mock_ax
        )
        
        assert start_coords == []
        assert end_coords == []

    def test_plot_spanned_block_outside_window(self, mock_ax, sample_window):
        """Test plot_spanned_block with coordinates outside window"""
        outside_coverages = {
            "chr1:500": 5,
            "chr1:2500": 10
        }
        start_coords, end_coords = plot_spanned_block(
            outside_coverages, 0.5, sample_window, "+", 1, 0, 0, mock_ax
        )
        
        # The function returns coordinates even for outside window
        assert isinstance(start_coords, list)
        assert isinstance(end_coords, list) 