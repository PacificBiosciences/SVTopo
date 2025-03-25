import pytest
import matplotlib.pyplot as plt
from svtopovz.svtopovz import (
    create_complex_sv_image,
    get_image_name,
    validate_sv_info
)

class TestImageGeneration:
    @pytest.fixture(autouse=True)
    def setup_method(self):
        plt.switch_backend('Agg')  # Use non-interactive backend for testing

    def test_get_image_name(self, sample_event_info):
        image_name = get_image_name(sample_event_info)
        assert isinstance(image_name, str)
        assert "chr1" in image_name
        assert "-" in image_name  # Should contain coordinate separators

    def test_validate_sv_info_valid(self, sample_sv_info):
        validate_sv_info(sample_sv_info)  # Should not raise exception

    def test_validate_sv_info_invalid(self):
        with pytest.raises(AssertionError):
            validate_sv_info({"invalid": "data"})

    def test_create_complex_sv_image(self, sample_event_info):
        result = create_complex_sv_image(
            sample_event_info,
            gene_records=None,
            bed_records=None,
            max_window_size=1000
        )
        assert isinstance(result, list)
        assert all(isinstance(x, str) for x in result)
        plt.close('all')

    def test_create_complex_sv_image_empty(self):
        result = create_complex_sv_image(
            [],
            gene_records=None,
            bed_records=None,
            max_window_size=1000
        )
        assert result == []
        plt.close('all')

    @pytest.mark.parametrize("window_size", [0.1, 1.0, 10.0])
    def test_different_window_sizes(self, sample_event_info, window_size):
        result = create_complex_sv_image(
            sample_event_info,
            None,
            None,
            window_size
        )
        assert result is not None
        plt.close('all') 