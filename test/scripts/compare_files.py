from pprint import pprint
from skimage.metrics import structural_similarity as ssim
from skimage import io
import argparse
import sys
import json

MATCHING_SSI = 0.98


class Colors:
    RED = "\033[91m"
    GREEN = "\033[92m"
    RESET = "\033[0m"


def setup():
    parser = argparse.ArgumentParser(
        description="Compare two images to determine if essentially identical",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-r",
        "--reference-file-path",
        help="the path to a reference file: json, image, or test",
    )
    parser.add_argument(
        "-t",
        "--test-file-path",
        help="the path to a test file: json, image, or text",
    )
    args = parser.parse_args()
    return args


def compare_plots(reference_image_path, test_image_path):
    """
    Compare two images and report an error if they're dissimilar
    """
    test_image = io.imread(test_image_path)
    reference_image = io.imread(reference_image_path)
    # window size has to be equal to or less than the smallest side of the input images
    win_size = int(min(test_image.shape[:2] + reference_image.shape[:2]))
    # must be odd
    if win_size % 2 == 0:
        win_size -= 1
    similarity_index, _ = ssim(
        test_image, reference_image, full=True, win_size=win_size, channel_axis=-1
    )
    print("SSI:", similarity_index)
    if similarity_index < MATCHING_SSI:
        print(
            "{}Test image {} does not match reference image {}: SSI={}{}".format(
                Colors.RED,
                test_image_path,
                reference_image_path,
                similarity_index,
                Colors.RESET,
            ),
            file=sys.stderr,
        )
        sys.exit(-1)
    else:
        print(
            "{}Test image {} matches reference image {}: SSI={}{}".format(
                Colors.GREEN,
                test_image_path,
                reference_image_path,
                similarity_index,
                Colors.RESET,
            ),
            file=sys.stdout,
        )


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj

def compare_jsons(reference_json_path, test_json_path):
    """
    Compare two json files from svtopo and error if they don't match
    """
    with open(reference_json_path, "r") as ref_fh:
        ref_json = json.load(ref_fh)
    with open(test_json_path, "r") as test_fh:
        try:
            test_json = json.load(test_fh)
        except json.decoder.JSONDecodeError:
            print(
                "{}Test json {} is empty or malformed {}".format(
                    Colors.RED, test_file_path, Colors.RESET
                ),
                file=sys.stderr,
            )
            sys.exit(-1)
    if not ordered(ref_json) == ordered(test_json):
        with open("tmp1", "w") as fh:
            pprint(ordered(ref_json), stream=fh)
        with open("tmp2", "w") as fh:
            pprint(ordered(test_json), stream=fh)
        print(
            "{}Test json {} does not match reference json {}{}".format(
                Colors.RED, test_file_path, reference_json_path, Colors.RESET
            ),
            file=sys.stderr,
        )
        sys.exit(-1)
    else:
        print(
            "{}Test json {} matches reference json {}{}".format(
                Colors.GREEN, test_file_path, reference_json_path, Colors.RESET
            ),
            file=sys.stdout,
        )


def compare_text_files(reference_text_path, test_text_path):
    """
    Compare two text files from svtopo and error if they don't match
    """
    with open(reference_text_path, "r") as ref_fh:
        ref_text = ref_fh.read().strip()
    with open(test_text_path, "r") as test_fh:
        test_text = test_fh.read().strip()
    
    if not ref_text == test_text:
        print(
            "{}Test text file {} does not match reference text file {}{}".format(
                Colors.RED, test_text_path, reference_text_path, Colors.RESET
            ),
            file=sys.stderr,
        )
        sys.exit(-1)    
    else:
        print(
            "{}Test text file {} matches reference text file {}{}".format(
                Colors.GREEN, test_text_path, reference_text_path, Colors.RESET
            ),
            file=sys.stdout,
        )


args = setup()
ref_file_path = args.reference_file_path
test_file_path = args.test_file_path
if ref_file_path.endswith("json"):
    compare_jsons(reference_json_path=ref_file_path, test_json_path=test_file_path)
elif ref_file_path.endswith("png"):
    compare_plots(reference_image_path=ref_file_path, test_image_path=test_file_path)
else:
    compare_text_files(reference_text_path=ref_file_path, test_text_path=test_file_path)
