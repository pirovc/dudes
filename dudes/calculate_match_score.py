import re
from typing import Union


def calculate_match_score(
    cigar_string: str,
    edit_distance: Union[None, int] = None,
    mismatches: Union[None, int] = None,
) -> int:
    """Calculate the match score of the read.

    The match score is calculated as follows:
        - (matches + insertions) - (Edit distance to the reference + deletions + insertions)
        -

    Args:
        cigar_string: matches, insertions and deletions are read from the CIGAR string
        edit_distance: Edit distance to the reference is defined as: Number of differences (mismatches
          plus inserted and deleted bases) between the sequence and reference
        mismatches: Number of mismatches

    Returns: int, match score

    """
    cig = {"I": 0, "D": 0, "M": 0}
    for val, ci in re.findall(r"(\d+)([IDM])", cigar_string):
        cig[ci] += int(val)
    if (edit_distance is None) and (mismatches is None):
        raise ValueError("One of the parameters edit_distance or mismatches is required.")
    if edit_distance is None:
        edit_distance = mismatches + cig["I"] + cig["D"]
    return (cig["M"] + cig["I"]) - (edit_distance + cig["D"] + cig["I"])
