from test.helper_funcs import RESOURCE_DIR

# read file
# parse uniprot accessions (format: only accession between pipes)
# build arrays:
#         first array: 4 columns:
#             - 0: int, 'RefID': matched reference ID in 'rl' dict
#             - 1: int, 'MatchPosStart': start position of match in reference sequence
#             - 2: int, 'MatchScore': score of the match
#             - 3: int, 'ReadID': ID of the read
#         second array: 2 columns
#             - 0: int, dudes internal id for the reference accession
#             - 1: int, LN attribute of reference
