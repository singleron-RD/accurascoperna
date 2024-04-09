#!usr/bin/env python3
"""
Convert 5p barcode to 3p barcode
Generate starsolo command for 3p and 3p5p
Input:
    protocol.json: whitelist and linker
Output:
    3p R1 and R2 fastq
    3p5p R1 and R2 fastq
"""
import argparse
import json
import os

from utils import parse_pattern


def get_protocol_dict(assets_dir):
    """
    Return:
    protocol_dict. Key: protocol name, value: values from protocol.json

    >>> protocol_dict = get_protocol_dict("./assets/")
    >>> protocol_dict["AccuraSCOPE-V1"]["pattern_dict_5p"]
    {'C': [slice(0, 6, None)], 'U': [slice(6, 11, None)]}
    """
    json_file = os.path.join(assets_dir, "protocols.json")
    protocol_dict = json.load(open(json_file))
    whitelist_dir = os.path.join(assets_dir, "whitelist")
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        for x in ["bc_3p", "bc_5p", "linker_3p"]:
            cur[x] = os.path.join(whitelist_dir, protocol, x)
        cur["pattern_dict_3p"] = parse_pattern(cur["pattern_3p"])
        cur["pattern_dict_5p"] = parse_pattern(cur["pattern_5p"])
    return protocol_dict

class Convert:
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--genomeDir', required=True)
    parser.add_argument('--fq1', required=True)
    parser.add_argument('--fq2', required=True)
    parser.add_argument('--assets_dir', required=True)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--thread', required=True)
    parser.add_argument('--ext_args_3p', )
    parser.add_argument('--ext_args_3p5p',)

    args = parser.parse_args()

    runner = Convert(args)
    runner.write_cmd()

