#!/usr/bin/env python3
"""
Convert 5p barcode to 3p barcode
Generate starsolo command for 3p and 3p5p
Input:
    protocol.json: whitelist and linker
Output:
    3p R1 and R2 fastq
    5p R1 and R2 fastq
    3p starsolo command
    3p5p starsolo command
"""
import argparse
import json
import os
import logging

import pyfastx
import utils
import sys

# stdout
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO, stream=sys.stdout)
logger = logging.getLogger(__name__)

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
            cur[x] = os.path.join(whitelist_dir, protocol, cur[x])
        cur["pattern_dict_3p"] = utils.parse_pattern(cur["pattern_3p"])
        cur["pattern_dict_5p"] = utils.parse_pattern(cur["pattern_5p"])
    return protocol_dict


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
    protocol_dict = get_protocol_dict(args.assets_dir)
    protocol = protocol_dict[args.protocol]
    bc_3p = utils.read_one_col(protocol["bc_3p"])
    bc_5p = utils.read_one_col(protocol["bc_5p"])
    bc_3p_mismatch_dict = utils.get_mismatch_dict(bc_3p)
    bc_5p_mismatch_dict = utils.get_mismatch_dict(bc_5p)
    # 3p 5p bc correspondance
    bc_5p_3p_dict = {p5:p3 for p5,p3 in zip(bc_5p, bc_3p)}

    pattern_dict_3p = protocol["pattern_dict_3p"]
    pattern_dict_5p = protocol["pattern_dict_5p"]
    outdict = {f"{prime}_{r}":f"{args.sample}_{prime}_{r}.fastq" for prime in ['3p','5p'] for r in ['R1','R2']}

    # add extra bp to force 5p UMI the same length as 3p UMI
    umi_3p_len = pattern_dict_3p["U"][0].stop - pattern_dict_3p["U"][0].start
    umi_5p_len = pattern_dict_5p["U"][0].stop - pattern_dict_5p["U"][0].start
    extra_len =  umi_3p_len - umi_5p_len
    extra_bp = ('ATCG' * 5)[:extra_len]
    logger.info(f"extra_len: {extra_len}, extra_bp: {extra_bp}")

    # open for write
    outdict = {k:open(v,'w') for k,v in outdict.items()}

    fq1_list = args.fq1.split(',')
    fq2_list = args.fq1.split(',')
    for fq1,fq2 in zip(fq1_list, fq2_list):
        fq1 = pyfastx.Fastx(fq1)
        fq2 = pyfastx.Fastx(fq2)
        for (name1, seq1, qual1), (name2,seq2,qual2) in zip(fq1, fq2):
            prime = 'invalid'
            bc_3p = utils.get_seq_str(seq1, pattern_dict_3p["C"])
            if bc_3p in bc_3p_mismatch_dict:
                prime = '3p'
                bc = bc_3p_mismatch_dict[bc_3p]
                umi = utils.get_seq_str(seq1, pattern_dict_3p["U"])
                bc_qual = utils.get_seq_str(qual1, pattern_dict_3p["C"])
                umi_qual = utils.get_seq_str(qual1, pattern_dict_3p["U"])
            else:
                bc_5p = utils.get_seq_str(seq1, pattern_dict_5p["C"])
                if bc_5p in bc_5p_mismatch_dict:
                    prime = '5p'
                    bc = bc_5p_3p_dict[bc_5p_mismatch_dict[bc_5p]]
                    umi = utils.get_seq_str(seq2, pattern_dict_5p["U"]) + extra_bp
                    bc_qual = utils.get_seq_str(qual1, pattern_dict_5p["C"])
                    umi_qual = utils.get_seq_str(qual1, pattern_dict_5p["U"]) + 'F' * extra_len
            if prime in ['3p','5p']:
                seq1 = bc+umi
                qual1 = bc_qual+umi_qual
                outdict[f"{prime}_R1"].write(utils.str_fq(name1, seq1, qual1))
                outdict[f"{prime}_R2"].write(utils.str_fq(name2, seq2, qual2))











