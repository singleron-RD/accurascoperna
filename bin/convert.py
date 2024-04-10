#!/usr/bin/env python3
"""
Convert p5 barcode to p3 barcode
Generate starsolo command for p3 and p3p5
Input:
    protocol.json: whitelist and linker
Output:
    p3 R1 and R2 fastq
    p5 R1 and R2 fastq
    p3 starsolo command
    p3p5 starsolo command
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
    >>> protocol_dict["AccuraSCOPE-V1"]["pattern_dict_p5"]
    {'C': [slice(0, 6, None)], 'U': [slice(6, 11, None)]}
    """
    json_file = os.path.join(assets_dir, "protocols.json")
    protocol_dict = json.load(open(json_file))
    whitelist_dir = os.path.join(assets_dir, "whitelist")
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        for x in ["bc_p3", "bc_p5", "linker_p3"]:
            cur[x] = os.path.join(whitelist_dir, protocol, cur[x])
        cur["pattern_dict_p3"] = utils.parse_pattern(cur["pattern_p3"])
        cur["pattern_dict_p5"] = utils.parse_pattern(cur["pattern_p5"])
    return protocol_dict


class Starsolo:
    def __init__(self, args):
        self.args = args
        fq1_list = args.fq1.split(",")
        fq2_list = args.fq2.split(",")
        fq1_number = len(fq1_list)
        fq2_number = len(fq2_list)
        if fq1_number != fq2_number:
            sys.exit('fastq1 and fastq2 do not have same file number!')

        self.read_command = 'cat'

        protocol = args.protocol.strip()
        protocol_dict = get_protocol_dict(args.assets_dir)[protocol]
        pattern = protocol_dict["pattern_p3"]
        whitelist_str = protocol_dict["bc_p3"]

        pattern_args = Starsolo.get_solo_pattern(pattern)
        if not whitelist_str:
            whitelist_str = args.whitelist if args.whitelist else 'None'
        self.cb_umi_args = pattern_args + f' --soloCBwhitelist {whitelist_str} '



    @staticmethod
    def get_solo_pattern(pattern) -> str:
        """
        Returns:
            starsolo_cb_umi_args
        """
        pattern_dict = utils.parse_pattern(pattern)
        if len(pattern_dict['U']) != 1 or len(pattern_dict['C']) > 1:
            sys.exit(f'Error: Wrong pattern:{pattern}. \n Solution: fix pattern so that UMI and barcode only have 1 position.\n')
        
        ul = pattern_dict["U"][0].start
        ur = pattern_dict["U"][0].stop
        umi_len = ur - ul

        solo_type = 'CB_UMI_Simple'
        l, r = pattern_dict["C"][0].start, pattern_dict["C"][0].stop
        cb_start = 1
        cb_len = r - l
        umi_start = cb_len + 1
        cb_str = f'--soloCBstart {cb_start} --soloCBlen {cb_len}'
        umi_str = f'--soloUMIstart {umi_start} --soloUMIlen {umi_len}'

        starsolo_cb_umi_args = " ".join([f'--soloType {solo_type} ', cb_str, umi_str, '--soloCBmatchWLtype Exact'])
        return starsolo_cb_umi_args

    def write_cmd(self, ext_args, fq1, fq2, prime):
        """
        If UMI+CB length is not equal to the barcode read length, specify barcode read length with --soloBarcodeReadLength.
        To avoid checking of barcode read length, specify soloBarcodeReadLength 0
        """
        prefix = f'{self.args.sample}.{prime}.'
        fn = f'{self.args.sample}.{prime}.starsolo_cmd.txt'
        cmd = (
            'STAR \\\n'
            f'{self.cb_umi_args} \\\n'
            f'--genomeDir {self.args.genomeDir} \\\n'
            f'--readFilesIn {fq2} {fq1} \\\n'
            f'--readFilesCommand {self.read_command} \\\n'
            f'--outFileNamePrefix {prefix} \\\n'
            f'--runThreadN {self.args.thread} \\\n'
            f'{self.args.common_args} \\\n'
            f'{ext_args} \n'
        )
        logger.info(cmd)
        with open(fn, 'w') as f:
            f.write(cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--genomeDir', required=True)
    parser.add_argument('--fq1', required=True)
    parser.add_argument('--fq2', required=True)
    parser.add_argument('--assets_dir', required=True)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--thread', required=True)
    parser.add_argument('--common_args', required=True)
    parser.add_argument('--p3_args', required=True)
    parser.add_argument('--p3p5_args', required=True)

    args = parser.parse_args()
    protocol_dict = get_protocol_dict(args.assets_dir)
    protocol = protocol_dict[args.protocol]
    bc_p3 = utils.read_one_col(protocol["bc_p3"])
    bc_p5 = utils.read_one_col(protocol["bc_p5"])
    bc_p3_mismatch_dict = utils.get_mismatch_dict(bc_p3)
    bc_p5_mismatch_dict = utils.get_mismatch_dict(bc_p5)
    # p3 p5 bc correspondance
    bc_p5_p3_dict = {p5:p3 for p5,p3 in zip(bc_p5, bc_p3)}

    pattern_dict_p3 = protocol["pattern_dict_p3"]
    pattern_dict_p5 = protocol["pattern_dict_p5"]
    out_fq_fn = {f"{prime}_{r}":f"{args.sample}_{prime}_{r}.fastq" for prime in ['p3','p5'] for r in ['R1','R2']}
    # starsolo cmd
    runner = Starsolo(args)
    p3_fq1,p3_fq2,p5_fq1,p5_fq2 = out_fq_fn['p3_R1'],out_fq_fn['p3_R2'],out_fq_fn['p5_R1'],out_fq_fn['p5_R2']
    runner.write_cmd(args.p3_args, p3_fq1, p3_fq2, 'p3')
    runner.write_cmd(args.p3p5_args, ",".join([p3_fq1,p5_fq1]), ",".join([p3_fq2,p5_fq2]), 'p3p5')

    # add extra bp to force p5 UMI the same length as p3 UMI, p5 bc qual length is the same as p3 bc
    umi_p3_len = pattern_dict_p3["U"][0].stop - pattern_dict_p3["U"][0].start
    umi_p5_len = pattern_dict_p5["U"][0].stop - pattern_dict_p5["U"][0].start
    extra_len =  umi_p3_len - umi_p5_len
    extra_bp = ('ATCG' * 5)[:extra_len]
    bc_p3_len = pattern_dict_p3["C"][0].stop - pattern_dict_p3["C"][0].start
    bc_p5_len = pattern_dict_p5["C"][0].stop - pattern_dict_p5["C"][0].start
    bc_extra_len = bc_p3_len - bc_p5_len
    logger.info(f"UMI extra_len: {extra_len}, UMI extra_bp: {extra_bp}")

    # open for write
    outdict = {k:open(v,'w') for k,v in out_fq_fn.items()}

    fq1_list = args.fq1.split(',')
    fq2_list = args.fq2.split(',')
    for fq1,fq2 in zip(fq1_list, fq2_list):
        fq1 = pyfastx.Fastx(fq1)
        fq2 = pyfastx.Fastx(fq2)
        for (name1, seq1, qual1), (name2,seq2,qual2) in zip(fq1, fq2):
            prime = 'invalid'
            bc_p3 = utils.get_seq_str(seq1, pattern_dict_p3["C"])
            if bc_p3 in bc_p3_mismatch_dict:
                prime = 'p3'
                bc = bc_p3_mismatch_dict[bc_p3]
                umi = utils.get_seq_str(seq1, pattern_dict_p3["U"])
                bc_qual = utils.get_seq_str(qual1, pattern_dict_p3["C"])
                umi_qual = utils.get_seq_str(qual1, pattern_dict_p3["U"])
            else:
                bc_p5 = utils.get_seq_str(seq1, pattern_dict_p5["C"])
                if bc_p5 in bc_p5_mismatch_dict:
                    prime = 'p5'
                    bc = bc_p5_p3_dict[bc_p5_mismatch_dict[bc_p5]]
                    umi = utils.get_seq_str(seq2, pattern_dict_p5["U"]) + extra_bp
                    bc_qual = utils.get_seq_str(qual1, pattern_dict_p5["C"]) + 'F' * bc_extra_len
                    umi_qual = utils.get_seq_str(qual1, pattern_dict_p5["U"]) + 'F' * extra_len
            if prime in ['p3','p5']:
                seq1 = bc+umi
                qual1 = bc_qual+umi_qual
                outdict[f"{prime}_R1"].write(utils.str_fq(name1, seq1, qual1))
                outdict[f"{prime}_R2"].write(utils.str_fq(name2, seq2, qual2))











