import os

def cal_overall(eva_list: list) -> dict:
    out_dict = dict()
    total_snv, total_indel, total_sv = 0, 0, 0
    total_block = 0
    total_event, total_se = 0, 0
    total_ghd_event, total_ghd = 0, 0
    for i in eva_list:
        total_snv += i['total_snv']
        total_indel += i['total_indel']
        total_sv += i['total_sv']
        total_block += i['total_block']
        total_event += i['total_phase']
        total_se += i['total_se']
        total_ghd_event += i['total_ghd']
        total_ghd += i['total_present']
    
    out_dict["total_phase"] = total_event
    out_dict["overall_snv"] = total_snv
    out_dict["overall_indel"] = total_indel
    out_dict["overall_sv"] = total_sv
    out_dict["overall_block_count"] = total_block
    out_dict["se_count"] = total_se
    out_dict["se_rate"] = total_se / total_event
    out_dict["total_ghd"] = total_ghd
    out_dict["total_ghd_event"] = total_ghd_event
    out_dict["ghd_percentage"] = total_ghd / total_ghd_event

    return out_dict


def write_results(eva_list: list, chrom: list, NG: tuple, prefix: str, sample_name: str, verbose: bool):
    # write per chromosome results
    with open(os.path.abspath(prefix + ".perchrom.tsv"), 'w') as f:
        f.write("Chromosome\t\
                Total Phased\t\
                Total SNV\t\
                Total INDEL\t\
                Total SV\t\
                Block count\t\
                Phased block length median\t\
                Total switch error\t\
                Switch error rate\t\
                Total GHD\t\
                Total GHD event\t\
                GHD percentage\n")
        for c in range(len(eva_list)):
            f.write(f"chr{chrom[c]}\t\
                    {eva_list[c]['total_phase']}\t\
                    {eva_list[c]['total_snv']}\t\
                    {eva_list[c]['total_indel']}\t\
                    {eva_list[c]['total_sv']}\t\
                    {eva_list[c]['total_block']}\t\
                    {eva_list[c]['length_median']}\t\
                    {eva_list[c]['total_se']}\t\
                    {eva_list[c]['total_se'] / eva_list[c]['total_phase']}\t\
                    {eva_list[c]['total_present']}\t\
                    {eva_list[c]['total_ghd']}\t\
                    {eva_list[c]['total_present'] / eva_list[c]['total_ghd']}\n")
    
    out_dict = cal_overall(eva_list)
    # write overall results
    with open(os.path.abspath(prefix + ".overall.tsv"), 'w') as g:
        g.write(f"Sample name\t\
                Total Phased\t\
                Total SNV\t\
                Total INDEL\t\
                Total SV\t\
                Total block\t\
                NG50\t\
                NG90\t\
                Switch Error count\t\
                Switch Error rate\t\
                Total GHD\t\
                Total GHD event\t\
                GHD percentage\n")
        g.write(f"{sample_name}\t\
                {out_dict['total_phase']}\t\
                {out_dict['overall_snv']}\t\
                {out_dict['overall_indel']}\t\
                {out_dict['overall_sv']}\t\
                {out_dict['overall_block_count']}\t\
                {NG[0]}\t\
                {NG[1]}\t\
                {out_dict['se_count']}\t\
                {out_dict['se_rate']}\t\
                {out_dict['total_ghd']}\t\
                {out_dict['total_ghd_event']}\t\
                {out_dict['ghd_percentage']}\n")
    
    if verbose:
        print(f"\t\tSample name: {sample_name}\r")
        print(f"\t\tTotal phased: {out_dict['total_phase']}\r")
        print(f"\t\tTotal SNV: {out_dict['overall_snv']}\r")
        print(f"\t\tTotal INDEL: {out_dict['overall_indel']}\r")
        print(f"\t\tTotal SV: {out_dict['overall_sv']}\r")
        print(f"\t\tTotal block: {out_dict['overall_block_count']}\r")
        print(f"\t\tNG50: {NG[0]}\r")
        print(f"\t\tNG90: {NG[1]}\r")
        print(f"\t\tSwitch Error rate: {out_dict['se_rate']}\r")
        print(f"\t\tGHD percentage: {out_dict['ghd_percentage']}\r")

