def cal_overall(eva_list: list) -> dict:
    out_dict = dict()
    total_event, total_se = 0, 0
    total_ghd_event, total_ghd = 0, 0
    for i in eva_list:
        total_event += i['total_phase']
        total_se += i['total_se']
        total_ghd_event += i['total_ghd']
        total_ghd += i['total_present']
    out_dict["se_rate"] = total_se / total_event
    out_dict["ghd_percentage"] = total_ghd / total_ghd_event

    return out_dict


def write_results(eva_list: list, chrom: list, NG: tuple, prefix: str, sample_name: str, verbose: bool):
    # write per chromosome results
    with open(prefix + ".perchrom.tsv", 'w') as f:
        f.write("Chromosome\tTotal phased\tPhased block length median\tTotal switch error\tSwith error rate\tghd distance\ttotal ghd\tghd percentage\n")
        for c in range(len(eva_list)):
            f.write(f"chr{chrom[c]}\t\
                    {eva_list[c]['total_phase']}\t\
                    {eva_list[c]['length_median']}\t\
                    {eva_list[c]['total_se']}\t\
                    {eva_list[c]['total_se'] / eva_list[c]['total_phase']}\t\
                    {eva_list[c]['total_present']}\t\
                    {eva_list[c]['total_ghd']}\t\
                    {eva_list[c]['total_present'] / eva_list[c]['total_ghd']}\n")
    
    out_dict = cal_overall(eva_list)
    # write overall results
    with open(prefix + ".overall.tsv", 'w') as g:
        g.write("Sample name\tNG50\tNG90\tSwitch Error rate\tGHD percentage\n")
        g.write(f"{sample_name}\t{NG[0]}\t{NG[1]}\t{out_dict['se_rate']}\t{out_dict['ghd_percentage']}\n")

    if verbose:
        print("Sample name\tNG50\tNG90\tSwitch Error rate\tGHD percentage\n")
        print(f"{sample_name}\t{NG[0]}\t{NG[1]}\t{out_dict['se_rate']}\t{out_dict['ghd_percentage']}\n")
